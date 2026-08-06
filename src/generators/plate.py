import numpy as np
from scipy.fft import dctn, dstn, idctn, idstn

from src.utils.damping import calc_damping
from src.utils.excitation import calc_exc
from src.utils.math import pow2floor
from src.utils.sav import calc_g_mod, calc_g_std, calc_psi


class Plate:
    def __init__(
            self,
            fs: int,
            dur: float,
            kappa: float,
            sigma0: float,
            sigma1: float,
            xe: float,
            ye: float,
            xo: float,
            yo: float,
            ratio: float,
            exc_factor: float = 1.0,
            out_factor: float = 1.0
        ):
        self._fs = int(fs)
        self._dur = dur
        self._kappa = kappa
        self._sigma0 = sigma0
        self._sigma1 = sigma1
        self._xe = xe
        self._ye = ye
        self._xo = xo
        self._yo = yo
        self._ratio = ratio
        self._exc_factor = exc_factor
        self._out_factor = out_factor

        self._k = 1.0 / self._fs
        self._num_samples = int(np.floor(self._dur * self._fs))
        self._energy = None
        self._power = None

    @classmethod
    def from_physical_parameters(
            cls,
            fs: int,
            dur: float,
            Lx: float,
            Ly: float,
            rho: float,
            H: float,
            E: float,
            nu: float,
            loss: list[list[float, float]],
            exc_pos: list[float, float],
            out_pos: list[float, float],
            **kwargs
        ):
        # derived parameters
        L = np.sqrt(Lx * Ly)
        D = E * H**3 / (12.0 * (1.0 - nu**2))

        # parameters of a scaled PDE
        kappa = np.sqrt(D / (rho * H)) / L**2
        xe = exc_pos[0] / L
        ye = exc_pos[1] / L
        xo = out_pos[0] / L
        yo = out_pos[1] / L
        ratio = Lx / Ly
        out_factor = H / np.sqrt(6.0 * (1.0 - nu**2))
        exc_factor = 1.0 / (rho * H) / L**2 / out_factor

        # damping parameters
        sigma0, sigma1 = calc_damping(loss, kappa)

        return cls(
            fs,
            dur,
            kappa,
            sigma0,
            sigma1,
            xe,
            ye,
            xo,
            yo,
            ratio,
            exc_factor=exc_factor,
            out_factor=out_factor,
            **kwargs
        )

    def calc_energy_total(self, st: float = 0.0) -> np.ndarray[tuple[int], float] | None:
        if (self._energy is not None) and (self._power is not None):
            st_idx = int(np.ceil(st * self._fs))
            energy_total = np.zeros(((self._num_samples - st_idx),))
            energy_total[0] = self._energy[st_idx]
            energy_total[1:] = self._energy[st_idx + 1:] - self._k * np.cumsum(self._power[st_idx:-1])
        else:
            energy_total = None
        return energy_total

    def __call__(
            self,
            exc_amp: float,
            exc_dur: float,
            exc_st: float,
            exc_type: int,
            num_repeat: int = 1
        ) -> np.ndarray[tuple[int], float]:
        raise NotImplementedError

    @property
    def k(self):
        return self._k

    @ property
    def num_samples(self):
        return self._num_samples

    @property
    def energy(self):
        return self._energy

    @property
    def power(self):
        return self._power


class LinearPlateModal(Plate):
    def __init__(
            self,
            fs: int,
            dur: float,
            kappa: float,
            sigma0: float,
            sigma1: float,
            xe: float,
            ye: float,
            xo: float,
            yo: float,
            ratio: float,
            exc_factor: float = 1.0,
            out_factor: float = 1.0,
            fmax: float | None = None,
            use_exact: bool = True
        ):
        super().__init__(fs, dur, kappa, sigma0, sigma1, xe, ye, xo, yo, ratio, exc_factor, out_factor)

        # calculate maximum frequency and wavenumber from stability condition
        if fmax is None:
            if use_exact:
                omega_max = np.pi * self._fs
            else:
                omega_max = 2.0 * self._fs
        else:
            omega_max = 2.0 * np.pi * fmax
        beta_max = np.sqrt(omega_max / self._kappa)
        self._num_modes_x = int(np.floor(
            np.sqrt((beta_max**2 / np.pi**2 - self._ratio) * self._ratio)
        ))
        self._num_modes_y = int(np.floor(
            np.sqrt((beta_max**2 / np.pi**2 - 1.0 / self._ratio) / self._ratio)
        ))

        # calculate wavenumbers
        self._mx, self._my = np.meshgrid(
            np.arange(start=1, stop=(self._num_modes_x + 1)),
            np.arange(start=1, stop=(self._num_modes_y + 1)),
            indexing="ij"
        )
        self._mx = self._mx.flatten(order="C")
        self._my = self._my.flatten(order="C")
        self._beta_x = self._mx * np.pi / np.sqrt(self._ratio)
        self._beta_y = self._my * np.pi * np.sqrt(self._ratio)
        self._beta = np.sqrt(self._beta_x**2 + self._beta_y**2)

        # discard unstable modes
        self._mx = self._mx[self._beta < beta_max]
        self._my = self._my[self._beta < beta_max]
        self._beta_x = self._beta_x[self._beta < beta_max]
        self._beta_y = self._beta_y[self._beta < beta_max]
        self._beta = self._beta[self._beta < beta_max]
        self._num_modes = self._beta.shape[0]

        # calculate modes
        omega = self._kappa * self._beta**2
        self._omega_sq = omega**2
        self._sigma = self._sigma0 + self._sigma1 * self._beta**2

        # modify modes to mitigate numerical dispersion
        if use_exact:
            exp_s = np.exp(-self._sigma * self._k)
            exp_2s = np.exp(-2.0 * self._sigma * self._k)
            cos_d = np.cos(np.sqrt(self._omega_sq - self._sigma**2) * self._k)
            self._omega_sq = 2.0 / self._k**2 * (1.0 - 2.0 * exp_s * cos_d + exp_2s) / (1.0 + exp_2s)
            self._sigma = 1.0 / self._k * (1.0 - exp_2s) / (1.0 + exp_2s)

        # precompute update matrices (q[n] = B*q[n-1] + C*q[n-2] + J*f[n])
        self._d = 1.0 + self._k * self._sigma
        self._B = (2.0 - self._k**2 * self._omega_sq) / self._d
        self._C = (self._k * self._sigma - 1.0) / self._d
        self._exc_basis = self._calc_basis(self._xe, self._ye)
        self._Je = (self._k**2 * self._exc_basis) / self._d
        self._Jo = self._calc_basis(self._xo, self._yo)

    def _calc_basis(self, x: np.ndarray[tuple[int], float], y: np.ndarray[tuple[int], float]) -> np.ndarray[tuple[int], float]:
        return 2.0 * np.sin(self._beta_x * x) * np.sin(self._beta_y * y)

    def _calc_energy(self, q: np.ndarray[tuple[int], float], q1: np.ndarray[tuple[int], float]) -> float:
        q_tm = (q - q1) / self._k
        return 0.5 * (
            q_tm.dot(q_tm)
            + q.dot(self._omega_sq * q1)
        )

    def _calc_power(self, q: np.ndarray[tuple[int], float], q2: np.ndarray[tuple[int], float], fe_point: float) -> float:
        q_t = (q - q2) / (2.0 * self._k)
        diss_power = np.negative(q_t.dot(2.0 * self._sigma * q_t))
        exc_power = q_t.dot(self._exc_basis) * fe_point
        return (diss_power + exc_power)

    def __call__(
            self,
            exc_amp: float,
            exc_dur: float,
            exc_st: float,
            exc_type: int,
            num_repeat: int = 1,
            init_amp: float = 0.0
        ) -> np.ndarray[tuple[int], float]:
        # calculate excitation
        fe_points = calc_exc(self._fs, self._num_samples, exc_amp, exc_dur, exc_st, exc_type, num_repeat)
        fe_points = self._exc_factor * fe_points
        self._num_samples = num_repeat * self._num_samples

        # initialise
        q2 = np.zeros((self._num_modes,), dtype=np.float64)
        q1 = np.zeros((self._num_modes,), dtype=np.float64)
        q2[0] = init_amp
        q1[0] = init_amp

        # allocate output arrays
        out = np.zeros((self._num_samples,), dtype=np.float64)
        self._energy = np.zeros((self._num_samples,), dtype=np.float64)
        self._power = np.zeros((self._num_samples,), dtype=np.float64)

        # main loop
        for n in range(self._num_samples):
            # update state
            q = self._B * q1 + self._C * q2 + self._Je * fe_points[n]
            out[n] = self._Jo.dot(q)

            # calculate energy
            self._energy[n] = self._calc_energy(q1, q2)
            self._power[n] = self._calc_power(q, q2, fe_points[n])

            # shift state
            q2 = q1.copy()
            q1 = q.copy()
        out = self._out_factor * out

        return out

    @property
    def num_modes(self):
        return self._num_modes


class VKPlatePSTD(LinearPlateModal):
    def __init__(
            self,
            fs: int,
            dur: float,
            kappa: float,
            sigma0: float,
            sigma1: float,
            xe: float,
            ye: float,
            xo: float,
            yo: float,
            ratio: float,
            exc_factor: float = 1.0,
            out_factor: float = 1.0,
            fmax: float | None = None,
            use_exact: bool = True,
            lambda0: float = 1e3,
            norm_type: int = 1,
            os_factor: float = 1.5
        ):
        super().__init__(fs, dur, kappa, sigma0, sigma1, xe, ye, xo, yo, ratio, exc_factor, out_factor, fmax, use_exact)
        self._lambda0 = lambda0
        self._norm_type = norm_type
        self._os_factor = os_factor

        # calculate 1D modes for the Airy function
        mx_1d = np.concatenate((
            np.arange(start=1, stop=(self._num_modes_x + 1)),
            np.zeros((self._num_modes_y,), dtype=np.int64)
        ))
        my_1d = np.concatenate((
            np.zeros((self._num_modes_x,), dtype=np.int64),
            np.arange(start=1, stop=(self._num_modes_y + 1))
        ))
        beta_x_1d = mx_1d * np.pi / np.sqrt(self._ratio)
        beta_y_1d = my_1d * np.pi * np.sqrt(self._ratio)
        beta_1d = np.sqrt(beta_x_1d**2 + beta_y_1d**2)

        # append 1D modes to existing 2D modes
        self._mx_ext = np.append(self._mx, mx_1d)
        self._my_ext = np.append(self._my, my_1d)
        self._beta_x_ext = np.append(self._beta_x, beta_x_1d)
        self._beta_y_ext = np.append(self._beta_y, beta_y_1d)
        self._beta_ext = np.append(self._beta, beta_1d)
        self._num_modes_ext = self._beta_ext.shape[0]

        # prepare transforms
        num_nodes_x = int(np.ceil(self._os_factor * self._num_modes_x)) + 1
        num_nodes_y = int(np.ceil(self._os_factor * self._num_modes_y)) + 1
        self._buf_s = np.zeros((num_nodes_x, num_nodes_y))
        self._buf_c = np.zeros((num_nodes_x, num_nodes_y))
        self._buf_c_ext = np.zeros((num_nodes_x, num_nodes_y))
        self._tf_factor = np.sqrt(num_nodes_x) * np.sqrt(num_nodes_y)

        # precompute constants
        self._nl_factor = self._k * self._kappa

    def _calc_energy(self, q: np.ndarray[tuple[int], float], q1: np.ndarray[tuple[int], float], psi: float) -> float:
        energy_lin = super()._calc_energy(q, q1)
        return (
            energy_lin
            + 0.5 * self._kappa**2 * psi**2
        )

    def _idst(self, q: np.ndarray[tuple[int], float]) -> np.ndarray[tuple[int, int], float]:
        self._buf_s[self._mx - 1, self._my - 1] = q
        return self._tf_factor * idstn(self._buf_s, type=2, norm="ortho")

    def _dst(self, U: np.ndarray[tuple[int, int], float]) -> np.ndarray[tuple[int], float]:
        return (1.0 / self._tf_factor) * dstn(U, type=2, norm="ortho")[self._mx - 1, self._my - 1]

    def _idct(self, q: np.ndarray[tuple[int], float]) -> np.ndarray[tuple[int, int], float]:
        self._buf_c[self._mx, self._my] = q
        return self._tf_factor * idctn(self._buf_c, type=2, norm="ortho")

    def _dct(self, U: np.ndarray[tuple[int, int], float]) -> np.ndarray[tuple[int], float]:
        return (1.0 / self._tf_factor) * dctn(U, type=2, norm="ortho")[self._mx, self._my]

    def _idct_ext(self, q: np.ndarray[tuple[int], float]) -> np.ndarray[tuple[int, int], float]:
        self._buf_c_ext[self._mx_ext, self._my_ext] = q
        return self._tf_factor * idctn(self._buf_c_ext, type=2, norm="ortho")

    def _dct_ext(self, U: np.ndarray[tuple[int, int], float]) -> np.ndarray[tuple[int], float]:
        return (1.0 / self._tf_factor) * dctn(U, type=2, norm="ortho")[self._mx_ext, self._my_ext]

    def _calc_potential(
            self,
            q: np.ndarray[tuple[int], float],
            calc_grad: bool = True
        ) -> tuple[float, np.ndarray[tuple[int], float] | None]:
        # calculate derivatives of displacement
        Uxx = np.negative(self._idst(self._beta_x**2 * q))
        Uyy = np.negative(self._idst(self._beta_y**2 * q))
        Uxy = self._idct(self._beta_x * self._beta_y * q)

        # calculate Airy function
        xi = np.negative(
            self._dct_ext(2.0 * (Uxx * Uyy - Uxy * Uxy))
        ) / self._beta_ext**4

        # calculate potential
        xi_d = self._beta_ext**2 * xi
        V = 0.25 * xi_d.dot(xi_d)

        if calc_grad:
            # calculate derivatives of Airy function
            Phi_xx = np.negative(self._idct_ext(self._beta_x_ext**2 * xi))
            Phi_yy = np.negative(self._idct_ext(self._beta_y_ext**2 * xi))
            Phi_xy = self._idst(self._beta_x * self._beta_y * xi[:self._num_modes])

            # calculate gradient
            grad_V = np.negative(
                self._dst(Uxx * Phi_yy + Uyy * Phi_xx - 2.0 * Uxy * Phi_xy)
            )
        else:
            grad_V = None

        return V, grad_V

    def __call__(
            self,
            exc_amp: float,
            exc_dur: float,
            exc_st: float,
            exc_type: int,
            num_repeat: int = 1,
            init_amp: float = 0.0
        ) -> np.ndarray[tuple[int], float]:
        # calculate excitation
        fe_points = calc_exc(self._fs, self._num_samples, exc_amp, exc_dur, exc_st, exc_type, num_repeat)
        fe_points = self._exc_factor * fe_points
        self._num_samples = num_repeat * self._num_samples

        # initialise
        q2 = np.zeros((self._num_modes,), dtype=np.float64)
        q1 = np.zeros((self._num_modes,), dtype=np.float64)
        q2[0] = init_amp
        q1[0] = init_amp
        q_half = 0.5 * (q1 + q2)
        V_half, _ = self._calc_potential(q_half, calc_grad=False)
        psi1 = calc_psi(V_half)

        # allocate output arrays
        out = np.zeros((self._num_samples,), dtype=np.float64)
        self._energy = np.zeros((self._num_samples,), dtype=np.float64)
        self._power = np.zeros((self._num_samples,), dtype=np.float64)

        # allocate arrays for testing
        self._psi_half = np.zeros((self._num_samples,), dtype=np.float64)
        self._V_half = np.zeros((self._num_samples,), dtype=np.float64)
        self._V_diff = np.zeros((self._num_samples,), dtype=np.float64)
        self._u_max = np.zeros((self._num_samples,), dtype=np.float64)

        # main loop
        for n in range(self._num_samples):
            # calculate nonlinearity
            q_half = 0.5 * (q1 + q2)
            p_half = (q1 - q2) / self._k
            V, grad_V = self._calc_potential(q1)
            V_half, _ = self._calc_potential(q_half, calc_grad=False)
            g = (
                calc_g_std(V, grad_V)
                + calc_g_mod(V_half, psi1, p_half, lambda0=self._lambda0, norm_type=self._norm_type)
            )
            a = 0.5 * self._nl_factor * g
            b = a / self._d

            # update state
            q = (
                self._B * q1
                + self._C * q2
                + self._Je * fe_points[n]
                + b * a.dot(q2)
                - self._nl_factor**2 * psi1 * (g / self._d)
            )
            q = q - b * (a.dot(q) / (1.0 + a.dot(b)))
            psi = psi1 + 0.5 * g.dot(q - q2)
            out[n] = self._Jo.dot(q)

            # calculate energy
            self._energy[n] = self._calc_energy(q1, q2, psi1)
            self._power[n] = self._calc_power(q, q2, fe_points[n])

            # save outputs for testing
            self._psi_half[n] = psi1
            self._V_half[n] = V_half
            self._V_diff[n] = (V - 0.25 * grad_V.dot(q1)) / pow2floor(V)
            self._u_max[n] = np.max(np.abs(self._idst(q)))

            # shift state
            q2 = q1.copy()
            q1 = q.copy()
            psi1 = psi
        out = self._out_factor * out

        return out

    @property
    def psi_half(self):
        return self._psi_half

    @property
    def V_half(self):
        return self._V_half

    @property
    def V_diff(self):
        return self._V_diff

    @property
    def u_max(self):
        return self._u_max
