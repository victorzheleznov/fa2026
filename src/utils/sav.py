import numpy as np

from src.utils.const import EPS


def calc_psi(V: float) -> float:
    return np.sqrt(2.0 * V + EPS)


def calc_g_std(V: float, grad_V: np.ndarray[tuple[int], float]) -> np.ndarray[tuple[int], float]:
    true_psi = calc_psi(V)
    return (grad_V / true_psi)


def calc_g_mod(
        V: float,
        psi: float,
        p: np.ndarray[tuple[int], float],
        lambda0: float = 1e3,
        norm_type: int = 1
    ) -> np.ndarray[tuple[int], float]:
    true_psi = calc_psi(V)
    if norm_type == 1:
        vec = np.sign(p)
    elif norm_type == 2:
        vec = p
    norm = vec.dot(p)
    return -lambda0 * (psi - true_psi) * vec / (norm + EPS)
