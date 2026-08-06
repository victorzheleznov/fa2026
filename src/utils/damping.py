import math


def calc_damping(loss: list[list[float, float]], kappa: float, gamma: float = 0.0) -> tuple[float, float]:
    def beta_sq(f: float) -> float:
        omega = 2.0 * math.pi * f
        return (-gamma**2 + math.sqrt(gamma**4 + 4.0 * kappa**2 * omega**2)) / (2.0 * kappa**2)

    beta_sq_0 = beta_sq(loss[0][0])
    t60_0 = loss[0][1]
    beta_sq_1 = beta_sq(loss[1][0])
    t60_1 = loss[1][1]

    sigma0 = 6.0 * math.log(10.0) / (beta_sq_1 - beta_sq_0) * (beta_sq_1 / t60_0 - beta_sq_0 / t60_1)
    sigma1 = 6.0 * math.log(10.0) / (beta_sq_1 - beta_sq_0) * (-1.0 / t60_0 + 1.0 / t60_1)

    return sigma0, sigma1
