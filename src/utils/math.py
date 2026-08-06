import numpy as np


def log2clip(x: np.ndarray[tuple[int], float]) -> np.ndarray[tuple[int], float]:
    return np.log2(np.clip(x, a_min=np.finfo(float).eps, a_max=None))


def pow2floor(x: np.ndarray[tuple[int], float]) -> np.ndarray[tuple[int], float]:
    return 2**np.floor(log2clip(x))
