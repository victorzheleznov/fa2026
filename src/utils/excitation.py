import numpy as np


def calc_exc(
        fs: int,
        num_samples: int,
        amp: float,
        dur: float,
        st: float,
        type: int,
        num_repeat: int = 1
    ) -> np.ndarray[tuple[int], float]:
    t_points = np.arange(start=0, stop=num_samples, step=1) / fs
    fe_points = 0.5 * amp * (
        1.0 - np.cos(type * np.pi * (t_points - st) / dur)
    ) * np.logical_and(
        (t_points >= st),
        (t_points <= st + dur)
    ).astype(int)
    fe_points = np.tile(fe_points, num_repeat)
    return fe_points
