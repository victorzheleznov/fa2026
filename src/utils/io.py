from pathlib import Path

import numpy as np
import soundfile as sf
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox

GOLDEN_RATIO = (5**0.5 - 1) / 2


def write_wav(
        file: str | Path,
        data: np.ndarray[tuple[int, ...], float],
        fs: int,
        normalise: bool = False,
        fade_in: bool = False,
        subtype: str = "DOUBLE"
    ):
    if normalise:
        data /= np.max(np.abs(data))
    if fade_in:
        end_idx = int(np.floor(5e-3 * fs))
        data[:end_idx] *= 0.5 * (1.0 - np.cos(np.pi * np.arange(end_idx) / (end_idx - 1)))
    sf.write(file, data, samplerate=fs, subtype=subtype)


def save_fig(
        dir: str | Path,
        name: str,
        fig: Figure,
        width: float,
        height: float = None,
        dpi: float = 300,
        format: str = "pdf"
    ):
    if height is None:
        height = width * GOLDEN_RATIO
    fig.set_size_inches(width, height)
    fig.savefig(dir / (name + "." + format), format=format, bbox_inches=Bbox([[0, 0], [width, height]]), dpi=dpi)
