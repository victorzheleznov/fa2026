from math import ceil, log

import librosa
import matplotlib.pyplot as plt
import numpy as np
from librosa.display import specshow
from scipy.signal.windows import get_window

LOG_TICKS = [62.5, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]
LOG_LABELS = ["62.5", "125", "250", "500", "1k", "2k", "4k", "8k", "16k"]


def plot_spec(
        x: np.ndarray[tuple[int], float],
        fs: int,
        win_name: str = "blackmanharris",
        win_len: int = 4096,
        overlap: float = 0.9,
        dynamic_range: float = 120,
        freq_range: tuple[float, float] = (20, 20e3),
        freq_axis: str = "linear",
        plot_bar: bool = False,
        axs = None
    ):
    if axs is None:
        fig, axs = plt.subplots(nrows=1, ncols=1, layout="constrained")
    else:
        fig = axs.get_figure()

    win = get_window(win_name, win_len, fftbins=True)
    fft_len = 2**(ceil(log(win_len) / log(2)))
    hop_len = round(win_len * (1.0 - overlap))

    stft = librosa.stft(
        x,
        n_fft=fft_len,
        hop_length=hop_len,
        win_length=win_len,
        window=win,
        center=False
    )
    mag = np.abs(stft)
    mag_db = 20.0 * np.log10(mag / np.max(mag))

    num_frames = mag.shape[1]
    t_points = np.arange(start=0, stop=num_frames, step=1) / fs * hop_len
    freqs = np.arange(start=0, stop=((fft_len // 2) + 1), step=1) * fs / fft_len

    img = specshow(
        mag_db,
        x_coords=t_points,
        y_coords=freqs,
        x_axis="time",
        y_axis=freq_axis,
        cmap="magma",
        vmin=-dynamic_range,
        vmax=0,
        ax=axs
    )
    axs.set_ylim(freq_range)
    axs.set_xlabel("Time [sec]")
    axs.set_ylabel("Frequency [Hz]")
    if plot_bar:
        fig.colorbar(img, label="[dB]")
    if freq_axis == "log":
        axs.tick_params(left=False)
        axs.set_yticks(LOG_TICKS, LOG_LABELS)

    return fig, axs


def plot_against_time(x: np.ndarray[tuple[int], float], fs: int, st: float = 0.0, axs = None, **kwargs):
    if axs is None:
        fig, axs = plt.subplots(nrows=1, ncols=1, layout="constrained")
    else:
        fig = axs.get_figure()

    st_idx = int(np.ceil(st * fs))
    t_points = np.arange(start=st_idx, stop=(st_idx + len(x)), step=1) / fs
    axs.plot(t_points, x, **kwargs)
    axs.set_xlabel("Time [sec]")

    return fig, axs


def plot_with_eps(x: np.ndarray[tuple[int], float], fs: int, st: float = 0.0, axs = None, plot_grid: bool = True, **kwargs):
    if axs is None:
        fig, axs = plt.subplots(nrows=1, ncols=1, layout="constrained")
    else:
        fig = axs.get_figure()

    kwargs.update({"linestyle": "", "marker": "o"})
    plot_against_time(x, fs, st=st, axs=axs, **kwargs)

    if plot_grid:
        ylim = axs.get_ylim()
        if np.max(np.abs(np.array(ylim))) < 1e-12:
            eps = np.finfo(float).eps
            eps_ticks = np.append(
                np.arange(start=0, stop=ylim[0], step=-eps),
                np.arange(start=0, stop=ylim[1], step=eps)
            )
            axs.yaxis.remove_overlapping_locs = False
            axs.yaxis.set_tick_params(which="minor", length=0)
            axs.set_yticks(eps_ticks, minor=True)
            axs.grid(which="minor", axis="y", visible=True)
            axs.grid(which="minor", axis="x", visible=False)
            axs.grid(which="major", axis="both", visible=False)
            axs.set_ylim(ylim)

    return fig, axs


def num2str(x: float, var: str | None = None, num_digits: int = 0) -> str:
    if x == 0.0:
        s = f"{x:.{num_digits}f}"
    else:
        sign = np.sign(x)
        pow = int(np.floor(np.log10(np.abs(x))))
        coef = np.abs(x) / 10**pow
        if np.abs(pow) > 1:
            s = f"10^{{{pow:d}}}"
            if coef != 1.0:
                s = f"{coef:.{num_digits}f} \\times " + s
            if sign < 0:
                s = "-" + s
        else:
            s = f"{x:.{num_digits}f}"

    if var is not None:
        s = var + "=" + s

    s = "$" + s + "$"
    return s
