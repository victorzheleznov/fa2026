# %%
# print commit hash
import subprocess
print(subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip())

# %%
# check git diff
print(subprocess.check_output(["git", "diff", "HEAD"]).decode("ascii"))

# %%
import sys
sys.path.append("../")

# %%
import copy
from pathlib import Path

import numpy as np

from src.generators import VKPlatePSTD
from src.utils.plot import plot_spec, plot_with_eps, plot_against_time, num2str
from src.utils.sav import calc_psi
from src.utils.const import AUDIO_RATE
from src.utils.math import pow2floor
from src.utils.io import save_fig, write_wav

# %%
%matplotlib ipympl
import matplotlib as mpl
import matplotlib.pyplot as plt
import IPython.display as ipd
import scienceplots
from tabulate import tabulate

# %%
FA2026_TEXT_WIDTH = 6.69
FA2026_COLUMN_WIDTH = 3.18
FA2026_TEXT_HEIGHT = 9.52

# %%
# default plate parameters
plate_kwargs = {
    "fs": AUDIO_RATE,
    "dur": 3,
    "kappa": 8,
    "sigma0": 1.3,
    "sigma1": 1e-4,
    "xe": 0.17,
    "ye": 0.42,
    "xo": 0.64,
    "yo": 0.79,
    "ratio": 1.1,
    "fmax": 17e3,
    "lambda0": 1e3
}

# default excitation parameters
exc_kwargs = {
    "exc_amp": 2e6,
    "exc_dur": 2e-3,
    "exc_st": 0,
    "exc_type": 2,
    "num_repeat": 1,
    "init_amp": 0
}

# %%
# output directory
out_dir = "../out/fa2026"
out_dir = Path(out_dir)
out_dir.mkdir(exist_ok=True, parents=True)

# %%
# define testing function
def test_plate(
        plate_kwargs: dict[str, float],
        exc_kwargs: dict[str, float],
        freq_range: tuple[float, float] = (0, AUDIO_RATE // 2)
    ):
    print(tabulate(plate_kwargs.items(), floatfmt="g", numalign="left", headers=["Plate", ""]), end="\n\n")
    print(tabulate(exc_kwargs.items(), floatfmt="g", numalign="left", headers=["Excitation", ""]), end="\n\n")

    plate = VKPlatePSTD(**plate_kwargs)
    out = plate(**exc_kwargs)

    print(f"Number of modes: {plate.num_modes}")
    print(f"Max. abs. displacement: {np.max(plate.u_max):g}")

    ipd.display(ipd.Audio(out, rate=AUDIO_RATE))

    plt.style.use(["default"])
    fig, _ = plot_spec(out, AUDIO_RATE, freq_range=freq_range)
    plt.show()

    return plate, out, fig

# %% [markdown]
# ## Energy Balance

# %%
# simulate plate
exc_kwargs_init_amp = copy.deepcopy(exc_kwargs)
exc_kwargs_init_amp.update({"exc_amp": 0, "init_amp": 5})
plate, out, fig = test_plate(plate_kwargs, exc_kwargs_init_amp, freq_range=(0, 5e3))
name = "energy"
write_wav(out_dir / (name + ".wav"), out, AUDIO_RATE, normalise=True, fade_in=True, subtype="PCM_24")
save_fig(out_dir, name, fig, width=(1920 / 300), height=(1440 / 300), dpi=300, format="png")

# %%
# calculate energy balance
energy_total = plate.calc_energy_total()
energy_total_floor = pow2floor(energy_total)
energy_balance_step = np.diff(energy_total)
energy_balance_step_rel = energy_balance_step / energy_total_floor[:-1]

# %%
# plot energy balance
plt.style.use(["science", "ieee", "std-colors"])
fig, axs = plt.subplots(nrows=1, ncols=2, layout="constrained")
fig.get_layout_engine().set(w_pad=0, wspace=0)
plot_with_eps(energy_balance_step_rel, AUDIO_RATE, axs=axs[0], plot_grid=False, markersize=1)
plot_with_eps(energy_balance_step_rel, AUDIO_RATE, axs=axs[1], markersize=1)
axs[0].set_xlim([0, 0.1])
axs[0].set_xticks([0, 0.05, 0.1], ["0", "50", "100"])
if np.max(np.abs(energy_balance_step_rel)) < 1.3e-13:
    axs[0].set_ylim([-1.3e-13, 1.3e-13])
    axs[0].set_yticks([-1e-13, 0, 1e-13])
axs[0].ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
axs[0].set_xlabel("Time [ms]")
axs[0].set_ylabel("Energy error", labelpad=0)
axs[1].set_xlim([0, 0.1])
axs[1].set_xticks([0, 0.05, 0.1], ["0", "50", "100"])
axs[1].set_ylim([-1.3e-15, 1.3e-15])
axs[1].set_yticks([-1e-15, 0, 1e-15])
axs[1].ticklabel_format(style="sci", axis="y", scilimits=(0, 0), useMathText=True)
axs[1].set_xlabel("Time [ms]")
save_fig(out_dir, "energy", fig, width=FA2026_COLUMN_WIDTH, height=(0.55 * FA2026_COLUMN_WIDTH), dpi=300, format="pdf")

# %% [markdown]
# ## Drift Control

# %%
# define function for drift calculation
def calc_drift_rel(plate: VKPlatePSTD):
    psi_V_half = calc_psi(plate.V_half)
    drift = plate.psi_half - psi_V_half
    drift_rel = drift / np.max(psi_V_half)
    return drift_rel

# %%
# simulate plate without and with drift control
exc_kwargs_repeat = copy.deepcopy(exc_kwargs)
exc_kwargs_repeat.update({"num_repeat": 3})
lambda0_list = [0, 1e3]
drift_rel_list = []
for lambda0 in lambda0_list:
    plate_kwargs_lambda0 = copy.deepcopy(plate_kwargs)
    plate_kwargs_lambda0.update({"lambda0": lambda0})
    plate, out, fig = test_plate(plate_kwargs_lambda0, exc_kwargs_repeat)
    drift_rel = calc_drift_rel(plate)
    drift_rel_list.append(drift_rel)
    name = f"drift_lambda0_{lambda0:g}"
    write_wav(out_dir / (name + ".wav"), out, AUDIO_RATE, normalise=True, subtype="PCM_24")
    save_fig(out_dir, name, fig, width=(1920 / 300), height=(1440 / 300), dpi=300, format="png")

# %%
# plot drift
plt.style.use(["science", "ieee", "std-colors"])
fig, axs = plt.subplots(nrows=1, ncols=1, layout="constrained")
fig.get_layout_engine().set(w_pad=0, wspace=0)
for lambda0, drift_rel in zip(lambda0_list, drift_rel_list):
    plot_against_time(drift_rel, AUDIO_RATE, axs=axs, label=num2str(lambda0, var="\\lambda_0"))
exc_times = exc_kwargs_repeat["exc_st"] + plate_kwargs["dur"] * np.arange(start=0, stop=exc_kwargs_repeat["num_repeat"])
xlim = [0, plate_kwargs["dur"] * exc_kwargs_repeat["num_repeat"]]
ylim = [1e-7, 1e1]
axs.vlines(exc_times, ymin=ylim[0], ymax=ylim[1], color="k", linestyle="dashed", linewidth=0.5)
axs.set_xlim(xlim)
axs.set_xticks(np.arange(xlim[0], xlim[1] + 1))
axs.set_ylim(ylim)
axs.set_yscale("log", nonpositive="mask")
axs.set_yticks([1e-7, 1e-5, 1e-3, 1e-1, 1e1])
axs.get_yaxis().set_major_formatter(mpl.ticker.LogFormatterSciNotation())
axs.set_ylabel("Drift", labelpad=0)
axs.grid()
axs.legend(loc="upper center", bbox_to_anchor=(0.5, 1.2), ncols=2, fontsize="small")
save_fig(out_dir, "drift", fig, width=FA2026_COLUMN_WIDTH, height=(0.55 * FA2026_COLUMN_WIDTH), dpi=300, format="pdf")

# %% [markdown]
# ## Spectrograms

# %%
# define HTML table
table = []
headers = [
    "$f_{\mathrm{amp}}$",
    "$\lambda_0$",
    "Output"
]

def create_row(exc_amp, lambda0, name):
    row = [
        num2str(exc_amp),
        num2str(lambda0),
        f"<audio src=\"audio/{name}.wav\" controls></audio>"
    ]
    return row

# %%
# simulate plate under different excitation amplitudes
plate_kwargs_kappa = copy.deepcopy(plate_kwargs)
plate_kwargs_kappa.update({"kappa": 60})
exc_amp_list = [2e5, 2e6, 2e7, 2e7]
lambda0_list = [1e3, 1e3, 1e3, 0]
out_list = []
for exc_amp, lambda0 in zip(exc_amp_list, lambda0_list):
    exc_kwargs_amp = copy.deepcopy(exc_kwargs)
    exc_kwargs_amp.update({"exc_amp": exc_amp})
    plate_kwargs_kappa.update({"lambda0": lambda0})
    plate, out, fig = test_plate(plate_kwargs_kappa, exc_kwargs_amp)
    out_list.append(out)
    name = f"spec_exc_amp_{exc_amp:g}_lambda0_{lambda0:g}"
    write_wav(out_dir / (name + ".wav"), out, AUDIO_RATE, normalise=True, subtype="PCM_24")
    save_fig(out_dir, name, fig, width=(1920 / 300), height=(1440 / 300), dpi=300, format="png")
    table.append(create_row(exc_amp, lambda0, name))

# %%
# print HTML table
table_str = tabulate(
    table,
    headers=headers,
    stralign=None,
    tablefmt="unsafehtml"
)
print(table_str)

# %%
# plot spectrograms
plt.style.use(["science", "ieee", "std-colors"])
fig, axs = plt.subplots(nrows=1, ncols=len(out_list), layout="constrained")
for i, (exc_amp, lambda0, out) in enumerate(zip(exc_amp_list, lambda0_list, out_list)):
    plot_spec(out, AUDIO_RATE, freq_range=(0, 4e3), axs=axs[i])
    exc_amp_str = num2str(exc_amp, var="f_{{\\mathrm{{amp}}}}").replace("$", "")
    lambda0_str = num2str(lambda0, var="\\lambda_0").replace("$", "")
    title = "$" + exc_amp_str + ",\\;" + lambda0_str + "$"
    axs[i].set_title(title, fontsize="small")
    axs[i].set_yticks([0, 1000, 2000, 3000, 4000], ["0", "1", "2", "3", "4"])
    if i == 0:
        axs[0].set_ylabel("Frequency [kHz]")
    else:
        axs[i].set_yticklabels([])
        axs[i].set_ylabel(None)
    axs[i].set_xticks(np.arange(0, plate_kwargs_kappa["dur"]))
save_fig(out_dir, "spec", fig, width=FA2026_TEXT_WIDTH, height=(0.55 * FA2026_COLUMN_WIDTH), dpi=300, format="pdf")

# %% [markdown]
## Additional Examples

# %%
# define HTML table
table = []
headers = [
    "$\kappa$",
    "$\sigma_0$",
    "$\sigma_1$",
    "$\eta$",
    "$x_{\mathrm{e}}$",
    "$y_{\mathrm{e}}$",
    "$x_{\mathrm{o}}$",
    "$y_{\mathrm{o}}$",
    "$f_{\mathrm{amp}}$",
    "$T_{\mathrm{e}}$",
    "Output"
]

def create_row(plate_kwargs, exc_kwargs, name):
    row = [
        num2str(plate_kwargs["kappa"]),
        num2str(plate_kwargs["sigma0"], num_digits=1),
        num2str(plate_kwargs["sigma1"]),
        num2str(plate_kwargs["ratio"], num_digits=1),
        num2str(plate_kwargs["xe"], num_digits=2),
        num2str(plate_kwargs["ye"], num_digits=2),
        num2str(plate_kwargs["xo"], num_digits=2),
        num2str(plate_kwargs["yo"], num_digits=2),
        num2str(exc_kwargs["exc_amp"]),
        num2str(exc_kwargs["exc_dur"]),
        f"<audio src=\"audio/{name}.wav\" controls></audio>"
    ]
    return row

# %%
# simulate plate for different parameters
np.random.seed(0)
plate_kwargs_example = copy.deepcopy(plate_kwargs)
exc_kwargs_example = copy.deepcopy(exc_kwargs)
dur_list     = [6,    5,    4   ]
kappa_list   = [10,   15,   30  ]
sigma0_list  = [1,    1.3,  1.5 ]
sigma1_list  = [1e-4, 2e-4, 4e-4]
ratio_list   = [1.4,  1.2,  1.1 ]
exc_amp_list = [1e6,  2e6,  3e6 ]
exc_dur_list = [3e-3, 2e-3, 3e-3]
xe_list = np.random.uniform(size=len(ratio_list)) * np.sqrt(ratio_list)
ye_list = np.random.uniform(size=len(ratio_list)) / np.sqrt(ratio_list)
xo_list = np.random.uniform(size=len(ratio_list)) * np.sqrt(ratio_list)
yo_list = np.random.uniform(size=len(ratio_list)) / np.sqrt(ratio_list)
for dur, kappa, sigma0, sigma1, ratio, xe, ye, xo, yo, exc_amp, exc_dur \
    in zip(
        dur_list,
        kappa_list,
        sigma0_list,
        sigma1_list,
        ratio_list,
        xe_list,
        ye_list,
        xo_list,
        yo_list,
        exc_amp_list,
        exc_dur_list
    ):
    plate_kwargs_example.update({
        "dur": dur,
        "kappa": kappa,
        "sigma0": sigma0,
        "sigma1": sigma1,
        "ratio": ratio,
        "xe": xe,
        "ye": ye,
        "xo": xo,
        "yo": yo
    })
    exc_kwargs_example.update({"exc_amp": exc_amp, "exc_dur": exc_dur})
    plate, out, fig = test_plate(plate_kwargs_example, exc_kwargs_example)
    name = f"example_kappa_{kappa:g}"
    write_wav(out_dir / (name + ".wav"), out, AUDIO_RATE, normalise=True, subtype="PCM_24")
    save_fig(out_dir, name, fig, width=(1920 / 300), height=(1440 / 300), dpi=300, format="png")
    table.append(create_row(plate_kwargs_example, exc_kwargs_example, name))

# %%
# print HTML table
table_str = tabulate(
    table,
    headers=headers,
    stralign=None,
    tablefmt="unsafehtml"
)
print(table_str)

# %%
