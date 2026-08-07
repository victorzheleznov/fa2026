<h2 style="font-size: 1.5em" align="center">
  Explicit and Stable Pseudospectral Time-Domain Method for the Föppl–von Kármán Equations
</h2>

<p style="font-size: 1.0em" align="center">
  Victor Zheleznov and Stefan Bilbao
</p>

<p style="font-size: 1.0em" align="center">
  Accompanying repository for the FA2026 paper
</p>

<div align="center">

  [![Sound Examples](https://img.shields.io/badge/Sound_Examples-blue)](https://victorzheleznov.github.io/fa2026/)
  [![arXiv](https://img.shields.io/badge/arXiv-2608.06139-b31b1b.svg)](https://arxiv.org/abs/2608.06139)
  
</div>



## Repository Contents

`cfg/plate.yaml` is a configuration file for a plate.

`notebooks/plot_fa2026.py` is a notebook for reproducing results in the paper.

`out/2026_07_16_fa2026/2026_07_16_plot_fa2026.html` is an archived notebook used for results in the paper.

`src/` includes source code for plate models and other utils.

`run_generator.py` is a script used for simulation.



## Instructions

[Python 3.11.9](https://www.python.org/downloads/release/python-3119/) was used for simulations.
The required packages are provided in the `requirements.txt` file. To setup the environment, use:
```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

To reproduce results in the paper, run the `plot_fa2026.py` notebook in a [VSCode's Python Interactive window](https://code.visualstudio.com/docs/python/jupyter-support-py). The result should match the archived `2026_07_16_plot_fa2026.html` notebook which can be opened in any web browser.

To run plate models for other simulation parameters, use the `run_generator.py` script by overriding values in the `plate.yaml` configuration file using the [Hydra interface](https://hydra.cc/docs/tutorials/basic/your_first_app/config_file/). A plate model from the paper with a stiffness parameter `kappa=60` can be simulated as:
```
python -m run_generator generator.kappa=60
```
Output for a linear plate model can be obtained as:
```
python -m run_generator generator.kappa=60 generator._target_=src.generators.LinearPlateModal
```
Simulation results will be saved within the `out/` folder.
