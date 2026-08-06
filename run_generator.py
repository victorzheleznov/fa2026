from pathlib import Path
from time import time

import hydra
import numpy as np
import resampy
from hydra.utils import instantiate
from omegaconf import OmegaConf

from src.utils.const import AUDIO_RATE
from src.utils.io import write_wav


@hydra.main(version_base=None, config_path="cfg", config_name="plate")
def main(cfg):
    # synthesise
    print("Initialising...")
    generator = instantiate(cfg.generator)
    print("Simulating...")
    start_time = time()
    out = generator(**cfg.excitation)
    elapsed_time = time() - start_time

    # create output directory
    save_dir = Path("out") / generator.__class__.__name__
    save_dir.mkdir(parents=True, exist_ok=True)

    # save
    OmegaConf.save(cfg, save_dir / "cfg.yaml", resolve=True)
    np.save(save_dir / "out.npy", out)
    write_wav(save_dir / "out.wav", out, cfg.generator.fs, normalise=True, subtype="PCM_24")
    if cfg.generator.fs != AUDIO_RATE:
        out_ar = resampy.resample(out, cfg.generator.fs, AUDIO_RATE, filter="kaiser_best")
        write_wav(save_dir / "out_ar.wav", out_ar, AUDIO_RATE, normalise=True, subtype="PCM_24")

    # print
    print(f"Generator: {generator.__class__.__name__}")
    print(f"Elapsed time: {elapsed_time} seconds")
    print(f"Real-time factor: {elapsed_time / cfg.generator.dur * 100:.2f}%")
    if hasattr(generator, "num_nodes"):
        print(f"Number of nodes: {generator.num_nodes}")
    elif hasattr(generator, "num_modes"):
        print(f"Number of modes: {generator.num_modes}")


if __name__ == "__main__":
    main()
