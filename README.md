# Pipelines and scripts used for analyses in the BEAM method paper.

Most inference runs on simulated or real datasets were done using the snakemake pipelines in `pipelines/`. Pipeline dependencies are in `inputs/` and `scripts/`. The `scripts/` directory also contains some standalone scripts.

The snakemake pipelines should be run from within an environment with `graphposterior` installed via pip, since many functions are called in the pipeline from that supplementary package. The pipelines will not run without `graphposterior` installed from [https://github.com/CshlSiepelLab/graphposterior](https://github.com/CshlSiepelLab/graphposterior).
