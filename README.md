# Supplement Simulation Scripts

This repository contains the simulation scripts prepared for the supplementary material of the "Preserving Rare Features in Big Data Regression: Balanced Subsampling" paper.

- `rareF_main.R`: main estimation wrapper
- `rareF_functions.R`: helper functions used by all designs
- `rareF_simu_HPC_supplement.R`: synthetic simulation driver
- `rareF_supplement.sh`: SLURM launcher for the synthetic simulation
- `rareF_Covertype_HPC_supplement.R`: CoverType simulation driver
- `rareF_Covertype_supplement.sh`: SLURM launcher for the CoverType simulation

## Synthetic Simulation

Upload the `scripts/` folder to the server and run:

```bash
cd scripts
sbatch rareF_supplement.sh
```

This submits a 54-task array covering the scenario grid used in the paper:

- `n_pl in {300, 1000}`
- `d_u in {2, 5, 10}`
- `d_r in {2, 5, 10}`
- `rho in {0.001, 0.01, 0.05}`

The script writes one raw `.Rdata` file per scenario to:

```text
raw_results/<mmddyyyy>/logic/
```

Expected output names look like:

```text
plt300_dcont10_drare10_prare0.01.Rdata
```

## CoverType Simulation

Keep `CoverType_gc_full.Rdata` in the same `scripts/` folder, then run:

```bash
cd scripts
sbatch rareF_Covertype_supplement.sh
```

This submits an 80-task array:

- `n_pl in {300, 1000}`
- `40` replicate groups for each pilot size

The script writes grouped raw `.Rdata` files to:

```text
raw_results/<mmddyyyy>/CoverType/
```

Expected output names look like:

```text
plt300_group1.Rdata
plt1000_group40.Rdata
```

## Expected Outputs

The synthetic rerun produces raw scenario files for the simulation study.

The CoverType rerun produces raw grouped files for the empirical study.

These raw `.Rdata` files can then be used locally to generate the supplementary summaries, tables, and figures.
