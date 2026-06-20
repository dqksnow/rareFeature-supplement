# Supplement Simulation Scripts

This repository contains the simulation scripts prepared for the "Preserving Rare Features in Big Data Regression: Balanced
Subsampling" paper.

- `rareF_main.R`: main implementation wrapper
- `rareF_functions.R`: helper functions used by all designs
- `rareF_simu_HPC_supplement.R`: synthetic simulation driver
- `rareF_supplement.sh`: SLURM launcher for the synthetic simulation
- `rareF_simu_combine_HPC_supplement.R`: combine-comparison simulation driver
- `rareF_combine_compare_supplement.sh`: SLURM launcher for the combine-comparison simulation
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

## Combine-Comparison Simulation

Upload the `scripts/` folder to the server and run:

```bash
cd scripts
sbatch rareF_combine_compare_supplement.sh
```

This submits a 27-task array for the main-text side-by-side aggregation comparison:

- `n_pl in {300, 1000, 10000}`
- `d_u = 5`
- `d_r in {2, 5, 10}`
- `rho in {0.001, 0.005, 0.4}`

The rerun fixes the `BL` pilot and uses the `R-Lopt(BL)` second-step rule, and
saves results for:

- second-step sample
- union sample
- pooled sample
- weighted estimator combination

The script writes one dated `.Rdata` file per scenario to:

```text
raw_results/<mmddyyyy>/combine_compare/
```

Expected output names look like:

```text
combineCompare_plt300_dcont5_drare10_prare0.001.Rdata
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

The combine-comparison rerun produces scenario files for the main-text
aggregation comparison.

The CoverType rerun produces raw files for the first real-data study.

The manuscript figures and result tables were made by running the replications
above and then reading the raw results from the output files. The raw output
files are not included in the repository due to their large file sizes. Here are
the settings corresponding to the main simulation tables and figures in the main
manuscript.

- Table 2: synthetic simulation with `N = 2e5` and
  `d_u = 5`. The displayed settings are:
  - `n_pl = 300`, `d_r = 2`, `rho = 0.01`, from
    `plt300_dcont5_drare2_prare0.01.Rdata`
  - `n_pl = 1000`, `d_r = 2`, `rho = 0.001`, from
    `plt1000_dcont5_drare2_prare0.001.Rdata`
  - `n_pl = 1000`, `d_r = 10`, `rho = 0.001`, from
    `plt1000_dcont5_drare10_prare0.001.Rdata`

- Table 3: synthetic simulation with `N = 2e5` and `d_u = 5`. The displayed
  settings are:
  - `n_pl = 300`, `d_r = 2`, `rho = 0.001`, from
    `plt300_dcont5_drare2_prare0.001.Rdata`
  - `n_pl = 1000`, `d_r = 10`, `rho = 0.01`, from
    `plt1000_dcont5_drare10_prare0.01.Rdata`

- Table 4: CoverType experiment with `n_pl in {300, 1000}` and second-step
  sample sizes `5000`, `10000`, and `50000`, generated from
  `plt300_group1.Rdata` through `plt300_group40.Rdata` and
  `plt1000_group1.Rdata` through `plt1000_group40.Rdata`.

- Table 5: CoverType experiment with the same settings and raw files as
  Table 4.

- Figure 1: synthetic simulation with `N = 2e5`, `n_pl = 300`,
  `d_u = 10`, and `rho = 0.01`.
  - Panel (A): `d_r = 2`, from `plt300_dcont10_drare2_prare0.01.Rdata`
  - Panel (B): `d_r = 10`, from `plt300_dcont10_drare10_prare0.01.Rdata`

- Figure 2: combine-comparison simulation with `N = 2e5`,
  `n_pl = 1000`, `d_u = 5`, `d_r = 2`, and `rho = 0.001`. This uses the
  `R-Lopt` second-step rule with the `BL` pilot, from
  `combineCompare_plt1000_dcont5_drare2_prare0.001.Rdata`.

- Figure 3:
  - Panel (A): CoverType experiment with `n_pl = 1000`, generated from
    `plt1000_group1.Rdata` through `plt1000_group40.Rdata`
  - Panel (B): IRIS Registry experiment. These data are restricted from public
    release, as mentioned in the manuscript. The implementation is very similar
    to the simulation and CoverType implementations in this repository.
