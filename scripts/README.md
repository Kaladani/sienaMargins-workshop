# The meta-analysis scripts

Run everything with the working directory set to the repository root (the folder
holding `index.qmd`), not from inside `scripts/`. Paths are relative to that
root: `data/`, `results/`, `output/`.

```r
setwd("path/to/workshop")     # or open the folder as an RStudio project
source("scripts/descriptives.R")
```

## Start with Chapter 4, not with these scripts

**Chapter 4** (`04-minimal-meta-analysis.qmd`) is the version to copy: three
classes, one specification, one marginal effect, pooled — the whole arc, with three
places marked for your own groups, specification and target. It is a `.qmd`, so you
can open it in RStudio and run the chunks one at a time.

The four scripts here are the same thing at full size: 17 classes, six
specifications, three ways of averaging, plus the descriptive counterparts. They
are worth **reading** once the small version makes sense — they show what the extra
bookkeeping is for — but they are not the easiest place to start editing. Chapters 2
and 3 cover what the quantities mean; this folder is about how they were computed at
scale.

## Start wherever you like

The results of each stage ship with the repository, so you do not have to run the
expensive parts to see the later ones work.

| start at | you need | cost | gets you |
|---|---|---|---|
| `estimation.R` | `data/` | hours | fits for 17 classes × 6 specifications |
| `descriptives.R` | `results/estimation_results.rds` | seconds | the model-free counterparts |
| `Postestimation.R` | `results/estimation_results.rds` | minutes to hours | per-class marginal effects |
| `meta_analysis.R` | the three `.rds` above | seconds | the pooled tables and plots in Chapter 5 |

So to try your own pooling or plotting, start at `meta_analysis.R`. To try your
own marginal-effect targets, start at `Postestimation.R`. Only change the
specifications themselves if you want to re-estimate.

## Run order

1. **`estimation.R`** — fits every class under six specifications, checks
   convergence, re-runs the classes that need it, and writes
   `results/estimation_results.rds`.

   One block per specification, each followed by its convergence check, so you
   can run a single model and stop. The notes about which classes needed a second
   pass are left in place.

2. **`descriptives.R`** — network descriptives per class and wave, the Jaccard
   index between waves, and the *descriptive differences*: the observable,
   uncontrolled counterpart of each marginal effect. No estimate enters these
   quantities; the script reads `estimation_results.rds` only to get the data
   objects and the list of classes back. Writes into `results/` and `output/`.

3. **`Postestimation.R`** — predicted probabilities, first and second differences
   per class, for each model. Writes per-part files into `results/postest_parts/`
   and the combined frame to `results/postest_all_df.rds`, and the same pair for
   the model-implied descriptive counterpart (`postest_descr_parts/`,
   `postest_descr_df.rds`).

   Defaults are set in a block at the top of the file — which models, which
   estimands, which conditionings. Edit those directly. The same settings can
   also be passed as environment variables if you prefer to run it from a shell:

   ```
   MODELS=0,2 RUN_DYNAMIC=FALSE Rscript scripts/Postestimation.R
   ```

   Useful because the estimands differ enormously in cost: the observed-wave
   (`static`) quantities take seconds, while the simulated (`dynamic`) and
   `accumulated` ones take minutes per model.

4. **`meta_analysis.R`** — pools the per-class results with `metafor`, and writes
   the tables and plots that Chapter 5 reads, into `output/`.

## What you need installed

```r
remotes::install_github("stocnet/rsiena@sienaMargins")
install.packages(c("metafor", "ggplot2", "sna"))
```

The `sienaMargins` branch is required: `predict()`, `marginalEffects()` and the
target-construction helpers do not exist on the release version.

## Layout

```
data/      where the friendship networks and covariates go — not included
scripts/   these scripts
results/   .rds produced by estimation, descriptives and postestimation
output/    tables (.csv) and plots (.png) produced by descriptives.R and
           meta_analysis.R
```

`estimation.R` is the only script that reads the raw data, and it fetches them
on first use via `download_data.R`. The other three work from `results/`, which
ships. See `data/README.md` for where the data come from and on what terms.
