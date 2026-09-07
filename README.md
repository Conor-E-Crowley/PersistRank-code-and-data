# PersistRank: persistence-based spatial conservation prioritization

PersistRank is a reproducible, restartable framework for estimating species
persistence under habitat loss and ranking cells according to their marginal
contribution to persistence. This repository contains the complete workflow,
the Madagascar application, and comparisons with four external Zonation
benchmarks. The application combines terrestrial mammals and birds represented
by binary potential-distribution models produced with point-process modeling
(PPM) or range bagging (RangeBag). The same code supports another named
application when its species, trait, distribution, land-cover, and study-area
inputs satisfy the documented contracts.

The workflow has two deliberately different kinds of products. Expensive
scientific states—posterior draws, persistence points, fitted functions, patch
and population-unit states, Stage 6 removal states, and exact benchmark
lookups—are durable and contract-validated. Persistence trajectories,
statistics, formatted tables, and most report objects are regenerated from
those states so presentation settings can change without repeating spatial
reconstruction.

Numbered R Markdown files are the public interactive entry points. Reusable
implementation lives in `R/`; human-readable application settings live in
`config/`; and the two native kernels live in `src/`. Zonation is not launched
by this repository: Stage 5.2 prepares its inputs, the user runs Zonation
externally, and later stages consume only its `rankmap.tif`.

In code, returned objects and column names sometimes use `pipeline` to identify
the PersistRank sequence. Those internal names are retained as part of the
current programmatic interface; throughout this guide, **PersistRank** refers to
that sequence and **benchmark** refers to a Zonation sequence.

## Quick start

Run R from the repository root so relative paths resolve correctly. The Rmds
are designed for interactive, top-to-bottom chunk execution; knitting is
optional. Start a new R session after changing a configuration, Rmd parameter,
or sourced R module.

The bundled Madagascar analysis already contains:

- a canonical 196-species application table and spatial products for the 188
  species with qualifying population units;
- completed PersistRank removal sequences for all five optimization curves;
- complete ABF, CAZ1, CAZ2, and CAZMAX benchmark libraries, each containing 216
  exact retained-cell states; and
- the scientific states needed to regenerate the report tables and figures
  without rerunning demographic simulation, spatial construction,
  prioritization, or benchmark reconstruction.

Stage 6 recovery checkpoints are not required for a completed run and are not
present in the bundled analysis. The durable removal ledgers, stage lookup
tables, rank surfaces, manifests, and benchmark libraries are the completed
scientific products.

Choose an entry point according to the intended task:

| Task | Entry point and operation | What it does |
|---|---|---|
| Check paths and selected settings | `8_run_pipeline.Rmd`, `mode: inspect` | Reads configuration and small manifests without starting scientific work |
| Inspect one completed stage | The corresponding numbered Rmd in `inspect` mode, where available | Reports readiness and artifact status without modifying the analysis |
| Regenerate the four-benchmark results | `8.1_report_results.Rmd`, then run `setup → results` | Recomputes trajectories, statistics, and tables in memory from completed states |
| Regenerate a completed PersistRank report | `7.1_stage_meta_and_priority_surface.Rmd` | Rebuilds summaries and figures for one selected curve |
| Regenerate the detailed comparison | `7.3_persist_cmp.Rmd`, `mode: report` | Rebuilds the selected PersistRank–Zonation comparison and figure |
| Resume selected scientific stages | `8_run_pipeline.Rmd`, `mode: resume` | Runs only the explicitly selected stages and reuses compatible handoffs |
| Reproduce everything from source inputs | Follow **Genuine from-scratch reproduction** in Section 9 | Fits models, constructs spatial states, runs PersistRank, and rebuilds benchmark libraries |

Do not use **Run All** without first reviewing the Rmd parameters. In
particular, Stage 1 defaults to `fit`, and Stage 8 defaults to `resume` with a
multi-stage selection. Stages 2, 5, 6, and first-time benchmark reconstruction
are computationally expensive. Use `inspect` wherever available before an
active operation.

### Public repository layout

| Path | Role |
|---|---|
| Numbered `*.Rmd` files | Public interactive entry points for the complete workflow |
| `R/` | Reusable workflow, validation, scientific-calculation, recovery, and figure modules |
| `src/` | Native C++ kernels loaded only by active Stage 2, Stage 6, or benchmark reconstruction |
| `config/` | Complete application configurations and the copy-and-edit template |
| `Data/Raw/` | Tabular, raster, workbook, and Zonation-template inputs |
| `mammal_ppm_bin/`, `mammal_rangebag_bin/`, `bird_ppm_bin/`, `bird_rangebag_bin/` | Binary potential-distribution rasters used by the Madagascar application |
| `Data/Models/` | Durable demographic and persistence-model artifacts |
| `Data/Applications/` | Application manifests, spatial states, PersistRank runs, and benchmark libraries |
| `Data/Cache/` | Regenerable query or compiled-code caches; not scientific authority |
| `Figures/` | Regenerable workflow figures |

## 1. Scientific scope and terminology

- A **quasi-extinction threshold** is the configured abundance below which a
  population is treated as not persisting. It is 500 individuals in the
  bundled `qe_500` application.
- A **persistence function** maps abundance above that threshold to the
  probability of persistence over the configured time horizon. The five fitted
  functions are posterior quantile functions `q025`, `q16`, `q50`, `q84`, and
  `q975`.
- **Area of habitat (AOH)** is the intersection of a species' aligned binary SDM
  with land-cover classes mapped from its suitable IUCN habitat labels.
- A **patch** is a rook-contiguous component of AOH that strictly exceeds the
  species-specific minimum patch-area threshold.
- A **population unit (PU)** is one or more retained patches connected within
  the species-specific dispersal distance. A PU is retained only when its total
  area strictly exceeds the area corresponding to the quasi-extinction
  abundance.
- **PU persistence** is evaluated from PU area, density, the population-area
  threshold, and fitted Gompertz coefficients. Under the default independent-PU
  model, **species persistence** is one minus the product of PU failure
  probabilities.
- An **optimization curve** is the fitted function whose coefficients Stage 6
  uses when scoring candidate removals. A completed spatial sequence can later
  be evaluated under all five functions.
- A **pruning iteration** removes one configured batch of frontier cells. An
  **ecological stage** contains up to the configured number of pruning
  iterations followed by fragmentation and distance-connectivity repair.
- A **benchmark method** is one of ABF, CAZ1, CAZ2, or CAZMAX. Each method has a
  separate rank map and exact-state library.
- An **exact target** is an integer retained-cell count. Within one method and
  removal run, that count is the identity of a reconstructed benchmark
  patch/PU state. Equal targets required by different curves share one file.

### Artifact levels

| Level | Identity | Main contents | Principal consumers |
|---|---|---|---|
| Demographic model | calibration data and Stage 1 sampler contract | four posterior tables and mammal/bird trait grids | Stage 2 |
| Persistence model | threshold, horizon, simulation grid, and Stage 3 span | mammal/bird persistence points and five fitted functions | Stage 4 onward |
| Spatial application/scenario | application name, species/SDM contract, threshold, land cover, study area, and patch rules | species table, patch rasters, lookup, connectivity, and one shared Stage 6 initialization | Stages 5.1–8.1 |
| Removal run | application scenario and removal schedule, with an independent result for each optimization curve | curve-specific Stage 6 state sequences, benchmark libraries, and reports | Stages 7–8.1 |

## 2. Chronological data flow

1. **Stage 1** turns raw mammal and bird demographic calibration tables into
   posterior draws for growth and environmental variation.
2. **Stage 2** propagates those draws through stochastic population simulations
   and estimates the abundance needed for configured persistence probabilities.
3. **Stage 3** fits the five reusable trait-dependent Gompertz persistence
   functions to the Stage 2 points.
4. **Stage 4** resolves application species, names, traits, SDMs, habitats, area
   thresholds, and Stage 3 coefficients into one canonical species table.
5. **Stage 5** intersects SDMs with suitable land cover, constructs patches and
   dispersal-connected PUs, and stores the common spatial starting state.
6. **Stages 5.1–5.3** optionally explain that spatial process, prepare external
   Zonation inputs, and report initial persistence.
7. **Stage 6** first converts the Stage 4–5 state into one curve- and
   schedule-neutral initialization, then selects one fitted function at a time
   to generate independent reverse-removal sequences for all requested curves.
8. **Stage 7.1** reports each completed Stage 6 sequence. In parallel, external
   Zonation uses Stage 5.2 features to produce four benchmark rank maps.
9. **Stage 7.2** reconstructs each Zonation sequence only at the exact
   retained-cell counts needed to compare it with completed Stage 6 curves.
10. **Stages 7.3 and 7.4** produce detailed and cross-curve comparisons;
    **Stage 8.1** performs the four-method reporting workflow.

Stage 8 is an orchestrator for a literal selection of Stages 2–7.1; it is not a
new scientific transformation and never runs Stage 1.

## 3. Software and system requirements

### R packages

Production code directly uses `data.table`, `Rcpp`, `dplyr`, `readr`, `readxl`,
`tibble`, `stringr`, `terra`, `sf`, `s2`, `units`, `igraph`, `Rfast`,
`fastmatch`, `ggplot2`, `tidyterra`, `cowplot`, `scales`, `rnaturalearth`,
`rredlist`, `rjags`, `coda`, `rmarkdown`, and `knitr`. `fasterRaster` is needed
only when that Stage 5 clumping backend is selected. The lightweight `yaml`
package reads the shared application configuration. Base packages such as
`stats`, `utils`, `tools`, `grid`, and `grDevices` ship with R.

The workflow never installs packages automatically. After installing the
system software below, missing R packages can be installed from an R session:

```r
required <- c(
  "coda", "cowplot", "data.table", "dplyr", "fastmatch", "ggplot2",
  "igraph", "knitr", "readr", "readxl", "Rcpp", "Rfast", "rjags",
  "rmarkdown", "rnaturalearth", "rredlist", "s2", "scales", "sf",
  "stringr", "terra", "tibble", "tidyterra", "units", "yaml"
)
install.packages(setdiff(required, rownames(installed.packages())))
```

Install `fasterRaster` as well when using the `fasterRaster` clumping backend.
Package installation can require the JAGS, compiler, geospatial, OpenMP, and
GRASS prerequisites described below. A successful package installation alone
does not confirm that these external runtimes are available.

### System software

- R plus a working C/C++ compiler toolchain.
- JAGS, discoverable by `rjags`, when Stage 1 runs in `fit` mode. Reuse mode
  validates saved posterior tables but still checks the configured workflow
  prerequisites.
- OpenMP for the Stage 2 simulator. Linux GCC installations commonly provide
  it; macOS normally needs an OpenMP-capable compiler and runtime. The Stage
  6/7.2 shared native kernel requires C++11 but does not declare OpenMP.
- Geospatial system libraries required by `terra`, `sf`, `s2`, and `units`.
- GRASS GIS plus `fasterRaster` only for `clump_backend: "fasterRaster"`.
  `clump_backend: "terra"` does not require GRASS.
- Zonation for the external ABF, CAZ1, CAZ2, and CAZMAX runs. The R workflow
  prepares inputs and consumes outputs but does not invoke the program.

### IUCN Red List access

Stage 4 can read previously resolved habitat records from the application/cache
or query the IUCN Red List API, depending on `iucn_mode`. A live query requires
an API key in the environment variable `IUCN_REDLIST_KEY`. Set it outside the
configuration and do not commit it to the repository. For example, in a shell
used to launch R:

```sh
export IUCN_REDLIST_KEY="your-key"
```

`cache_only` never requests new records. `cache_or_query` uses available cached
records and queries only unresolved species. `refresh` requests current records
for every selected species. Stage 4 stops before querying if a required key is
absent or invalid.

### Installation check

From the repository root, the following read-only check loads the configuration
infrastructure, validates the complete YAML schema, resolves configured paths,
and counts the four standard distribution-raster folders without opening their
rasters:

```r
source(file.path("R", "project_loader.R"))
project_source(c(
  "R/project_utils.R", "R/project_artifacts.R", "R/project_paths.R",
  "R/project_config.R"
))
cfg <- read_project_config("config/madagascar.yml")
inspect_project_config(cfg)
```

This check does not test JAGS, spatial libraries, GRASS, native compilation, or
Zonation. Their availability is validated when an operation that needs them is
started.

`Data/Cache/Rcpp` is a platform-specific, regenerable compilation cache. It is
not a scientific input, is not referenced by manifests as scientific identity,
and should not be transferred between operating systems. Stage 2 compiles its
simulator only for active simulation; Stage 6 and benchmark reconstruction load
their native kernels only for active spatial work. Inspection and completed
report paths deliberately avoid those costs where the workflow supports it.
`Data/Cache/iucn_habitat_cache.rds` is likewise a regenerable query accelerator;
the resolved rows and fingerprints installed by Stage 4, rather than the cache
file itself, become the application authority.

## 4. Execution model

Each numbered Rmd begins by sourcing one workflow module, reading its small set
of operational options through `project_rmd_params()`, then merging them with
the configuration selected by `config_file`. `config/madagascar.yml` is the
complete bundled scientific/input configuration; `config/example.yml` is a
fully enumerated template for another application. Copy the example, change
its values, and point each Rmd used for that application at the copy. The
configuration supplies scientific identity, input paths, spatial settings,
removal schedule, and report defaults. The Rmd supplies only the operation:
mode, selected stages/curves/method, machine-specific threads/backend, or which
report to generate. An Rmd cannot silently override a scientific configuration
key.

`R/project_config.R` rejects missing and unknown configuration
keys. `inspect_project_config()` resolves every configured input path and
counts the four standard SDM folders without opening tables or rasters,
computing checksums, loading spatial/native packages, or writing files. Stage
8 displays this inspection alongside its execution plan. Once durable work
exists, its manifests and checkpoints remain authoritative at resume
boundaries; editing a configuration does not make incompatible saved state
silently reusable.

Run an Rmd's setup or first execution chunk after starting a new R session or
changing the configuration or its Rmd options. The documents are designed for
interactive chronological chunk execution; knitting is not required. Their
ordinary `html_document` declarations provide a common R Markdown format, not a
separate batch execution contract.

Public relative paths are resolved from the repository root. Supported inputs
also accept absolute paths and `~`. `R/project_loader.R` sources each normalized
module at most once per session into the global environment required by the
Rmds. After changing an R module interactively, start a fresh R session or call
`project_source("R/file.R", reload = TRUE)` before rerunning dependent code.

Configuration, inspection, and reporting are intentionally separated from
expensive execution. `inspect` modes read only the small metadata needed for
status. `resume` modes validate native recovery state before continuing.
Report-only stages read authoritative Stage 6 or exact benchmark states and do
not create a second scientific cache.

## 5. Authoritative inputs

### Raw files

| Path | Format and required content | Used by | Role |
|---|---|---|---|
| `Data/Raw/mammal_rmax.txt` | tab-delimited mammal table containing `log10(M)`, `log10(rm)`, `Order`, and `Species` | Stage 1 | Mammal maximum-growth calibration; bats, enumerated marine/aquatic species, and invalid rows are excluded before fitting. |
| `Data/Raw/sigma.csv` | CSV containing at least `Class`, `Genus`, `Species`, `Vr`, and `Mass` | Stage 1 | Mammal and bird environmental-variation calibration; standard deviation is `sqrt(Vr)`. |
| `Data/Raw/bird_growth_niel_lebreton.csv` | CSV containing `Species` and `lambda` | Stage 1 | Bird growth calibration; growth is converted to `log(lambda)` and matched to generation length. |
| `Data/Raw/bird_data.txt` | tab-delimited EltonTraits bird records containing `Scientific`, `Diet-5Cat`, body mass, and diet fields | Stages 1 and 4 | Supplies bird diet branches for Stage 1 and bird mass/diet traits for the application species table. |
| `Data/Raw/mammal_data.txt` | tab-delimited EltonTraits mammal records containing `Scientific`, body mass, and diet fractions | Stage 4 | Supplies mammal mass and diet traits. |
| `Data/Raw/cobi13486-sup-0004-tables4.xlsx` | Bird et al. generation-length workbook | Stages 1 and 4 | Supplies generation lengths used in bird demographic fitting and species-specific Stage 3 coefficient prediction. |
| `Data/Raw/iucn_bird_synonyms.csv` | IUCN bird synonym export containing `scientificName`, `genusName`, and `speciesName` | Stage 1 | Resolves historical bird growth/variation names to generation-length and diet records. |
| `Data/Raw/simple_summary.csv` | application species inventory with scientific name, taxonomic hierarchy, and Red List category | Stage 4 | Defines the candidate Madagascar species set and its canonical order. |
| `Data/Raw/synonyms.csv` | general synonym table with accepted scientific, genus, and species names | Stage 4 | Provides ordinary synonym candidates for SDM, trait, and generation-length matching. |
| `Data/Raw/input_synonyms.csv` | curated pairs `scientificName` and `synonym_scientificName` | Stages 1 and 4 | Supplies explicit fallback/override mappings when general sources do not resolve an input name uniquely. |
| `Data/Raw/synonyms_reasoning` | UTF-8 Markdown-style notes and supporting links | Human provenance | Records the reasoning behind selected curated historical-name mappings; production code does not read it. |
| `Data/Raw/random_effects.csv` | class/order/family/species random-effect table | Stage 4 | Adds taxonomic effects to density calculations; missing matches receive zero. |
| `Data/Raw/esacci_2022_pfts.tif` | categorical ESA CCI land-cover raster with a valid CRS and resolution | Stages 4, 5, and 5.1 | Recorded in the application contract and converted to habitat masks inside the selected study area. |
| `Data/Raw/settings.z5.txt` | Zonation settings template containing exactly one feature-list pointer | Stage 5.2 | Copied with only that pointer replaced to reference the generated feature list. |

### Binary species-distribution rasters

The default folder-based SDM inventory contains 210 single-layer GeoTIFFs in
four directories:

| Directory | Taxon | Method | Included rasters |
|---|---|---|---:|
| `mammal_ppm_bin/` | mammals | PPM | 43 |
| `mammal_rangebag_bin/` | mammals | RangeBag | 64 |
| `bird_ppm_bin/` | birds | PPM | 61 |
| `bird_rangebag_bin/` | birds | RangeBag | 42 |

Filenames use `<Genus>_<species>_bin.tif`. Presence cells must equal `1`; other
values are treated as absence/no data by Stage 5. Stage 4 inventories the
selected taxon/method folders, resolves each retained species to one raster,
and records the absolute source path and match provenance. Stage 5 then projects
or nearest-neighbour resamples each selected raster to the common land-cover
template if required.

For a portable or nonstandard collection, set `sdm_index_file` to a CSV that
identifies scientific name, taxon, method, and raster path. Relative raster
paths inside that CSV resolve beside the index file. A vector study area may be
supplied with `study_area_mode: "vector"`, `study_area_file`, and optionally
`study_area_layer`; bounds mode uses the four ROI coordinates, while
`full_raster` uses the land-cover extent.

## 6. Canonical hierarchy

```text
Data/
├── Raw/
├── Cache/
│   ├── iucn_habitat_cache.rds
│   └── Rcpp/
├── Models/
│   ├── Demography/{Outputs,Checkpoints}/
│   └── Persistence/qe_<threshold>[_horizon_<years>yr]/{Stage2,Stage3,Logs}/
└── Applications/<application>/
    ├── application_manifest.rds
    ├── application_inputs.rds
    └── qe_<threshold>[_horizon_<years>yr]/
        ├── scenario_manifest.rds
        ├── Species/
        ├── Spatial/
        │   ├── Patches/
        │   ├── Checkpoints/
        │   └── priority_initialization.rds
        ├── Logs/stage6_initialization.log
        ├── Zonation/
        │   ├── Patches_binary/
        │   ├── feature_list.txt
        │   ├── settings.z5.txt
        │   ├── ABF/rankmap.tif
        │   ├── CAZ1/rankmap.tif
        │   ├── CAZ2/rankmap.tif
        │   └── CAZMAX/rankmap.tif
        └── Runs/remove_<cells>_stage_<iterations>/
            ├── run_manifest.rds
            ├── run_state.rds
            ├── Curves/<curve>/Run/
            │   ├── checkpoints/
            │   ├── removal_events.csv
            │   ├── removal_order.tif
            │   ├── rankmap.tif
            │   ├── patch_lookup_tables/
            │   └── Analysis/stage_meta.csv
            └── BenchmarkLookups/<METHOD>/
                ├── lookup_manifest.rds
                ├── checkpoint.rds
                └── Lookups/retained_cells_<N>.rds
```

Model figures are written under `Figures/Models/`. Application figures mirror
application, threshold, and run under
`Figures/Applications/<application>/<threshold-tag>/<run-tag>/`, with
`Spatial`, `Priority`, and `Comparisons` children.

## 7. Artifact ownership and invalidation

| Owner | Durable identity | Invalidated by |
|---|---|---|
| Stage 1 | demographic calibration | demographic data, priors, chain settings |
| Stages 2–3 | threshold and horizon model | threshold, horizon, simulation grid, fitted points |
| Stage 4 | named application species contract | species, names, SDMs, study area, IUCN or trait inputs |
| Stage 5 | threshold-specific spatial source | application contract, land cover, habitat, patch/PU rules |
| Stage 6 initialization | application + threshold-specific Stage 4–5 spatial source | species/spatial source, Stage 4 coefficients, patch/PU rules; not curve or schedule |
| Stage 6 curve execution | shared initialization + run schedule + optimization curve | shared-initialization fingerprint, schedule, selected coefficient pair |
| Stage 7.2 | run + benchmark method + rank-map/spatial contract + exact count | rank map, Stage 5 source, threshold contract, schedule |
| Stages 7.3/7.4/8.1 | no derived scientific cache | regenerated when report settings change |

Changing a focal species, point style, central statistic, late-stage target, or
persistence reporting threshold never invalidates benchmark lookups.

The shared Stage 6 initialization identity is exactly the threshold scenario's
selected taxa/SDM source and analysis contract together with the contents of
`Species/species_table.csv`, `Spatial/all_patch_lookup.rds`,
`Spatial/all_connectivity.rds`, and `Spatial/Patches/`. Those sources determine
the retained species, thresholds, all five coefficient pairs, grid cells,
patches, PUs and graphs. The scenario recovery record fingerprints those inputs
and `Spatial/priority_initialization.rds`; run manifests reference the latter
fingerprint. Optimization curve, cells removed per iteration, iterations per
stage, maximum stages, checkpoint/logging controls, report settings, and the
removal-run directory are deliberately absent from initialization identity.
They affect only curve execution or operations. Thus different curves and
compatible `remove_<cells>_stage_<iterations>` runs resolve to the same path,
while any Stage 4/5 content change fails the shared recovery check.

## 8. Chronological stage manual

Every entry uses the same fields. **Returned** means an object in the current R
session; only **Durable outputs** survive that session. Paths are canonical
application-storage paths unless a standalone low-level context is injected.

### Stage 1 — demographic models

- **Entry point and purpose:** `1_demographic_models.Rmd` calls
  `run_stage1()`. It constructs the mammal and bird demographic allometries
  from which Stage 2 draws growth-rate and environmental-variation values.
- **Inputs:** `mammal_rmax.txt`, `sigma.csv`,
  `bird_growth_niel_lebreton.csv`, `bird_data.txt`, the CoBi
  generation-length workbook, `iucn_bird_synonyms.csv`,
  `input_synonyms.csv`, and the configured bird diet branches. Each file path
  comes from the `inputs` section of the selected configuration. These are
  calibration data, not the application species inventory.
- **Configuration:** `fit` or `reuse`, verbosity, random seed, JAGS chains,
  adaptation/iteration/thinning counts, coefficient and residual-SD priors,
  R-hat and effective-sample-size limits, trait-grid resolution, mammal mass
  bounds, and the bird diet branches with separate intercepts.
- **Ordered calculation:**

  1. Read and type-check every calibration table. Mammal growth rows with
     invalid measurements, Chiroptera, and the enumerated marine/aquatic taxa
     are excluded; mammal variation values are converted from variance to
     standard deviation and invalid or excluded rows are removed.
  2. Normalize bird names, apply the curated synonym mapping, join generation
     lengths, and assign diet branches by deterministic precedence. Ambiguous
     or incomplete matches fail before fitting.
  3. Construct logarithmically spaced mammal-mass and bird-generation-length
     prediction grids. These grids fix the traits at which Stage 2 will later
     simulate persistence.
  4. Fit four log-scale JAGS regressions: mammal maximum growth versus body
     mass; mammal environmental variation versus body mass; bird maximum growth
     versus generation length; and bird environmental variation versus
     generation length with the configured diet intercept structure.
  5. Validate every monitored coefficient against the configured R-hat ceiling
     and effective-sample-size floor. A fit that does not meet the convergence
     contract is not installed.
  6. Summarize posterior coefficients and generate posterior prediction draws
     over the two trait grids. Write the complete output set transactionally so
     a failed fit cannot partially replace a valid model.
- **Returned:** calibration row/diet summaries, posterior draw count and branch
  names, durable output/figure paths, installed manifest, and elapsed time. The
  large calibration, posterior, fit, and plot objects are intentionally not
  retained in the Rmd session.
- **Durable outputs:** four posterior-draw CSVs, mammal and bird trait
  grids, and the demographic manifest under
  `Data/Models/Demography/`, plus allometry and posterior figures under
  `Figures/Models/Demography/`.
- **Downstream handoff:** Stage 2 validates the demographic manifest and reads
  the grids and posterior tables. No spatial stage reads the raw calibration
  files directly.
- **Recovery, invalidation, and cost:** `reuse` never invokes JAGS and succeeds
  only for complete compatible posterior tables. It still reconstructs and
  transactionally writes the configured trait grids, figures, and demographic
  manifest so the installed Stage 1 handoff is complete. A change to calibration rows, diet grouping,
  priors, sampler identity, convergence rules, or prediction-grid identity
  requires Stage 1 and all dependent stages to be rebuilt. Fitting is
  compute-intensive and requires JAGS; reuse does not refit.
- **Rmd order:** `setup` calls the reusable workflow once; `summary` displays
  only the compact returned summaries and installed figures.

### Stage 2 — persistence simulation

- **Entry point and purpose:** `2_persistence_simulation.Rmd` calls
  `run_stage2_model()`, which runs `run_stage2()` and publishes its recovery
  boundary after active work. For each mammal-mass or bird-generation-length grid row, it
  estimates the starting abundance required to reach a requested probability
  of remaining above the quasi-extinction threshold for the full horizon.
- **Inputs:** the installed Stage 1 manifest, posterior tables and trait grids;
  the threshold/horizon contract; and the Stage 2 numerical settings. Stage 2
  never substitutes unvalidated raw demographic rows for this handoff.
- **Configuration:** `inspect`, `resume`, or `restart`; demographic uncertainty
  contract; the exact five curve labels; posterior-draw and stochastic-
  replicate counts; chunk size; base seed; anchor and reporting probabilities;
  population cap and growth buffer; K-search start/tolerance/maximum/rounding;
  and optional OpenMP thread count.
- **Ordered calculation:**

  1. Validate Stage 1 provenance and deterministically select posterior rows.
     Construct all residual-uncertainty streams, when requested, from the base
     seed so a resumed block receives the same draws as an uninterrupted run.
  2. On the first active simulation call, compile or load
     `src/simulate_persist_probs_cpp.cpp`, verify the native interface, and create a
     common-random-number context. Inspection does not compile the kernel.
  3. For every trait row, evaluate abundance at the configured anchor
     probabilities. Bracket each target and search K with the configured
     rounding and relative-tolerance rules, reusing identical stochastic
     innovations across competing K values.
  4. Fit a temporary shifted-Gompertz inversion to the anchor results and use
     it to propose K for the complete regular probability grid and any added
     probabilities. The simulator, rather than the temporary fit, supplies the
     reported persistence values.
  5. After each previously unseen K evaluation, atomically update the compact
     checkpoint. After a complete trait/curve block, promote its rows to the
     partial table so resumption can skip completed work.
  6. Validate uniqueness, monotonic ordering, coverage, provenance and schema
     across all blocks. Only then atomically promote complete mammal and bird
     point tables and remove superseded recovery files.
- **Returned:** inspection or completion status, block progress, validated
  table summaries, native/runtime diagnostics, and the active log path.
- **Durable outputs:** final mammal and bird persistence-point CSVs under
  `Data/Models/Persistence/<persistence-model>/Stage2/`; while incomplete,
  partial tables and per-block checkpoint/state files in the same model tree;
  and append-only `Logs/stage2.log` operational timings.
- **Downstream handoff:** Stage 3 consumes only the two complete point tables
  and their model/provenance contract. Partial files are recovery artifacts,
  never scientific inputs to Stage 3.
- **Recovery, invalidation, and cost:** `resume` reuses every compatible K
  evaluation and completed block; `restart` replaces only Stage 2-derived model
  products. Threshold, horizon, uncertainty identity, probability grid,
  posterior/replicate counts, seed or K-search semantics invalidate Stage 2.
  Thread count and checkpoint granularity affect execution but not an otherwise
  numerically identical scientific contract. This is normally the most
  simulation-intensive model stage.
- **Rmd order:** `run` performs setup/inspection or simulation and publishes the
  model boundary; `branches` and `artifacts` display its compact result.

### Stage 3 — persistence curves

- **Entry point and purpose:** `3_persistence_curves.Rmd` calls
  `run_stage3_model()`, which runs `run_stage3()`, records its model boundary,
  and installs the persistence-model manifest. It turns the discrete Stage 2 simulations into five compact,
  reusable trait-to-persistence functions.
- **Inputs:** the final mammal and bird Stage 2 point tables and their
  threshold, horizon, curve, trait-grid and simulation provenance.
- **Configuration:** `fit` or `reuse`, LOESS span, threshold, horizon, and the
  independent `write_figures` presentation switch.
- **Ordered calculation:**

  1. Validate both point tables as one complete Stage 2 handoff; reject mixed
     model identities, missing trait/curve/probability combinations, duplicate
     rows, or nonpositive abundance values.
  2. For each trait-grid row and each of the five curves, fit the exact
     two-coefficient shifted-Gompertz relationship between starting abundance
     and persistence probability.
  3. Transform the positive coefficients to logarithmic scale and fit LOESS
     models over log10 body mass for mammals and log10 generation length for
     birds, retaining bird diet branches where required.
  4. Generate fitted coefficient predictions and diagnostic figures, validate
     coefficient positivity, curve coverage, schema and threshold/horizon
     identity, then atomically install the reusable model and manifest.
- **Returned:** the input-table summary, fitted-model path, optional figure
  paths, threshold/horizon, bird branches and elapsed time. Prediction grids
  are constructed only when figures request them and are not retained.
- **Durable outputs:** `persistence_curve_models.rds`, its persistence-model
  manifest, and optional model/diagnostic figures under the canonical
  persistence model and figure directories.
- **Downstream handoff:** Stage 4 predicts and stores all five coefficient pairs
  for every application species. Later stages use those frozen species-table
  coefficients, not a refit in the reporting session.
- **Recovery, invalidation, and cost:** compatible `reuse` validates rather than
  refits. Changed Stage 2 points, LOESS span, threshold/horizon or model schema
  require Stage 3 and dependent application artifacts to be rebuilt; changing
  only `write_figures` regenerates presentation files. Fitting is moderate and
  does not invoke the stochastic simulator.
- **Rmd order:** `run`, `inputs`, and `figures`.

### Stage 4 — species table

- **Entry point and purpose:** `4_build_species_table.Rmd` calls
  `run_stage4()`. It resolves the named taxa, SDMs, traits, habitat classes and
  persistence coefficients into the canonical row-per-species application
  contract.
- **Inputs:** `simple_summary.csv`, `synonyms.csv`, `input_synonyms.csv`, mammal
  and bird trait sources, random effects,
  generation lengths, the four SDM collections or an explicit SDM index, the
  Stage 3 persistence model, land-cover and study-area definitions, and frozen,
  cached or newly queried IUCN habitat rows according to `iucn_mode`.
- **Configuration:** `figure` or `build`; application, threshold and horizon;
  minimum patch abundance; taxa and SDM-method selectors; all input paths;
  IUCN cache/query policy and pause; and bounds, vector or full-raster study
  area.
- **Ordered calculation:**

  1. Read the selected species inventory and construct ordered scientific-name
     candidates from the original name, general synonyms and curated overrides.
     Record the reasoning used for accepted mappings.
  2. Discover the standard SDM directories or validate the explicit index.
     Select exactly one raster/method record for each retained taxon and reject
     ambiguous or missing method identity.
  3. Select one trait record per taxon; attach bird generation length; derive
     the diet branch; and calculate density, home range and dispersal from the
     retained allometries and random effects.
  4. Convert abundance requirements to area requirements: minimum patch area is
     `minimum_patch_abundance / density`, and minimum PU area is
     `quasi_extinction_abundance / density`, in the canonical area units.
  5. Predict and attach both shifted-Gompertz coefficients for every one of the
     five persistence functions from the Stage 3 model.
  6. Resolve suitable IUCN habitat classes from frozen rows or the application
     cache; query only when the selected policy permits it. Validate scientific
     names, SDM files, traits, areas, coefficient coverage and habitat rows as a
     single application table.
  7. Atomically install the resolved habitat rows and application/scenario
     manifests, then publish the single canonical `species_table.csv`. The
     application manifest owns input, SDM and study-area fingerprints; there is
     no duplicate species-table metadata artifact.
- **Returned:** build action, row and class/SDM counts, number of queried IUCN
  species, output path and elapsed time. Newly resolved habitat rows are passed
  directly to the lifecycle commit and are not duplicated in the public result.
- **Durable outputs:** the application input inventory and manifests,
  IUCN habitat cache when queries update it, and
  `Species/species_table.csv`. `figure` mode reads only that table, skips
  build-input checksums and writes only the replaceable `area_curve.png`.
- **Downstream handoff:** Stage 5 consumes the species table and frozen input
  contract; Stages 5.3–8.1 use its five coefficient pairs and demographic/spatial
  thresholds. The SDM path recorded here is the path Stage 5 opens.
- **Recovery, invalidation, and cost:** resolved cache rows support offline
  reuse, while `refresh` deliberately requests fresh IUCN data. Any change to
  selected species, accepted names, traits, SDM identity, habitat rows,
  threshold/horizon, minimum patch abundance or study-area/application contract
  invalidates Stage 4 onward. Computation is moderate; permitted web queries can
  make the stage network-sensitive.
- **Rmd order:** `setup`, `run`, and `summary`.

### Stage 5 — patches and connectivity

- **Entry point and purpose:** `5_build_patches_and_connectivity.Rmd` calls
  `run_stage5()`. It converts every selected SDM into threshold-filtered habitat
  patches, population units (PUs), and compact connectivity data shared by all
  optimization curves.
- **Inputs:** the completed Stage 4 species table and manifests, every recorded
  SDM, the PFT land-cover raster, IUCN habitat classes, and the configured study
  area. The stage validates all of these before loading the spatial runtime.
- **Configuration:** `inspect`, `resume`, or `restart`; clumping backend
  (`auto`, `fasterRaster`, or `terra`); optional GRASS installation; verbosity;
  land-cover path; and study-area mode/file/layer/bounds.
- **Ordered calculation for each species:**

  1. Crop or mask the land-cover raster to construct the canonical study-area
     context, map IUCN habitat codes to suitable PFT classes, and create a binary
     habitat mask.
  2. Open and align the chosen SDM to the study-area grid. Intersect positive
     occurrence with suitable land cover to obtain area of habitat (AOH).
  3. Label rook-connected AOH cells, calculate geodesic cell area by patch, and
     retain only patches whose area is **strictly greater than** the species'
     minimum patch area.
  4. Polygonize the surviving patches and connect any two patches within the
     species' dispersal distance. Connected components form provisional PUs.
  5. Sum patch area within each provisional PU and retain only PUs whose total
     area is **strictly greater than** the minimum PU area. Remove all cells in
     rejected PUs.
  6. Renumber surviving patches and PUs into canonical deterministic IDs,
     rebuild the final patch raster, and encode within-PU distance edges as
     compressed-sparse-row connectivity.
  7. Write a species checkpoint containing its raster, lookup rows,
     connectivity and status. Species with no surviving spatial unit are
     recorded explicitly and do not receive a misleading nonempty raster.
  8. After every species is complete, validate cross-species IDs and schemas,
     then atomically promote the per-species rasters, aggregate lookup,
     connectivity object and metadata. Successful promotion supersedes the
     temporary checkpoints.
- **Returned:** a compact action/status summary, selected and checkpoint-reuse
  counts, retained/filtered/error rows, published patch/PU counts, backend,
  elapsed time and active runtime-log path. Completed per-species objects,
  aggregate lookup and connectivity are released after publication.
- **Durable outputs:** one canonical patch-ID raster per spatially retained
  species under `Spatial/Patches/`, aggregate all-patch lookup and connectivity
  artifacts, spatial metadata and finalization records; incomplete work retains
  per-species checkpoints; active work appends operational timings to
  `Logs/stage5.log`.
- **Downstream handoff:** Stage 5.2 converts the patch rasters to Zonation
  features; Stage 5.3 computes baseline persistence; Stage 6 initializes the rasters,
  lookup and connectivity; Stage 7.2 reconstructs benchmark states from the
  same immutable spatial source.
- **Recovery, invalidation, and cost:** `inspect` is read-only. `resume` skips
  compatible species checkpoints and final products; `restart` replaces the
  Stage 5 scenario output. Changes to the Stage 4 contract, SDM pixels,
  land-cover pixels, study-area geometry, habitat mapping, thresholds or any
  spatial rule invalidate Stage 5 onward. This stage is I/O-, geometry- and
  memory-intensive; the optional fasterRaster backend requires GRASS.
- **Rmd order:** `setup`, `run_pipeline`, and `summary`.

### Stage 5.1 — one-species spatial explanation

- **Entry point and purpose:**
  `5.1_single_species_aoh_patches_pu_process.Rmd` calls
  `run_stage51()` to explain, rather than supply, Stage 5 processing for one
  selected species.
- **Inputs and configuration:** the focal species' Stage 4 row, SDM, land cover,
  habitat classes and study area, plus the same clumping backend/GRASS controls
  used by Stage 5.
- **Ordered calculation:** repeat the focal species' aligned SDM, suitable-
  habitat intersection, rook clumping, strict patch-area filtering,
  dispersal-based PU formation, strict PU-area filtering and final surviving
  unit construction; retain the intermediate layers needed for the panel
  explanation.
- **Returned:** focal species, mapped area, intermediate/final unit counts,
  source raster, figure path and elapsed time; intermediate spatial layers are
  released after the figure is written.
- **Durable output:** one atomically replaceable process figure. It does not
  replace or amend Stage 5 rasters, lookups or connectivity.
- **Downstream, recovery, invalidation and cost:** presentation only. Regenerate
  after changing the focal species, spatial source or figure implementation;
  no scientific stage consumes it. A null focal species records this optional
  stage as not applicable; a supplied name must match exactly one validated
  species-table row before spatial initialization begins. Cost is one-species
  spatial processing.
- **Rmd order:** `run` and `summary`.

### Stage 5.2 — Zonation inputs

- **Entry point and purpose:**
  `5.2_build_binary_patch_rasters_zonation_feature_list.Rmd` calls
  `run_stage52()` to prepare the only inputs sent to the external Zonation
  benchmark program.
- **Inputs:** the complete Stage 5 patch-raster set, aggregate lookup and
  spatial/application manifests, plus
  `Data/Raw/settings.z5.txt`.
- **Configuration and modes:** `inspect` validates readiness without writing;
  `inputs` rebuilds the binary feature directory and text inputs as one atomic
  set; `feature_list` opens no source patch rasters and rewrites only
  `feature_list.txt` and the settings file from already existing binary
  features.
- **Ordered calculation:** validate the Stage 5 output set and its species
  ordering; map every positive patch ID to feature value 1 while preserving the
  canonical geometry and NA/background contract; write one binary feature per
  spatial species; create the ordered feature list; derive the settings file
  from the supplied template; and ensure the four method directories exist.
- **Returned:** mode-specific readiness, feature counts, ordered paths and
  validation summaries.
- **Durable outputs:** binary feature rasters, `feature_list.txt`, the generated
  settings file and `ABF`, `CAZ1`, `CAZ2`, `CAZMAX` directories beneath the
  scenario's `Zonation/` directory.
- **Downstream, recovery, invalidation and cost:** Zonation reads the binary
  features/list/settings; R later consumes only each method's `rankmap.tif`.
  Rebuild `inputs` when Stage 5 geometry/species change; use `feature_list` when
  only the text handoff needs refreshing. Raster conversion is moderate;
  inspection and feature-list refresh are light and avoid the spatial runtime
  where implemented.
- **Rmd order:** `setup`, `build_zonation_inputs`, and `summary`; the selected
  mode determines what the build chunk actually opens or writes.

### Stage 5.3 — initial persistence

- **Entry point and purpose:** `5.3_initial_species_persistence.Rmd` calls
  `run_stage53()` to report the unpruned baseline under one selected persistence
  curve.
- **Inputs and configuration:** the Stage 4 species table, Stage 5 aggregate
  patch lookup and a selected curve/application/scenario.
- **Ordered calculation:** aggregate initial patch area to each PU; evaluate PU
  persistence with the species' selected shifted-Gompertz coefficients; combine
  independent PU failure risks into species persistence; retain species without
  spatial units in the fixed denominator; then summarize coverage, taxon
  groups, extremes and associations among traits, area and persistence.
- **Returned:** PU-, species- and assemblage-level tables, summary statistics
  and figure data, all in the current session.
- **Durable output:** only `initial_persistence.png`; the detailed tables are
  deliberately not cached.
- **Downstream, recovery, invalidation and cost:** report-only and inexpensive.
  No later stage reads the figure or in-memory tables. Regenerate after changing
  the selected curve or its Stage 4/5 inputs.
- **Rmd order:** `setup`, `coverage`, `pu_summary`, `species_summary`,
  `correlations`, `species_extremes`, `initial_persistence_figure`, and
  `summary`.

### Stage 6 — PersistRank prioritization

- **Entry point and purpose:** `6_spatial_prioritization_pipeline.Rmd` calls
  the Stage 6 workflow. The same Rmd deliberately exposes two artifact
  boundaries without adding another entry point: shared initialization
  (conceptual Stage 6.1) and curve-specific optimization (conceptual Stage
  6.2).
- **Inputs:** initialization reads the Stage 4 species table and the completed
  Stage 5 patch rasters, aggregate patch lookup and connectivity. Curve
  execution reads only the installed shared initialization,
  its selected coefficient pair, and the run/checkpoint state for that curve.
- **Configuration:** `inspect`, `initialize`, `run`, or `resume`; optimization
  curve for all modes except `initialize`; cells removed per iteration;
  iterations per ecological stage; optional maximum stage; ecological
  diagnostic cadence; checkpoint cadence, retention and explicit resume stage.
- **Shared initialization (Stage 6.1):** `initialize` selects the spatially
  retained species and stores their curve-neutral biological parameters plus
  all five named coefficient pairs. It converts every patch raster once into
  canonical cell/patch arrays, records cell areas and initial alive-species
  counts, builds rook adjacency, per-patch cell indexes and compact CSR PU
  connectivity graphs. The lifecycle transaction fingerprints the Stage 4–5
  sources, validates the current schema and installs exactly one immutable
  `Spatial/priority_initialization.rds` for the scenario. It does not select an
  optimization curve, inspect a removal schedule, create a run, or remove a
  cell.
- **Curve execution (Stage 6.2):** `run` or `resume` validates the shared
  initialization at the scenario recovery boundary, maps the requested curve
  deterministically to its named alpha/beta columns in memory, and records a
  compact coefficient identity. The optimized Stage 6 calculation receives the
  selected coefficients through its `a_pred`/`b_pred` interface; the large
  spatial state remains in the shared scenario-level initialization.
- **Ordered curve calculation:**

  1. Load and validate the shared initialization, select the coefficient pair,
     defer-load the native prioritization kernel,
     and initialize from Stage 0 or a compatible stage checkpoint. Construct the
     frontier of eligible habitat cells and the per-species patch/PU state.
  2. Compute each frontier cell's marginal loss on a log-persistence scale,
     including the configured redundancy aggregation. Select a deterministic
     batch by score, using canonical cell ID to break ties.
  3. Mark selected cells removed, decrement patch areas by their exact cell
     areas, and cascade the strict patch and PU area thresholds. When a patch or
     PU becomes nonviable, remove all habitat it can no longer contribute.
  4. Repair each affected PU graph, then apply rook-fragmentation repair within
     patches and dispersal-distance connectivity repair among patches. These
     repairs can trigger further patch/PU losses, so the cascade continues until
     the ecological state is stable.
  5. Record every direct and cascading removal event with iteration/stage and
     reason. After the configured iterations, commit a positive ecological
     stage and atomically write its canonical patch lookup.
  6. At checkpoint cadence, store the complete resumable state and retain only
     the configured recent checkpoints. Continue until no removable habitat or
     `max_stages` is reached.
  7. Finalize the terminal retained layer, removal-event ledger, removal-order
     raster and normalized rank map, validate their mutual ordering/counts and
     install the run/curve manifest.
- **Returned:** inspection/initialization/run status, initialized species,
  patch and PU counts, completed stages, retained and removed counts,
  workload/hotspot/timing profiles, checkpoint information and runtime-log
  path.
- **Durable outputs:** the scenario owns one immutable
  `Spatial/priority_initialization.rds` and an operational
  `Logs/stage6_initialization.log`. Each curve/run independently owns resumable
  stage checkpoints, a positive-stage lookup CSV for each committed stage, the
  removal-event ledger, terminal retained layer, removal-order and rank
  rasters, curve/run manifest entries, and
  `Runs/<run>/Logs/stage6_<curve>.log`.
- **Downstream handoff:** Stage 7.1 reads the Stage 6 state index for PersistRank
  reports; Stage 7.2 derives exact benchmark targets from all selected curve
  trajectories; Stages 7.3, 7.4, and 8.1 evaluate the fixed PersistRank states.
- **Recovery, invalidation and cost:** the scenario manifest owns the
  initialization recovery record and checksum. Any relevant Stage 4/5 output,
  selector or coefficient change invalidates it; a changed curve, removal batch
  size, compatible schedule, maximum stage, checkpoint policy or report setting
  does not. New or pre-existing compatible schedule states import a missing
  initialization record from the scenario manifest, without overwriting a
  record already committed by that schedule. Each curve manifest records the
  shared fingerprint plus its selected coefficient identity, and `resume`
  starts from that curve's newest or explicitly selected compatible checkpoint
  without recomputing committed lookups. Initialization mode does not resolve
  or validate curve,
  schedule, stage-limit, logging or checkpoint controls because none belongs to
  the shared artifact. Only the current scenario-level initialization schema is
  supported. Native scoring plus repeated graph/spatial
  repairs make curve execution very expensive, while initialization's
  raster/index build is paid only once per compatible scenario.
- **Rmd order:** `setup`, `execute`, and `summary`.

### Stage 7.1 — PersistRank reports

- **Entry point and purpose:** `7.1_stage_meta_and_priority_surface.Rmd` calls
  `run_stage71()` to reconstruct, validate and visualize one completed Stage 6
  trajectory.
- **Inputs and configuration:** the validated shared initialization, one
  curve's positive-stage lookup CSVs and removal-event/order/rank surfaces,
  and the selected application/run/curve.
- **Ordered calculation:** build the Stage 6 state index; use the shared
  initialization for Stage 0 and the lookup CSVs for positive stages; select
  the curve coefficients from that initialization; evaluate the
  configured curve and all five uncertainty curves; preserve the complete Stage
  4 species denominator by assigning zero to absent spatial species; validate
  cell removal counts against all removal surfaces; and derive stage, priority,
  persistence and uncertainty summaries.
- **Returned:** stage metadata, species/assemblage trajectories, surface data and
  figure paths.
- **Durable outputs:** atomically replaced `stage_meta.csv` and four PersistRank
  figures. The metadata is human-readable report output, not a state handoff.
- **Downstream, recovery, invalidation and cost:** no scientific consumer reads
  the report files. They can be regenerated from Stage 6 without rerunning
  native kernels; curve/run scientific changes require Stage 6, while report
  implementation changes require only Stage 7.1. Cost is moderate and chiefly
  state evaluation/plotting.
- **Rmd order:** `setup`, `calculate`, displayed analysis chunks, and `summary`.

### Stage 7.2 — exact benchmark lookup library

- **Entry point and purpose:** `7.2_build_rank_lut.Rmd` calls `run_stage72()`.
  It reconstructs one Zonation method's ecologically valid patch/PU state at
  exactly the retained-cell counts required to compare it with selected Stage 6
  trajectories.
- **Inputs:** the shared Stage 6 initialization, completed Stage 6 state indexes
  for the selected curves, the immutable Stage 5 patch rasters and spatial
  metadata, and one method's external `rankmap.tif`.
- **Configuration:** `inspect`, `resume`, or method-scoped `restart`; selected
  optimization curves; benchmark method; and shared application/scenario/run
  selectors. Presentation settings do not enter lookup identity.
- **Ordered calculation:**

  1. Read the shared initialization once plus all selected Stage 6 state indexes
     and derive the sorted union of exact positive retained-cell counts.
     Validate their common initial domain, grid, cell ordering, species source
     and threshold contract.
  2. Validate the rank map against Stage 6 geometry and domain, fingerprint it,
     order domain cells from greatest to least rank with cell ID as the stable
     tie-breaker, and verify the initial count.
  3. Start from the Stage 5 state or a compatible recovery checkpoint. Traverse
     targets in decreasing retained-cell order, applying the next rank losses.
  4. After each rank loss, apply the same strict patch/PU threshold cascade,
     rook-fragmentation repair and dispersal-distance connectivity repair used
     to define valid PersistRank states.
  5. At an exact target, validate the canonical patch table and atomically write
     an immutable `retained_cells_<N>.rds`. Existing compatible exact files are
     never overwritten during resume.
  6. Checkpoint traversal every five unique targets and at the final target.
     When every required target exists, finalize the method library manifest and
     remove the no-longer-needed checkpoint.
- **Returned:** readiness, required/existing/missing counts, created/reused
  targets, phase/workload/timing summaries, checkpoint status and active log
  path.
- **Durable outputs:** one method library manifest, immutable exact lookup RDS
  files, a temporary reconstruction checkpoint while incomplete, and
  append-only `Logs/benchmark_<method>_reconstruction.log` operational data.
- **Downstream handoff:** Stages 7.3, 7.4, and 8.1 read exact lookup files;
  they never reopen a rank map or reconstruct spatial state once the library is
  complete.
- **Recovery, invalidation and cost:** `resume` skips committed exact files and
  continues from a compatible checkpoint. If the target sequence grows, it can
  restart traversal from the initial state while still reusing immutable exact
  files. Only standalone `restart` deletes the selected method library. A changed
  rank-map checksum, geometry, initial domain, ordering, Stage 5 source,
  threshold contract or run target plan invalidates that method. Reconstruction
  is very expensive once; reading a complete library is light.
- **Rmd order:** `setup`, `inspect`, `lookup-library`, and `artifact-summary`.

### Stage 7.3 — detailed comparison

- **Entry point and purpose:** `7.3_persist_cmp.Rmd` calls `run_stage73()` to
  compare one Stage 6 optimization curve with one completed Zonation method in
  a detailed seven-panel report.
- **Inputs and configuration:** one Stage 6 state index, its shared
  initialization, one complete exact benchmark library, curve/method,
  persistence threshold, exactly two focal species, point styles and
  assemblage central statistic.
- **Ordered calculation:** read paired PersistRank and benchmark states at identical
  retained-cell counts; aggregate patch areas to PUs; evaluate PU, species and
  fixed-denominator assemblage trajectories; calculate focal-species and
  assemblage summaries; compose the seven panels; return the requested plain
  tables; then release large trajectory objects.
- **Returned:** detailed statistics, plain tables, figure metadata and compact
  report summaries.
- **Durable output:** one replaceable detailed-comparison figure. Detailed
  trajectories and tables remain in memory only.
- **Downstream, recovery, invalidation and cost:** Stage 8.1 reuses the same
  evaluation/figure implementation. A complete library makes this a moderate
  report computation with no raster reconstruction. Curve, method or report
  settings rerun Stage 7.3 only unless the referenced scientific artifacts are
  themselves stale.
- **Rmd order:** `setup`, `inspect`, `calculate`, figure chunks and table chunks.

### Stage 7.4 — one-method cross-curve report

- **Entry point and purpose:** `7.4_cross_curve_priority_comparison.Rmd` calls
  `run_stage74()` to compare selected optimization curves with one completed
  Zonation method under all five evaluation functions.
- **Inputs:** every selected Stage 6 state index, the single shared
  initialization containing all five coefficient pairs, and one complete exact
  benchmark library.
- **Configuration:** `inspect` or `report`, optimization-curve subset,
  main curve, benchmark method, late-stage removal targets, persistence
  threshold and the shared application/run selectors.
- **Ordered calculation:** validate every selected curve and the complete
  method library; evaluate each optimization trajectory under all five
  persistence functions; pair PersistRank and benchmark states at exact retained
  counts; compute late-stage, sequence-integrated, retained-area, direction and
  cross-curve stability statistics; prepare plain tables; and release large
  trajectories before returning.
- **Returned:** target plan, paired statistics, plain tables, configuration
  summary, elapsed time and report-log path.
- **Durable outputs:** no scientific or report table cache. Active reports append
  only operational timings to `Logs/benchmark_<method>_report.log`.
- **Downstream, recovery, invalidation and cost:** Stage 8.1 uses the same
  cross-curve calculation. Results are regenerated without reconstruction;
  only lookup identity changes require Stage 7.2. Report computation is
  moderate and handles saved states rather than rank rasters.
- **Rmd order:** `setup`, `inspect`, `results`, `primary-results`,
  `robustness-matrices`, and `artifact-summary`.

### Stage 8 — selective orchestrator

- **Entry point and purpose:** `8_run_pipeline.Rmd` calls `run_stage8()` to run
  an exact, canonically ordered selection from Stages 2–7.1. Stage 1 is
  intentionally outside Stage 8 and must already be installed when Stage 2 is
  selected.
- **Inputs:** the selected shared configuration and reusable handoffs for
  omitted prerequisites. The configuration chunk also returns cheap resolved-
  path and standard-SDM-folder diagnostics before execution.
- **Parameters:** inspect/resume/restart, `stages`, optional `restart_from`,
  Stage 3 figure policy, and forwarded application/run settings.
- **Ordered calculation:** validate the complete configuration, resolve active
  parameters and input records only for the selected stages;
  put selected stages in canonical order; require core selections from Stages
  2–6 to be contiguous; validate every omitted prerequisite as an installed
  handoff; load the native lifecycle only when an active selected stage needs
  it; validate or build the shared Stage 6 initialization once before entering
  the curve loop; pass that checksum-validated recovery proof to every curve so
  their execution does not rescan all source rasters; execute each curve
  sequentially with independent checkpoint recovery; release its mutable result
  state before the next curve; and
  optionally regenerate Stage 7.1 reports.
- **Returned:** plan, per-stage status/results, artifact summary, and optional
  in-session Stage 5.3/7.1 presentation objects.
- **Durable outputs:** only the ordinary artifacts owned by invoked Stages
  2–7.1. Stage 8 creates no duplicate scientific format or orchestration cache.
- **Downstream, recovery, invalidation and cost:** Stage 8.1 consumes completed
  Stage 6 runs, not a Stage 8 result. Each invoked stage keeps its native
  recovery behavior; restart is limited to safe Stage 2, 4 and 6 replacement
  boundaries. Cost depends entirely on the literal selection.
- **Rmd order:** `configuration`, `run`, `plan`, `status`, `artifacts`, optional
  Stage 5.3/7.1 summaries, and `figures`.

### Stage 8.1 — four-method results

- **Entry point and purpose:** `8.1_report_results.Rmd` uses the public Stage
  8.1 operations to resume ABF/CAZ1/CAZ2/CAZMAX libraries and report them
  together.
- **Inputs:** completed Stage 6 curves, their one shared initialization, four
  rank maps for reconstruction, and exact libraries plus the model/species
  contract for reporting.
- **Parameters:** inspect/resume, curve subset/main curve, late targets,
  persistence threshold, focal species, and figure styles.
- **Ordered calculation:** inspection reads only small manifests; resume derives
  the shared union of exact targets and builds missing libraries method by
  method; reporting evaluates all PersistRank curves once, then loads and releases
  one method's exact lookups at a time; report tables/figures are regenerated
  from those fixed states; optional detailed ABF work runs separately so its
  large trajectories need not coexist with the four-method results.
- **Returned:** readiness, per-method resume status, cross-method statistics and
  plain tables, optional figure metadata, elapsed time and log paths.
- **Durable outputs:** reconstruction writes only Stage 7.2 exact lookups,
  manifests and temporary checkpoints. Reporting appends operational
  `Logs/stage81_report.log` and may replace the detailed ABF figure; it does not
  cache trajectories or tables.
- **Downstream, recovery, invalidation and cost:** non-destructive resume has no
  restart path. Fresh-session `setup → results` reads complete lookups without
  rank rasters or kernels. Missing lookup reconstruction is very expensive;
  completed-library reporting is moderate.
- **Rmd order:** `setup`, `inspect`, `lookup-libraries`, `results`, presentation
  table chunks, `main-figure`, and `artifact-summary`.

## 9. Execution recipes, selection, and restart boundaries

`stages` is literal: omitted prerequisites are validated and reused, never run
implicitly. Multiple selected core stages must be contiguous from Stage 2
through Stage 6. Optional Stages 5.1, 5.2, 5.3, and 7.1 may be omitted or run
later. For a calculation-only model-to-priority run, select Stages 2–6 and set
`stage3_figures: false`.

The configuration file is always checked as one complete, typo-free schema.
Stage 8 then resolves and validates filesystem inputs or runtime controls only
when their owning stage is selected. For example, a Stage 6-only resume uses
the stored application contract and does not fingerprint current Stage 4 input
paths, initialize the Stage 5 backend, or construct Stage 2 execution state.
Likewise, `priority_curves` is operationally relevant only for Stages 5.3, 6,
and 7.1. Existing manifests remain authoritative for omitted handoffs.

| Change | `restart_from` | Restart scope |
|---|---|---|
| Threshold, horizon, K-search or simulation identity | `stage_2` | persistence model and dependent threshold scenario |
| Species, SDMs, IUCN, traits, land cover, study area, patch rules | `stage_4` | named threshold scenario |
| Removal schedule, requested optimization curves, Stage 6 limits | `stage_6` | selected run/curve products |
| Optional products or reporting settings | null | select and rerun only their owning stage |

Stage 8 modes are `inspect`, `resume`, and deliberately scoped `restart`.
Stage 8.1 has no restart mode and cannot delete the four benchmark libraries.

The recipes below describe operational order, not expected wall time: Stages 2,
5, 6 and first-time benchmark reconstruction can each be long-running.

### Genuine from-scratch reproduction

1. Install the R, JAGS, compiler and geospatial requirements in Section 3.
2. Copy `config/example.yml` for a new application, fill every input and
   scientific setting, select it with `config_file` in the numbered Rmds, and
   run `inspect_project_config(read_project_config("config/your-file.yml"))`.
   For Madagascar, review `config/madagascar.yml` directly.
3. Set Stage 1's Rmd `mode` to `fit` and run `1_demographic_models.Rmd`.
   Do not proceed until all four fits satisfy the convergence contract and the
   demographic manifest is installed.
4. Run `2_persistence_simulation.Rmd` in `resume` mode. `resume` is also the
   correct first invocation because it creates new work when no compatible
   checkpoint exists and preserves recoverable work after interruption.
5. Run `3_persistence_curves.Rmd` in `fit` mode and validate the installed model.
6. Configure and run Stages 4 and 5. Use Stage 5 `inspect` before and after the
   active `resume`; inspect does not load rasters or initialize GRASS.
7. Run Stage 5.2, execute all four external Zonation methods, and put each
   `rankmap.tif` in its canonical method directory.
8. Run Stage 6 once with `mode: initialize`. After the shared initialization is
   installed, select each desired curve in turn and use `run` for a fresh state
   or `resume` after interruption. A compatible alternative removal schedule
   reuses the same scenario-level initialization.
9. Build each Stage 7.2 method library, or use Stage 8.1 `resume` to complete all
   four. Generate Stages 7.1, 7.3, 7.4, and 8.1 reports only after their
   scientific handoffs validate.

Stages 2–7.1 may instead be driven through Stage 8 after Stage 1 is complete.
A literal from-scratch Stage 8 selection is `2, 3, 4, 5, 5.1, 5.2, 5.3, 6,
7.1`, with optional stages omitted as needed; Stage 8 will not fit Stage 1 or run
external Zonation.

### Resume the bundled Madagascar analysis with Stage 8

1. Set `mode: inspect` and keep `stages` limited to the work of interest. Review
   the returned handoff, scenario, curve and checkpoint statuses.
2. Set `mode: resume` with the same literal selection. Omitted prerequisites are
   validated and reused. Selecting only Stage 6, for example, validates or
   builds the installed Stage 4/5-derived initialization once and ignores
   unselected Stage 4 input paths and Stage 5 backend settings. All requested
   curves then reuse that artifact.
3. Use `restart` only with an explicit `restart_from` from the table above. A
   normal interrupted computation should use `resume`, not restart.

### Prepare and run external Zonation

1. Run Stage 5.2 with `mode: inputs` to create the complete binary feature set.
2. Run the external Zonation program four times with the generated feature list
   and settings, selecting ABF, CAZ1, CAZ2 and CAZMAX.
3. Copy only the resulting rank raster required by R to the corresponding
   canonical `rankmap.tif` path in Section 10. Preserve its pixels, geometry and
   numeric ordering; do not resample it.
4. Run Stage 7.2 `inspect` for each method before reconstruction. Ancillary
   Zonation logs and summary tables can be retained for audit but are not inputs
   to any R stage.

### Build or resume one Stage 7.2 library

Set `rank_method` and `optimization_curves`, run
`setup → inspect → lookup-library → artifact-summary`, and use `mode: resume` for
both first-time and interrupted construction. Exact files already committed are
reused. Use method-scoped `restart` only after intentionally replacing its rank
map or another lookup-identity input.

### Complete all four libraries with Stage 8.1

Run `setup → inspect`, then `lookup-libraries` with `mode: resume`. Stage 8.1
derives the union of required exact targets, resumes ABF/CAZ1/CAZ2/CAZMAX one at
a time and reports created versus reused counts. After all libraries are
complete, a fresh session can run `setup → results` without rank rasters, native
kernels or a spatial runtime. Run `main-figure` separately when the detailed ABF
figure is wanted.

### Regenerate reports only

- Stage 7.1 needs the shared initialization plus a completed Stage 6 curve and
  lookup sequence.
- Stages 7.3 and 7.4 need the relevant complete exact library.
- Stage 8.1 `results` needs all four complete libraries; it evaluates PersistRank
  curves once and benchmark methods sequentially.

These paths read saved patch tables and coefficients. They do not need external
rank rasters, `terra`, GRASS, OpenMP kernels or reconstruction state. Their
tables are returned in memory and their figures are replaceable.

### Configuration policy

Scientific settings and input paths have exactly one human-edited source:
`config/madagascar.yml` or a copied configuration based on
`config/example.yml`. Its seven required sections are `demography`,
`persistence`, `application`, `inputs`, `spatial`, `priority`, and `reporting`;
`schema_version` must be 1. Every numbered Rmd exposes `config_file` plus only
its operation selectors. The precedence rule is therefore unambiguous:

1. the selected configuration supplies scientific, path, application, spatial,
   schedule, and default report values;
2. the Rmd supplies only allowed operational choices;
3. an installed manifest/checkpoint governs whether saved work is compatible
   for reuse or resume.

`project_config_params()` selects only the keys a stage is allowed to receive
and rejects an Rmd attempt to override a scientific key. Stage 8 receives all
settings it may need, including the same Stage 2/3 settings used by the
standalone Rmds. Scientific defaults therefore come from the selected shared
configuration rather than from Stage 8.
Use a standalone Rmd when its specialized operation is clearer, not because it
has different scientific defaults.

User-supplied relative input paths resolve from the repository root. Supported
path parameters also accept absolute paths and `~`. After editing an R module
interactively, either start a fresh R session or explicitly call
`project_source("R/file.R", reload = TRUE)` before rerunning dependent code.

## 10. External Zonation handoff

Stage 5.2 prepares four method directories: `ABF`, `CAZ1`, `CAZ2`, and
`CAZMAX`. These labels identify four separately configured external Zonation
runs. R treats the label as an opaque benchmark identity and does not infer or
change the Zonation method settings.

The binary feature rasters, feature list and generated settings file are the
external inputs common to these runs. The only Zonation result consumed by R is
`rankmap.tif`:

```text
Data/Applications/<application>/qe_<threshold>[_horizon_<years>yr]/Zonation/ABF/rankmap.tif
Data/Applications/<application>/qe_<threshold>[_horizon_<years>yr]/Zonation/CAZ1/rankmap.tif
Data/Applications/<application>/qe_<threshold>[_horizon_<years>yr]/Zonation/CAZ2/rankmap.tif
Data/Applications/<application>/qe_<threshold>[_horizon_<years>yr]/Zonation/CAZMAX/rankmap.tif
```

Stage 7.2 requires rank-map geometry, CRS, extent, resolution, NA/background
domain and initial retained-cell count to equal the Stage 5/6 source. It orders
rank values decreasingly with canonical cell ID as the deterministic tie-break.
The rank-map fingerprint/checksum is stored in the library identity; a changed
checksum or ordering stops before an incompatible exact lookup is written.
Other Zonation logs, curves and summary tables are external products, not R
workflow inputs.

### Benchmark comparison definitions

All comparisons pair PersistRank and Zonation states at the exact integer number
of retained cells specified by the PersistRank trajectory. Stage 7.2 reconstructs
the benchmark state at each required count; later stages do not approximate a
target by reading a nearby benchmark state.

The configured late-stage targets are percentages of cells removed. For each
target and optimization curve, the reporting code selects the available
PersistRank stage with the smallest absolute difference from that target and
reports the stage's actual removal percentage. It then calculates the mean
modeled persistence across the fixed species set, counts species below the
configured descriptive persistence threshold, and identifies species below the
threshold under only one ranking.

For a complete-sequence comparison, let
\(S_s^{\mathrm{PR}}(x)\) and \(S_s^{\mathrm{Z}}(x)\) be persistence for species
\(s\) under PersistRank and one benchmark at removal percentage \(x\). The
initial state is \(x=0\), and \(x_{\max}\) is the last evaluated stage. Values
between stages are treated as linear. The species-level signed mean difference
is

\[
\Delta_s = \frac{1}{x_{\max}}
\int_0^{x_{\max}}
\left[S_s^{\mathrm{PR}}(x)-S_s^{\mathrm{Z}}(x)\right]\,dx .
\]

The trapezoidal rule is applied to the evaluated stage points. Dividing by the
evaluated removal span leaves the result on the persistence scale. Positive
values favor PersistRank, and negative values favor the benchmark.

Cell counts do not guarantee that two rankings retain the same amount of
habitat relevant to a species. The retained-area comparison therefore treats
persistence as a function of the total area of population units that remain
above the species-specific unit threshold. For each ranking, the calculation:

1. adds the origin `(area, persistence) = (0, 0)`;
2. removes nonfinite observations and negative areas;
3. averages persistence when multiple states have the same retained area;
4. restricts both curves to their shared retained-area support;
5. combines both support endpoints and every area breakpoint from both curves;
6. linearly interpolates each curve at these combined knots without
   extrapolation; and
7. integrates the signed difference with the trapezoidal rule and divides by
   the shared area span.

Because the combined knots contain every breakpoint of both piecewise-linear
relationships, the trapezoidal calculation is the exact integral under the
stated linear interpolation. If the initial qualifying population-unit area is
\(A_{s,0}\), the bundled analysis has the common domain \([0,A_{s,0}]\) and the
normalized comparison is

\[
\bar\delta_s = \frac{1}{A_{s,0}}
\int_0^{A_{s,0}}
\left[S_s^{\mathrm{PR}}(A)-S_s^{\mathrm{Z}}(A)\right]\,dA .
\]

For both complete-sequence metrics, an absolute difference no greater than
`1e-9` is classified as nearly equal. The cross-curve sensitivity calculation
repeats the sequence comparison under all five evaluation functions while
holding each cell-removal order fixed. It reports species that favor PersistRank
under all five functions, favor the benchmark under all five, are nearly equal
under all five, or show at least one actual reversal in the favored ranking.
The first three counts and the reversal count are separate summaries; together
they are not intended to partition every species.

## 11. Stage 8.1 lookup, resume, and report workflow

`stage81_config(params)` reads small manifests and validates selectors.
`stage81_target_plan(config)` reads the shared Stage 6 initialization once and
derives consumers and unique exact targets from removal ledgers and patch lookup
filenames.
`inspect_stage81(config, plan)` returns one row per benchmark with required,
existing, missing, checkpoint, and compatibility fields.
`resume_stage81_lookups(config, plan)` builds only missing targets.
`build_stage81_results(config, plan)` regenerates trajectories/statistics/tables.
`build_stage81_figure(config, plan)` separately performs detailed ABF work.
There is deliberately no monolithic Run-All wrapper: inspection, spatial
lookup completion, report regeneration, and the optional detailed figure have
different costs and dependency requirements.

The reconstruction checkpoint is written every five unique targets and always
at the last target. A changed target sequence archives/replaces only the
checkpoint, starts traversal from the initial state, and skips writes for exact
targets already committed. A complete library has no checkpoint.

## 12. Individual-chunk execution

Every Rmd's first executable chunk sources one workflow module, resolves that
document's operational YAML through `project_rmd_params()`, and combines it
with `config_file` through `rmd_project_config()`. The chunk is named `setup` in
most report/spatial documents, `run` in Stages 2–4 and 5.1, and `configuration`
in Stage 8. Run that first chunk after editing the configuration or Rmd options,
or after starting a new R session.
Stage-specific calculation chunks then read required files directly;
presentation chunks name their prerequisite.

For Stage 8.1:

- readiness only: `setup → inspect`;
- missing lookup construction: `setup → lookup-libraries`;
- cheap report regeneration: `setup → results`;
- tables: run any presentation chunk after `results`;
- detailed figure only: `setup → main-figure` after ABF lookups are complete.

Stage 7.2 follows `setup → inspect → lookup-library`. Stage 7.4 follows
`setup → inspect → results → desired table chunks`. Stages 1–7.1 and 8
follow their visible chronological chunk order; no document needs objects from
another Rmd session.

## 13. Complete configuration and Rmd-option reference

The tables cover every key in the shared configuration and every operational
Rmd option. Unless explicitly described as an Rmd option, the key belongs in
the named configuration file, not in a numbered Rmd. “Identity” states what
persistent artifact the value helps identify. “Rerun” is the earliest required
boundary after a scientific change; `report` means no durable scientific
artifact is invalidated.

### Operation and Stage 1 parameters

| Parameter | Type / allowed values | Default | Meaning | Identity | Rerun |
|---|---|---|---|---|---|
| `schema_version` | integer; configuration metadata | 1 | selects the supported shared-configuration schema | parsing contract | update configuration only when schema changes |
| `config_file` | YAML path; Rmd option | `config/madagascar.yml` | selects the one complete scientific/input configuration used by that Rmd | configuration selection | rerun setup |
| `mode` | Stages 1/3 `reuse\|fit`; Stage 2 `inspect\|resume\|restart`; Stage 4 `figure\|build`; Stage 5 `inspect\|resume\|restart`; Stage 5.2 `inspect\|inputs\|feature_list`; Stage 6 `inspect\|initialize\|run\|resume`; Stage 7.2 `inspect\|resume\|restart`; Stages 7.3/7.4 `inspect\|report`; Stage 8 `inspect\|resume\|restart`; Stage 8.1 `inspect\|resume` | Stage 1 `fit`; 2 `inspect`; 3 `reuse`; 4 `figure`; 5/5.2/6/7.2/7.3/7.4 `inspect`; 8/8.1 `resume` | selects operation, never scientific identity | none | operation only |
| `verbose` | logical | `false` | detailed console progress in Stages 1, 2, 4, 5 | none | none |
| `bird_sigma_separate_intercepts` | any proper subset of `FruiNect`, `Invertebrate`, `Omnivore`, `PlantSeed`, `VertFishScav` | `VertFishScav` | bird residual-variation branches; unselected categories form `Other` | Stage 1 model | Stage 1+ |
| `seed` | non-negative integer | 123 | Stage 1 sampler seed | Stage 1 model | Stage 1+ |
| `n_chains` | integer ≥ 2 in `fit` and ≥ 1 in `reuse` | 3 | JAGS chain count | Stage 1 model | Stage 1+ |
| `n_adapt` | positive integer | 5000 | adaptation draws per chain | Stage 1 model | Stage 1+ |
| `n_iter` | positive integer divisible by `thin` | 5000 | sampled draws per chain | Stage 1 model | Stage 1+ |
| `thin` | positive integer dividing `n_iter` | 10 | posterior thinning interval | Stage 1 model | Stage 1+ |
| `coefficient_prior_mean` | finite number | 0 | regression-coefficient prior mean | Stage 1 model | Stage 1+ |
| `coefficient_prior_sd` | positive number | 100 | regression-coefficient prior SD | Stage 1 model | Stage 1+ |
| `residual_sd_min` | non-negative number below max | 0.01 | residual-SD lower prior bound | Stage 1 model | Stage 1+ |
| `residual_sd_max` | number above min | 5 | residual-SD upper prior bound | Stage 1 model | Stage 1+ |
| `rhat_max` | number ≥ 1 | 1.01 | convergence acceptance ceiling | validation contract | Stage 1 validation |
| `ess_min` | positive number | 400 | effective-sample-size acceptance floor | validation contract | Stage 1 validation |
| `trait_grid_points` | integer ≥ 5 | 31 | posterior prediction-grid resolution | Stage 1 grids | Stage 1+ |
| `mammal_mass_min_g` | positive number | 2 | mammal prediction-grid lower mass | Stage 1 grids | Stage 1+ |
| `mammal_mass_max_g` | number above min | 4,750,000 | mammal prediction-grid upper mass | Stage 1 grids | Stage 1+ |

The default Stage 1 schedule reproduces the Hilbers et al. (2017) MCMC
settings: three chains, 5,000 adaptation iterations per chain, 5,000 sampled
iterations per chain, thinning by 10, and no additional post-adaptation burn-in.
The combined posterior therefore contains 1,500 draws per fitted model.

### Stage 2–3 model parameters

| Parameter | Type / allowed values | Default | Meaning | Identity | Rerun |
|---|---|---|---|---|---|
| `demographic_uncertainty` | `coefficients_only\|posterior_predictive` | `coefficients_only` | posterior components propagated in simulation | Stage 2 model | Stage 2+ |
| `curves` | exact canonical five-curve vector | `q50,q16,q025,q84,q975` | Stage 2 output coefficient targets/order | Stage 2 model | Stage 2+ |
| `persistence_horizon_years` | positive integer | 1000 | persistence time horizon | persistence model/scenario/lookups | Stage 2+ |
| `quasi_extinction_abundance` | positive integer | 500 | abundance threshold for persistence | persistence model/scenario/lookups | Stage 2+ |
| `population_cap_factor` | number ≥ 1 | 1.1 | simulation cap relative to starting abundance | Stage 2 model | Stage 2+ |
| `growth_rate_buffer` | number in (0,1] | 0.8 | conservative K-search growth buffer | Stage 2 model | Stage 2+ |
| `n_posterior_draws` | positive integer | 1500 | posterior draws evaluated per trait | Stage 2 model | Stage 2+ |
| `replicates_per_draw` | positive integer | 2500 | stochastic trajectories per draw | Stage 2 model | Stage 2+ |
| `posterior_chunk_size` | positive integer ≤ `n_posterior_draws` | 25 | recovery/work unit size | checkpoint policy | resume Stage 2 |
| `base_seed` | non-negative integer | 123 | deterministic CRN seed | Stage 2 model | Stage 2+ |
| `anchor_probabilities` | at least three unique probabilities in (0,1) | `.025,.16,.5,.84,.975` | fitted anchor probability rows | Stage 2 points | Stage 2+ |
| `persistence_grid_min` | probability | .025 | regular reporting grid lower bound | Stage 2 points | Stage 2+ |
| `persistence_grid_max` | probability above min | .975 | regular reporting grid upper bound | Stage 2 points | Stage 2+ |
| `persistence_grid_step` | positive number | .025 | regular reporting grid interval | Stage 2 points | Stage 2+ |
| `additional_persistence_probabilities` | unique probabilities | .99 | extra reported targets | Stage 2 points | Stage 2+ |
| `k_search_start` | null or positive integer | null (next power of two above the threshold) | explicit initial K-search value | Stage 2 model | Stage 2+ |
| `k_search_relative_tolerance` | number in (0,1) | .01 | K-search stopping tolerance | Stage 2 model | Stage 2+ |
| `k_search_max` | positive number greater than effective `k_search_start` | 1e13 | hard K-search ceiling | Stage 2 model | Stage 2+ |
| `k_rounding` | `round\|ceiling` | `round` | conversion of numerical K to abundance | Stage 2 model | Stage 2+ |
| `threads` | null or positive integer | null | OpenMP threads; null uses runtime default | none if numerically equivalent | resume Stage 2 |
| `loess_span` | number in (0,1] | .5 | Stage 3 coefficient smoother span | Stage 3 model | Stage 3+ |

### Application and spatial parameters

| Parameter | Type / allowed values | Default | Meaning | Identity | Rerun |
|---|---|---|---|---|---|
| `application` | canonical lowercase name | `madagascar` | selects/owns an application hierarchy | application | new/select application |
| `minimum_patch_abundance` | positive integer | 10 | density-derived minimum patch abundance | application/spatial source/lookups | Stage 4+ |
| `cells_to_remove_per_iteration` | positive integer | 1000 | pruning batch size | removal run/lookups | Stage 6+ |
| `pruning_iterations_per_stage` | positive integer | 50 | iterations per committed ecological stage | removal run/lookups | Stage 6+ |
| `taxa` | `mammals\|birds\|both` | `both` | retained taxa and fitted coefficient coverage; downstream stages obtain it from the application manifest | application | Stage 4+ |
| `sdm` | `ppm\|rangebag\|ppm_rangebag` | `ppm_rangebag` | retained SDM method(s) | application | Stage 4+ |
| `sdm_index_file` | null or CSV path | null | portable explicit SDM inventory; null uses folders | application inputs | Stage 4+ |
| `iucn_mode` | `cache_only\|cache_or_query\|refresh` | `cache_or_query` | habitat-row resolution policy | resolved application inputs | Stage 4 if rows change |
| `iucn_pause_seconds` | non-negative number | 2 | delay between IUCN requests | none | none unless service output changes |
| `sdm_parent_dir` | directory path | `.` | root of standard four SDM folders | application inputs | Stage 4+ |
| `mammal_rmax_file` | tabular path | `Data/Raw/mammal_rmax.txt` | mammal growth calibration | demographic model | Stage 1+ |
| `sigma_file` | CSV path | `Data/Raw/sigma.csv` | mammal/bird environmental-variation calibration | demographic model | Stage 1+ |
| `bird_growth_file` | CSV path | `Data/Raw/bird_growth_niel_lebreton.csv` | bird growth calibration | demographic model | Stage 1+ |
| `mammal_traits_file` | tabular path | `Data/Raw/mammal_data.txt` | mammal trait source used by Stage 4 | application inputs | Stage 4+ |
| `bird_traits_file` | tabular path | `Data/Raw/bird_data.txt` | bird diet calibration in Stage 1 and bird traits in Stage 4 | demographic/application inputs | Stage 1+ if demographic branch changes; otherwise Stage 4+ |
| `bird_generation_lengths_file` | workbook path | `Data/Raw/cobi13486-sup-0004-tables4.xlsx` | bird generation lengths for Stages 1 and 4 | demographic/application inputs | Stage 1+ |
| `iucn_bird_synonyms_file` | CSV path | `Data/Raw/iucn_bird_synonyms.csv` | bird calibration-name resolution | demographic inputs | Stage 1+ |
| `species_list_file` | CSV path | `Data/Raw/simple_summary.csv` | selected species inventory | application inputs | Stage 4+ |
| `synonyms_file` | CSV path | `Data/Raw/synonyms.csv` | general synonym mappings | application inputs | Stage 4+ |
| `curated_synonyms_file` | CSV path | `Data/Raw/input_synonyms.csv` | explicit application synonym overrides | application inputs | Stage 4+ |
| `random_effects_file` | CSV path | `Data/Raw/random_effects.csv` | taxonomic random effects used in density | application inputs | Stage 4+ |
| `landcover_file` | raster path | `Data/Raw/esacci_2022_pfts.tif` | habitat/area raster source | spatial source | Stage 5+ |
| `zonation_settings_template` | text path | `Data/Raw/settings.z5.txt` | Stage 5.2 template whose feature-list pointer is replaced | external handoff | Stage 5.2 only |
| `study_area_mode` | `bounds\|vector\|full_raster` | `bounds` | spatial crop/mask representation | application/spatial source | Stage 4/5+ |
| `study_area_file` | null or vector path | null | vector mask in vector mode | application/spatial source | Stage 4/5+ |
| `study_area_layer` | null or layer name | null | optional vector layer | application/spatial source | Stage 4/5+ |
| `roi_xmin` | finite longitude | 43.18 | bounds-mode west edge | application/spatial source | Stage 4/5+ |
| `roi_xmax` | longitude above xmin | 50.56 | bounds-mode east edge | application/spatial source | Stage 4/5+ |
| `roi_ymin` | finite latitude | -25.64 | bounds-mode south edge | application/spatial source | Stage 4/5+ |
| `roi_ymax` | latitude above ymin | -11.89 | bounds-mode north edge | application/spatial source | Stage 4/5+ |
| `clump_backend` | `auto\|fasterRaster\|terra` | standalone Stages 5/5.1 `auto`; Stage 8 `fasterRaster` | patch clumping implementation recorded in metadata | spatial output provenance | Stage 5+ if result differs |
| `grass_dir` | null or directory path | null | machine-specific GRASS installation | none | rerun active Stage 5 operation |
| `process_figure_species` | scientific name or null | `Cryptoprocta ferox` | Stage 5.1 explanatory focal species | figure only | Stage 5.1 |

### Prioritization, benchmark, and report parameters

| Parameter | Type / allowed values | Default | Meaning | Identity | Rerun |
|---|---|---|---|---|---|
| `curve` | `q50\|q16\|q025\|q84\|q975` | `q50` | Stage 6 optimization or Stage 5.3/7.1 focus | selected curve run/report | Stage 6 for optimization; report otherwise |
| `max_stages` | positive integer or `Inf` | 100 in Stage 6; `Inf` in Stage 8 | execution stop limit | none | resume Stage 6 |
| `ecology_log_every_iterations` | null or positive integer | 10 | optional diagnostic cadence | none | none |
| `checkpoint_every_stages` | null or positive integer | 1 | Stage 6 recovery cadence | none | next checkpoint |
| `checkpoint_keep` | positive integer | 2 | recent Stage 6 checkpoints retained | none | none |
| `resume_stage` | null or completed stage integer | null | explicit Stage 6 recovery point | none | resume Stage 6 |
| `optimization_curves` | null or unique canonical subset | null | completed curves included in benchmark target/report plan | target plan, not lookup identity beyond exact counts | missing counts only/report |
| `rank_method` | `abf\|caz1\|caz2\|cazmax` | `abf` | selected benchmark library | method library | Stage 7.2 if missing |
| `main_curve` | selected canonical curve | `q50` | report focus and detailed figure curve | none | report only |
| `late_stage_removal_percentages` | increasing unique values in (0,100) | 90,95,99 | nearest-stage summary targets | none | report only |
| `persistence_threshold` | number in (0,1) | .75 | low-persistence classification cutoff | none | report only |
| `focal_species` | exactly two unique scientific names | `Pterocles personatus`; `Galidia elegans` | detailed panels/results | none | report/figure only |
| `assemblage_point_style` | `solid\|hollow\|none` | `solid` | assemblage panel markers | none | figure only |
| `focal_point_style` | `solid\|hollow\|none` | `none` | focal panels’ markers | none | figure only |
| `assemblage_central_stat` | `mean\|median` | `mean` | assemblage central line | none | figure only |
| `stages` | unique nonempty subset of `stage_2`, `stage_3`, `stage_4`, `stage_5`, `stage_5_1`, `stage_5_2`, `stage_5_3`, `stage_6`, `stage_7_1`, subject to core contiguity | `stage_2,stage_3,stage_4,stage_5,stage_5_1,stage_5_3,stage_6,stage_7_1` | exact Stage 8 operations; Stage 5.2 is not selected by default | none | selected stages only |
| `restart_from` | null or `stage_2\|stage_4\|stage_6` | null | safe replacement scope in restart mode | none | selected boundary |
| `stage3_figures` | logical | `true` | Stage 8 Stage 3 presentation policy | none | Stage 3 figures only |
| `write_figures` | logical | `true` | standalone Stage 3 presentation policy | none | Stage 3 figures only |
| `priority_curves` | unique canonical curve vector | all five | curves used by selected Stage 5.3, 6, or 7.1 work; ignored otherwise | set of curve runs | Stage 6 for new/changed curves |

`checkpoint_interval_targets` is deliberately not configurable. Benchmark recovery
uses the internal value five and always checkpoints the final target.

## 14. Complete artifact reference

| Stage / authoritative pattern | Contents and writer | Immediate reader / lifecycle |
|---|---|---|
| `Data/Models/Demography/Outputs/mammal_growth_posterior.csv` | Stage 1 posterior draws for the mammal growth allometry | Stage 2; durable |
| `Data/Models/Demography/Outputs/mammal_environmental_variation_posterior.csv` | Stage 1 posterior draws for mammal environmental variation | Stage 2; durable |
| `Data/Models/Demography/Outputs/bird_growth_posterior.csv` | Stage 1 posterior draws for bird growth | Stage 2; durable |
| `Data/Models/Demography/Outputs/bird_environmental_variation_posterior.csv` | Stage 1 posterior draws for bird variation/diet branches | Stage 2; durable |
| `Data/Models/Demography/Outputs/mammal_mass_grid.csv`, `bird_generation_length_grid.csv` | Stage 1 prediction/simulation trait grids | Stage 2; durable |
| `Data/Models/Demography/model_manifest.rds` and `Figures/Models/Demography/{demographic_calibration.png,SI/S1_demographic_posteriors.png}` | Stage 1 contract plus replaceable diagnostics | manifest is authoritative; figures are human-facing |
| `Data/Models/Persistence/<model>/Stage2/persistence_points_{mammals,birds}.csv` | Complete Stage 2 trait × curve × probability points | Stage 3; durable only when both validate |
| `Data/Models/Persistence/<model>/Stage2/persistence_points_{mammals,birds}.partial.csv` and `Stage2/Checkpoints/stage2_{mammals,birds}.rds` | Stage 2 completed-block and current-block recovery state | Stage 2 resume; temporary/superseded at completion |
| `Data/Models/Persistence/<model>/Stage3/persistence_curve_models.rds` and `Data/Models/Persistence/<model>/model_manifest.rds` | Stage 3 fitted five-function model and persistence-model authority | Stage 4 and lifecycle validation; durable |
| `Figures/Models/Persistence/<model>/{gompertz_wolff.png,SI/S3_gompertz_parameters.png,SI/S3_abundance_persistence.png}` | Stage 3 fit, coefficient and combined abundance/persistence diagnostics | users; replaceable |
| `Data/Applications/<application>/application_inputs.rds`, `application_manifest.rds` | Stage 4 resolved/fingerprinted application inputs and application identity | Stages 4–8.1; durable authority |
| `<scenario>/scenario_manifest.rds`, `Species/species_table.csv` | Stage 4 threshold scenario and canonical row-per-species contract | Stages 5–8.1; durable authority |
| `<scenario>/Spatial/Patches/<Genus_species>.tif` | Stage 5 canonical patch-ID raster for each retained spatial species | Stages 5.2, 6 and 7.2; durable |
| `<scenario>/Spatial/all_patch_lookup.rds`, `all_connectivity.rds`, `stage5_build_metadata.rds` | Stage 5 patch/PU areas, CSR connectivity and source/build contract | lookup: Stages 5.2–7.2; connectivity: Stages 6–7.2; durable authority |
| `<scenario>/Spatial/Checkpoints/manifest.rds` and `completed/<NNNN>_<Genus_species>/{result.rds,<Genus_species>.tif}` | Per-species Stage 5 recovery records while the final set is incomplete | Stage 5 resume; temporary/superseded |
| `<scenario>/Spatial/priority_initialization.rds` | Stage 6 immutable curve-/schedule-neutral species, cell, patch, PU, graph and all-five-coefficient starting state | all Stage 6 curves and Stages 7–8.1; one scenario-level durable authority |
| `<scenario>/Zonation/Patches_binary/<Genus_species>.tif`, `feature_list.txt`, `settings.z5.txt` | Stage 5.2 binary feature set and external configuration | Zonation; durable regenerable handoff |
| `<scenario>/Zonation/<METHOD>/rankmap.tif` | External ABF/CAZ1/CAZ2/CAZMAX cell ranking | Stage 7.2 only; external scientific input |
| `<run>/run_manifest.rds`, `run_state.rds` | application lifecycle record for schedule and completed curves, including shared-initialization fingerprint and each selected coefficient identity | Stage 6 and downstream handoff validation |
| `<curve>/Run/checkpoints/priority_checkpoint_stage_<NNNN>.rds` | complete Stage 6 state at a committed ecological stage | Stage 6 resume; recent files retained by policy |
| `<curve>/Run/patch_lookup_tables/stage_patch_lookup_stage_<NNNN>.csv` | positive-stage canonical patch/PU state | Stages 7.1–8.1; immutable once committed |
| `<curve>/Run/removal_events.csv`, `removal_order.tif`, `rankmap.tif` | removal ledger including terminal retained event, integer removal layer, and normalized PersistRank surface | Stages 7.1–8.1 and validation; durable |
| `<curve>/Run/Analysis/stage_meta.csv` | Stage 7.1 human-readable stage summary | humans only; durable but non-authoritative |
| `<run>/BenchmarkLookups/<METHOD>/lookup_manifest.rds` | Stage 7.2 method/rank/source/target library authority | Stages 7.2–8.1; durable |
| `<METHOD>/Lookups/retained_cells_<N>.rds` | immutable exact benchmark patch table | Stages 7.3–8.1; durable and shared across curves at equal counts |
| `<METHOD>/checkpoint.rds` and `checkpoint.rds.bak` | mutable Stage 7.2 traversal recovery and its transaction backup | reconstruction resume; absent after clean completion |
| `Figures/Applications/<application>/<scenario>/<run>/{Spatial,Priority,Comparisons}/*.png` | Stage 4/5.1/5.3/7.1/7.3/8.1 presentation products | users; atomically replaceable |
| Returned PU/species trajectories, statistics, target plans, formatted tables and matrices | Stage 5.3 and Stages 7.3–8.1 report functions | current R session only; deliberately not cached |

Each exact lookup contains only `schema`, `retained_cells`, and a canonical
`patch_table`. Method and rank-map identity live once in `lookup_manifest.rds`.

### Runtime diagnostic logs

Runtime logs contain only structured performance diagnostics for expensive
operations. Each active invocation appends a `runtime_session_start` line,
timed events, and a `runtime_session_end` line. Every line is appended and the
file is closed before the identical line is printed to the console. Failure to
append is fatal, so expensive work never continues with an incomplete selected
diagnostic. A start line without a matching end line indicates process
termination or interruption. Logs never participate in scientific manifests,
checksums, compatibility, checkpoint identity, or restart decisions; they may
be archived or deleted safely and are never truncated by the workflow.

| Canonical path | Creating operation | Event families |
|---|---|---|
| `Data/Models/Persistence/<model>/Logs/stage2.log` | active Stage 2 simulation, including the Stage 8 Stage 2 route | session boundaries; setup, curve, and trait timing |
| `Data/Applications/<application>/<scenario>/Logs/stage5.log` | active Stage 5 resume/restart | session boundaries; setup and attempted noncached-species timing |
| `Data/Applications/<application>/<scenario>/Logs/stage6_initialization.log` | active shared Stage 6 initialization, standalone or through Stage 8 | session boundaries; source reads, index/graph construction, validation and publication timings |
| `Data/Applications/<application>/<scenario>/Runs/<run>/Logs/stage6_<curve>.log` | Stage 6 curve run/resume, including Stage 8 | session boundaries; phase/stage/checkpoint timings, workloads, hotspots, and sampled ecological diagnostics |
| `Data/Applications/<application>/<scenario>/Runs/<run>/Logs/benchmark_<method>_reconstruction.log` | missing-target reconstruction through Stage 7.2 or Stage 8.1 | session boundaries; setup/resume/checkpoints, union targets, phase/progress timings, fragmentation/predicate/geometry workloads |
| `Data/Applications/<application>/<scenario>/Runs/<run>/Logs/benchmark_<method>_report.log` | active standalone Stage 7.4 report | session boundaries; PersistRank curves, exact-lookup evaluation, and report timing |
| `Data/Applications/<application>/<scenario>/Runs/<run>/Logs/stage81_report.log` | active Stage 8.1 multi-method report | session boundaries; PersistRank curves, one lookup evaluation and method timing per benchmark, and total report timing |

Here `<method>` is one of `abf`, `caz1`, `caz2`, or `cazmax`. Inspection,
artifact summaries, target planning, figures, and a complete no-op benchmark
resume do not open a runtime-log session or modify an existing log. Ordinary
configuration, progress, table, warning, and verbose scientific-debug output
remains console-only. Logging has no configuration/Rmd setting or runtime mode.

## 15. Complete source-file map

This inventory covers all 16 numbered Rmd entry points, all 134 R
modules, and both native C++ sources. Repeated species rasters and generated
exact-target files are data artifacts and are documented by pattern instead.

### Dependency loading

The dependency direction is numbered Rmd → workflow → stage modules → shared
application context/state → contracts, paths, and utilities. Each workflow
declares its complete ordered dependency set through `project_source()` from
`R/project_loader.R`. The loader uses base R, loads definitions into the global
environment for Rmd compatibility, and parses each normalized path at most once
per session. After an interactive edit, start a fresh R session or call
`project_source("R/file.R", reload = TRUE)` explicitly.

Ordinary production modules are definition-only: they never source other
modules or choose dependencies from the current workspace. Stage 6 is the
intentional exception: `load_stage6_modules()` keeps its spatial runtime and
compiled kernels deferred. Benchmark reconstruction likewise declares its
spatial modules in one lazy loader, so inspection and completed-report paths do
not load `terra`, `sf`, or native kernels.

Canonical storage is layered deliberately: `application_storage.R` constructs
paths without touching the filesystem, `storage_artifacts.R` owns relocatable
fingerprints and atomic metadata writes, and `model_lifecycle.R` owns reusable
demography/persistence validation and manifests. This keeps report-only loading
independent of model lifecycle operations.

Every numbered workflow begins with the same ordered infrastructure layers:
foundational utilities, artifact helpers, Rmd adaptation, runtime logging, and
transactions. Merely sourcing these files attaches no package and performs no
filesystem work. Optional GRASS/fasterRaster discovery is isolated in
`R/patch_runtime.R`; only active Stage 5 and Stage 5.1 paths load its definitions,
so Stage 8 inspection remains free of patch-runtime initialization.

### Configuration flow

`R/project_config.R` reads the selected complete configuration, rejects unknown,
missing, or duplicate keys, selects the keys owned by one stage, and overlays
only that Rmd's allow-listed operational options. `config/madagascar.yml`
contains the bundled values; `config/example.yml` is the complete copy-and-edit
template. `application_run_identity()` then derives the shared application,
threshold, model, scenario, and removal-schedule paths once. After an
application is created, its manifests are authoritative for taxa, SDM
selection, minimum patch abundance, source fingerprints, and completed
optimization curves. Stage configs remain ordinary lists and validate only
their own active operation. Configuration inspection never opens raster cells,
loads spatial/native packages, contacts IUCN, or writes files.

- `config/madagascar.yml`: authoritative bundled scientific settings and input
  paths for the Madagascar application.
- `config/example.yml`: complete copy-and-edit configuration for a new
  application; placeholders are intentionally not runnable until replaced.

Stages 1–3 follow the same internal boundary: configuration is pure, preflight
and artifact handoff live in an operational input/artifact module, and expensive
work begins only through the existing driver. Stage 2's native runtime module is
safe to source and compiles only when an active simulation explicitly requests
it.

### Execution lifecycle

Stage workflows own operational preflight and runtime loading. In particular,
Stage 5 performs preflight before attaching spatial packages, Stage 5.2 keeps
inspection and feature-list mode free of `terra`, and Stage 5.3 validates its
inputs before reading artifacts or creating figures. Numbered Rmds therefore
remain thin adapters: they construct context/config, invoke one public workflow
operation, and display its returned result.

`R/application_lifecycle.R` is the sole translator between stage results and
durable recovery boundaries. The same translators are used by standalone Rmds
and Stage 8, so output sets, upstream artifacts, summaries, scenario completion,
and curve completion have one owner. `R/stage8_plan.R` contains only read-only
handoff planning; it does not load Stage 4–7 science, spatial packages, or native
kernels.

### Numbered entry points

Each Rmd sources exactly the one workflow shown here. The workflow then asks the
project loader for its dependencies in fixed order.

| Rmd | Single workflow source | Responsibility |
|---|---|---|
| `1_demographic_models.Rmd` | `R/stage1_workflow.R` | Stage 1 model driver and summaries |
| `2_persistence_simulation.Rmd` | `R/stage2_workflow.R` | Stage 2 inspect/simulate driver |
| `3_persistence_curves.Rmd` | `R/stage3_workflow.R` | Stage 3 fit/reuse and figures |
| `4_build_species_table.Rmd` | `R/stage4_workflow.R` | Stage 4 build/figure driver |
| `5_build_patches_and_connectivity.Rmd` | `R/stage5_workflow.R` | Stage 5 inspect/resume/restart |
| `5.1_single_species_aoh_patches_pu_process.Rmd` | `R/stage51_workflow.R` | one-species explanatory figure |
| `5.2_build_binary_patch_rasters_zonation_feature_list.Rmd` | `R/stage52_workflow.R` | Zonation feature preparation |
| `5.3_initial_species_persistence.Rmd` | `R/stage53_workflow.R` | baseline persistence report |
| `6_spatial_prioritization_pipeline.Rmd` | `R/stage6_workflow.R` | shared initialization and one-curve Stage 6 execution façade |
| `7.1_stage_meta_and_priority_surface.Rmd` | `R/stage71_workflow.R` | PersistRank report driver |
| `7.2_build_rank_lut.Rmd` | `R/stage72_workflow.R` | one-method lookup construction |
| `7.3_persist_cmp.Rmd` | `R/stage73_workflow.R` | detailed comparison report |
| `7.4_cross_curve_priority_comparison.Rmd` | `R/stage74_workflow.R` | report-only cross-curve comparison |
| `8_run_pipeline.Rmd` | `R/stage8_workflow.R` | exact selected Stage 2–7.1 orchestration |
| `8.1_report_results.Rmd` | `R/stage81_workflow.R` | four-method lookup/report workflow |

### Shared contracts, storage, and presentation

- `R/project_loader.R`: deterministic, source-once workflow dependency loader.
- `R/project_utils.R`: foundational assertions, scalar/parameter validation,
  optional paths, package checks, species identifiers, and console logging.
- `R/project_artifacts.R`: file/directory checks, writable output preparation,
  lightweight artifact status, and cached-artifact freshness validation.
- `R/project_rmd.R`: Knit/Run-All parameter resolution and chunk-order guards.
- `R/project_runtime_log.R`: state-preserving append-before-print operational
  runtime logging.
- `R/project_transactions.R`: atomic file-set and directory-replacement
  transactions with rollback.
- `R/project_paths.R`: repository-level paths and normalized public input paths.
- `R/project_config.R`: versioned shared YAML schema, stage-specific key
  selection, operation merge rules, and lightweight input/SDM-folder inspection.
- `R/analysis_contract.R`: five-curve, abundance, schema, threshold, area, and
  benchmark enumerations.
- `R/demographic_contract.R`: canonical mammal/bird demographic branch contract.
- `R/application_storage.R`: deterministic model/application/scenario/run/curve,
  benchmark, and figure path construction.
- `R/storage_artifacts.R`: relocatable file fingerprints, artifact records,
  comparisons, and atomic metadata persistence.
- `R/model_lifecycle.R`: demographic- and persistence-model validation,
  recovery state, contracts, and manifest installation.
- `R/application_state.R`: application boundary recovery records and atomic
  manifest validation/writes.
- `R/application_lifecycle.R`: shared result-to-boundary translators, atomic
  boundary commits, scenario finalization, and priority-curve completion.
- `R/application_context.R`: shared application/scenario construction, input
  contracts, completed-run manifest validation, and handoff validation used by
  standalone Rmds, Stage 8, and reports without depending on Stage 8.
- `R/application_rmd.R`: thin standalone-Rmd adapter for parameter conversion;
  lifecycle decisions remain in `R/application_lifecycle.R`.
- `R/report_config.R`: shared lightweight focal-species, threshold, point-style,
  and assemblage-summary presentation settings.
- `R/figure_utils.R`: manuscript-neutral themes, labels, and atomic figure saving.
- `R/bird_generation_lengths.R`: bird scientific-name normalization and
  generation-length matching.

### Stages 1–3

- `R/stage1_config.R`: pure Stage 1 parameter, path, and model validation.
- `R/stage1_artifacts.R`: Stage 1 posterior contracts, validated reads, and
  operational preflight.
- `R/stage1_bird_matching.R`: deterministic EltonTraits/IUCN matching for the
  bird environmental-variation calibration data.
- `R/stage1_data.R`: calibration input reading, joins, exclusions, and final
  demographic-table validation.
- `R/stage1_models.R`: JAGS model fitting, diagnostics, and posterior summaries.
- `R/stage1_figures.R`: Stage 1 allometry figures.
- `R/stage1_posterior_figures.R`: Stage 1 posterior-coefficient figures.
- `R/stage1_outputs.R`: fit transactions plus reuse-mode grid and figure writes.
- `R/stage1_workflow.R`: deterministic Stage 1 dependency loader/orchestrator.
- `R/stage2_config.R`: pure Stage 2 identity, numerical settings, and paths.
- `R/stage2_runtime.R`: deferred OpenMP/C++ loading, native contracts, CRN
  contexts, and runtime interface checks.
- `R/stage2_inputs.R`: operational preflight, Stage 1 posterior handoff, and
  deterministic simulation-input preparation.
- `R/stage2_checkpoints.R`: compact Stage 2 simulation checkpoint contract.
- `R/stage2_simulation.R`: CRN simulation, K search, validation, and output commit.
- `R/stage2_workflow.R`: deterministic Stage 2 loader/orchestrator;
  `run_stage2_model()` also publishes the shared recovery boundary.
- `R/stage3_config.R`: pure Stage 3 model/figure identity and LOESS settings.
- `R/stage3_inputs.R`: operational preflight, package loading, Stage 2 point
  handoff, and provenance validation.
- `R/stage3_gompertz.R`: exact shifted-Gompertz fitting and fit diagnostics.
- `R/stage3_smoothing.R`: positive LOESS smoothing and the single reusable
  prediction path invoked only when figures require predictions.
- `R/stage3_models.R`: saved-model construction, provenance, and reuse
  validation.
- `R/stage3_figures.R`: shared formatting and Wolff/Gompertz calibration figure.
- `R/stage3_parameter_figures.R`: Gompertz parameter/LOESS presentation.
- `R/stage3_persistence_figures.R`: abundance/persistence heatmap panels and their combined manuscript figure.
- `R/stage3_outputs.R`: validated Stage 3 model/figure transaction.
- `R/stage3_workflow.R`: deterministic Stage 3 loader/orchestrator;
  `run_stage3_model()` publishes state and the reusable model manifest.
- `src/simulate_persist_probs_cpp.cpp`: OpenMP persistence simulator used only
  when Stage 2 simulation begins.

### Stages 4–5.3

- `R/species_table_config.R`: Stage 4 parameters, SDM input mode, paths, and preflight.
- `R/species_table_inputs.R`: species, trait, model, and SDM-index ingestion.
- `R/species_name_candidates.R`: scientific-name candidates, priorities,
  ambiguity handling, and deterministic selection.
- `R/species_name_diagnostics.R`: name-resolution and exclusion reporting.
- `R/species_table_names.R`: high-level raster, trait, and generation-length
  resolution coordination.
- `R/species_trait_covariates.R`: diet, density, dispersal, and published trait
  coefficients.
- `R/species_gompertz_parameters.R`: Stage 3 model-contract validation and
  Gompertz coefficient attachment.
- `R/species_table_traits.R`: high-level trait-enrichment coordination.
- `R/species_table_iucn.R`: IUCN habitat normalization and shared cache transaction.
- `R/species_table_build.R`: canonical species-table validation and atomic writes.
- `R/species_table_figures.R`: Stage 4 area-curve report.
- `R/stage4_workflow.R`: deterministic Stage 4 build/figure orchestration.
- `R/patch_config.R`: Stage 5 parameters, study area, backend, and paths.
- `R/patch_runtime.R`: lazily loaded GRASS discovery, environment setup, and
  fasterRaster initialization used only by active Stage 5/5.1 work.
- `R/patch_contract.R`: Stage 5 patch/connectivity schemas and validation.
- `R/patch_species_inputs.R`: Stage 4 species/SDM selection for patch processing.
- `R/patch_habitat.R`: raster context, study-area crop/mask, and habitat masks.
- `R/patch_processing.R`: patch/PU/connectivity scientific processing.
- `R/patch_process_figure.R`: one-species explanatory map construction.
- `R/stage5_checkpoints.R`: validated per-species Stage 5 recovery files.
- `R/stage5_workflow.R`: full Stage 5 inspect/resume/restart orchestration and
  compact post-publication summaries.
- `R/stage51_workflow.R`: single-read Stage 5.1 explanatory workflow with
  spatial initialization deferred until the focal row is confirmed.
- `R/stage52_config.R`: Stage 5.2 operation, paths, and read-only preflight.
- `R/stage52_workflow.R`: deterministic Stage 5.2 Zonation-input workflow;
  validates the lookup/raster contract without loading unused connectivity.
- `R/stage53_config.R`: Stage 5.3 selectors, paths, and preflight.
- `R/stage53_analysis.R`: Stage 5.3 species/PU persistence input and calculation.
- `R/stage53_figures.R`: Stage 5.3 plots and atomic SI-figure persistence.
- `R/stage53_workflow.R`: preflight-first Stage 5.3 sequencing and result assembly.

### Stage 6 spatial prioritization

- `R/stage6_config.R`: shared initialization path plus curve/run parameters and
  canonical result paths.
- `R/stage6_workflow.R`: inspect, shared initialization, deferred runtime
  loading, coefficient selection and curve-run dispatch.
- `R/priority_run_config.R`: canonical selectors, shared initialization
  metadata, all-five-coefficient contract and compact per-curve identity.
- `R/priority_inputs.R`: immutable scenario-level prioritization-initialization
  construction, source fingerprinting and validation.
- `R/priority_outputs.R`: stage lookup, ledger, raster, and final transactions.
- `R/priority_checkpoints.R`: Stage 6 checkpoint validation/retention.
- `R/priority_pipeline.R`: high-level pruning-stage loop and compact completion
  summary; large final mutable state is released after durable publication.
- `R/priority_runtime_kernels.R`: source-aware deferred Rcpp kernel loading.
- `R/priority_logging.R`: structured Stage 6 phase, workload, hotspot, and
  checkpoint profiling events written through the shared runtime logger.
- `R/priority_ecology_logging.R`: optional ecological diagnostics from compact state.
- `R/priority_graph_store.R`: mutable/lazy graph storage abstraction.
- `R/priority_csr_graph.R`: deterministic CSR component algorithms.
- `R/priority_pruning_frontier.R`: live frontier maintenance.
- `R/priority_pruning_scoring.R`: patch/PU score caches and candidate scoring.
- `R/priority_pruning_graph.R`: pruning-triggered PU graph repair.
- `R/priority_pruning_iteration.R`: one pruning iteration and its profiling.
- `R/priority_pruning_stage.R`: repeated iterations forming one ecological stage.
- `R/priority_fragmentation_patches.R`: incremental patch/index repair.
- `R/priority_fragmentation_graph.R`: provisional PU graph construction.
- `R/priority_fragmentation_stage.R`: fragmentation transaction for changed species.
- `R/priority_distance_geometry.R`: distance predicate and geometry preparation.
- `R/priority_distance_graph.R`: distance-invalid edge filtering.
- `R/priority_distance_stage.R`: connectivity repair and state commit.
- `src/priority_stage6_kernels.cpp`: compiled scoring/component kernels shared by
  Stages 6 and benchmark reconstruction.

### Stages 7–8.1 and shared benchmarks

- `R/stage7_config.R`: standalone Stage 7.1 configuration.
- `R/stage7_contracts.R`: sole authoritative Stage 6 state-index reader.
- `R/stage7_persistence.R`: PU and independent-PU species persistence science.
- `R/stage71_report.R`: five prepared Stage 7.1 report tables.
- `R/stage71_figures.R`: Stage 7.1 persistence figure builders and themes.
- `R/stage7_removal_order_map.R`: removal-order figure preparation.
- `R/stage71_workflow.R`: Stage 7.1 evaluation, metadata export, and figures.
- `R/stage7_rank_inputs.R`: rank-raster inputs and initial reconstruction state.
- `R/stage7_rank_incremental.R`: benchmark incremental spatial transition adapter.
- `R/benchmark_config.R`: shared application context, reconstruction config,
  report settings, and canonical benchmark paths.
- `R/benchmark_targets.R`: shared consumer/unique-target plan.
- `R/benchmark_library.R`: exact-file/manifest contracts, light/strict validation,
  atomic lookup writes, inspection, finalization, and one-method restart.
- `R/benchmark_reconstruction.R`: decreasing-target traversal and checkpointing.
- `R/benchmark_evaluation.R`: patch → PU area → persistence evaluation.
- `R/benchmark_statistics.R`: late-stage and sequence statistics, exact
  combined-knot retained-area integration, direction classifications, and
  cross-function stability summaries.
- `R/benchmark_tables.R`: format-neutral presentation tables.
- `R/benchmark_workflow.R`: stable shared inspect/resume/evaluate public operations.
- `R/stage72_config.R`: Stage 7.2 one-method selectors and canonical paths built
  from the shared benchmark configuration.
- `R/stage72_workflow.R`: Stage 7.2 inspection/resume dispatch and the sole
  destructive one-method restart boundary.
- `R/stage73_config.R`: detailed one-curve/method report settings and paths.
- `R/stage73_comparison.R`: detailed PU/species trajectory construction.
- `R/stage73_statistics.R`: detailed checkpoint, sequence, exact retained-area,
  zero-crossing, redundancy, uncertainty, and focal-species statistics.
- `R/stage73_figure_data.R`: shared stage selection, uncertainty bands, and
  single-pass focal/assemblage figure-data preparation.
- `R/stage73_focal_figures.R`: focal panels, headers, strips, legends, and grid.
- `R/stage73_assemblage_figures.R`: assemblage panels, legend, and section.
- `R/stage73_report.R`: concise plain report-table preparation.
- `R/stage73_figures.R`: final combined figure dimensions and composition.
- `R/stage73_report_pipeline.R`: single-pass in-memory Stage 7.3 report assembly.
- `R/stage73_workflow.R`: report-only detailed orchestration and memory release.
- `R/stage74_config.R`: Stage 7.4 cross-curve selectors and report-only
  configuration derived from the shared benchmark context.
- `R/stage74_workflow.R`: Stage 7.4 input validation, evaluation, plain-table
  assembly, runtime logging and memory release.
- `R/stage8_config.R`: selection-aware Stage 8 parameter validation and lightweight
  contexts; omitted-stage settings are not resolved or validated.
- `R/stage8_plan.R`: lightweight, read-only boundary and handoff planning.
- `R/stage8_workflow.R`: exact selected Stage 2–7.1 sequencing using shared
  lifecycle translators.
- `R/stage81_config.R`: four-method report/reconstruction configuration.
- `R/stage81_workflow.R`: independent Stage 8.1 inspect, lookup-resume, report,
  and optional detailed-figure operations.

## 16. Resume, restart, checkpoint, and interruption behavior

- Application recovery records fingerprint both their published outputs and
  the upstream artifacts that justified the boundary. Reuse requires both sets
  to remain unchanged. Stage 8 planning stays lightweight; active execution
  performs this checksum validation immediately before a recorded boundary is
  actually reused. Standalone Rmd and Stage 8 execution then pass that proof to
  Stage 6, avoiding duplicate hashing and per-curve source-raster timestamp
  scans. A direct programmatic `run_stage6()` call, which has no lifecycle
  proof, retains the conservative source-freshness scan.
- Stage 2 checkpoints posterior simulation chunks.
- Stage 5 checkpoints validated per-species spatial products.
- Stage 6 checkpoints completed pruning stages and keeps the configured recent
  count. Its shared initialization is owned by the scenario and is never
  deleted by a curve Run restart or duplicated into a curve directory.
- Stage 7.2 checkpoints the mutable incremental benchmark state every five
  unique targets and at the final target. Exact lookup writes are atomic and
  never overwrite compatible files.
- Stage 8.1 `resume` loops over four Stage 7.2 libraries; it has no destructive
  operation. A second resume performs zero rank-raster reads and kernel loads.
- Stage 7.4/results-only execution performs zero checkpoint and lookup writes.
- Runtime logs append independently across starts, resumes, and reports. Restart
  operations do not delete them, and interrupted sessions remain visible as a
  start record without a matching end record.
- A changed rank map or scientific contract is an error. Deliberately use
  standalone Stage 7.2 `restart` for that one method.

Only the current schemas are supported: application recovery state is
schema 3, Stage 5 checkpoints/manifests are schema 2, the shared Stage 6
initialization is schema 5, and Stage 6 checkpoints are schema 4. Artifacts with
other schemas are rejected with a regeneration instruction; the workflow does
not reinterpret incompatible layouts.

## 17. Configuration examples

Start by copying the complete template; the fragments below show only the
sections to change and are not valid replacements for the full file.

```r
file.copy("config/example.yml", "config/africa.yml")
```

Different application identity and inputs in `config/africa.yml`:

```yaml
application:
  application: africa
  minimum_patch_abundance: 10
  taxa: both
  sdm: ppm_rangebag
inputs:
  species_list_file: inputs/africa_species.csv
  sdm_index_file: inputs/africa_sdm_index.csv
  # retain and edit every other required input key from config/example.yml
```

Set `config_file: "config/africa.yml"` in each numbered Rmd used for this
application.

Different threshold/horizon: set `quasi_extinction_abundance: 50` and
`persistence_horizon_years: 500`, then start Stage 8 at `stage_2`.

Different removal schedule: set `cells_to_remove_per_iteration: 500` and
`pruning_iterations_per_stage: 100`; this creates a separate run folder.

Curve subset is an operational Rmd choice, not a shared scientific setting. For
Stage 8.1, edit its Rmd front matter:

```yaml
optimization_curves:
  value: ["q16", "q50", "q84"]
main_curve: "q50"
```

Vector study area: use `study_area_mode: "vector"`, set
`study_area_file`, and optionally `study_area_layer`. Bounds mode uses the four
ROI coordinates. A portable SDM collection should use `sdm_index_file`; paths
inside that CSV may be relative to the CSV location.

## 18. Figure and table inventory

- Stage 1: demographic fit/diagnostic figures.
- Stage 3: fitted persistence-function figures.
- Stage 4: the optional `area_curve.png`; build-mode species/trait/parameter
  summaries are returned as ordinary tables rather than additional figures.
- Stage 5.1: optional single-species process figure.
- Stage 5.3: initial persistence figure/tables.
- Stage 7.1: removal order, PersistRank persistence, curve uncertainty, species
  persistence.
- Stage 7.3/8.1 main figure: seven-panel ABF/main-curve comparison plus unique
  focal-species tables.
- Stage 7.4/8.1: late-stage, sequence, retained-area, robustness, direction,
  and stability tables for each benchmark.

Tables are ordinary data frames/data tables with readable labels and scientific
values. No table contains HTML markup.

## 19. Reproducibility and atomic writes

Manifests record schemas, scientific contracts, checksums, and artifact paths.
The shared storage layer records a project-relative path for files inside the
repository and readers use that path after relocation. Provenance tables can
also retain the absolute source path that was active when an input-derived
artifact was built. Exact lookup manifests record method, rank-map path and
checksum, application, threshold, removal run, Stage 5 source fingerprint, and
reconstruction schema. New files are serialized to a same-directory staging
path, validated, and promoted transactionally. Compatibility uses content and
contracts, never timestamps.

## 20. Server execution and output transfer

Run Stage 8 through Stage 6 on the server, then transfer the entire relevant
`Data/Models`, `Data/Applications/<application>`, and matching
`Figures/Applications/<application>` trees. Canonical artifact records use
project-relative resolution when available, so completed states can move with
the project hierarchy. Do not transfer `Data/Cache/Rcpp`; it is platform-specific
and regenerated locally. Preserve exact filenames and method folders.
The transferred `Logs` directories are optional operational history rather than
scientific input. Keeping them preserves server bottleneck diagnostics; omitting
them does not invalidate any artifact.

### Paths after cloning or moving the project

Scientific manifests and completed report inputs use relocatable fingerprints
and project-relative artifact records. Completed Stage 7–8.1 reporting can
therefore use the bundled scientific states after the repository is moved.
Two files that hand data to machine-specific software contexts contain absolute
paths by design:

- `Species/species_table.csv` records the source distribution-raster path used
  when Stage 4 builds the application; and
- `Zonation/feature_list.txt` and the generated `settings.z5.txt` identify the
  local feature-list and raster paths supplied to Zonation.

No path repair is needed to regenerate reports from completed Stage 5, Stage 6,
and benchmark artifacts. Before rebuilding Stage 5 from distribution rasters on
a different machine, run Stage 4 in `build` mode with the local configured input
directories so the canonical species table records the current raster paths;
then follow the reported lifecycle/restart boundary for Stage 5. Before running
Zonation, run Stage 5.2 in `feature_list` mode to rewrite only the feature list
and generated settings file from the binary rasters currently on disk. That
operation does not open or modify the rasters. Do not manually edit scientific
manifests or their fingerprints.

## 21. Troubleshooting

- **No completed curves:** finish Stage 6 and ensure the run manifest marks the
  curve complete.
- **Missing Stage 6 state:** check the scenario-level
  `Spatial/priority_initialization.rds`, `removal_events.csv`, and contiguous
  `stage_patch_lookup_stage_####.csv` files. Re-run `initialize` only when the
  shared artifact is absent, corrupt, or stale; otherwise resume the affected
  curve independently.
- **Incompatible lookup manifest:** verify application, threshold, schedule,
  Stage 5 source, and rank map; rebuild only that method with Stage 7.2.
- **Checkpoint archived/replaced:** the required target sequence changed;
  completed exact lookups remain reusable.
- **Geometry mismatch:** regenerate or correctly export the Zonation rank map
  on the Stage 5/6 grid.
- **OpenMP header/runtime missing:** install an OpenMP toolchain and clear only
  `Data/Cache/Rcpp` before recompiling.
- **GRASS initialization failure:** set `grass_dir` or use `clump_backend:
  "terra"` deliberately.
- **Runtime-log write failure:** confirm the process can create and append the
  documented `Logs` directory and file. The active operation stops before it
  can continue without the selected diagnostic.
- **Runtime session has no end record:** the process was terminated or
  interrupted before normal/error unwinding. Inspect the last checkpoint and
  resume through the owning stage; never edit scientific artifacts to repair a
  log.
- **Large runtime logs:** archive or delete them manually between invocations if
  desired. This does not alter reuse or compatibility, and the next active
  operation creates/appends its canonical path automatically.
- **Presentation chunk error:** run the named `results` or `calculate` chunk in
  the current session.

## 22. Objects intentionally regenerated rather than cached

The repository deliberately does **not** persist Stage 7.3/7.4/8.1 PersistRank or
benchmark persistence trajectories, PU trajectories, species trajectories,
summary statistics, target-plan CSVs, lookup-index CSVs, report metadata,
formatted tables, robustness matrices, or HTML fragments. Target plans and
lookup indexes are cheap views of authoritative Stage 6 states and exact lookup
files. Regenerating these objects is what allows report settings to change
without repeating benchmark spatial reconstruction.
