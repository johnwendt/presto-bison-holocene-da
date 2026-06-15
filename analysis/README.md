# Analysis — bison body size vs. Holocene temperature

This folder holds the analysis that motivated and used the reconstruction in this
repository: a test of whether the Holocene diminution of North American bison body
size was driven by temperature (see `../README_NOTES.md` for the scientific framing
and headline result).

## Quick start (portable)

```r
# from this analysis/ directory:
Rscript reproduce.R
```

`reproduce.R` regenerates the models and figures from the committed, pre-joined
dataset `data/bison_holocene_with_bioclim.csv` — **no large data required**. It
needs only `dplyr`, `ggplot2`, `tidyr`, `lme4`, `patchwork` (and optionally
`MuMIn` for R² and `lmerTest` for p-values).

## Contents

| Path | What |
|---|---|
| `reproduce.R` | Portable entry point: models + the five body-size figures from the joined CSV |
| `make_recon_map.R` | 6 ka temperature-anomaly map with body-size localities; needs the reconstruction NetCDF (`recon-26968202830` release asset, via `PRESTO_NC`) |
| `scripts/01_prep_bodysize.R` | Read Martin et al. (2018) supplement; build the Holocene specimen table |
| `scripts/02_extract_refit.R` | Sample this reconstruction's temperature at each specimen (lat/lon/age) |
| `scripts/03_bioclim_compare.R` | Add Greenland + TraCE paleo-bioclim temperatures; fit/compare models |
| `scripts/04_timeconflation_fig.R` | Body mass vs. age and two temperatures (age axis reversed) |
| `data/bison_holocene_with_bioclim.csv` | Joined per-specimen data (mass, coords, age, all three temperatures, paleo precip) |
| `data/bison_holocene_specimens.csv` | Holocene specimen table (pre-extraction) |
| `figures/` | `recon_map_6ka` (reconstruction + localities), `eda_mass_vs_age`, `fig_3panel_temp`, `fig_timeconflation`, `fig_model_comparison`, `fig_extracted_temp` (local temperature ±ensemble SD), `fig_age_confound` |
| `results/` | Model-comparison and coefficient tables |
| `SLIDES.md` | 6-slide workshop talk (methods + results) |

> The `scripts/` are provenance — they reference external raw data and have
> machine-specific paths (set via `PALEOBIOCLIM` etc.). Use `reproduce.R` to
> rebuild outputs without that data.

## Data sources

- **Body size:** Martin, J.M., Mead, J.I., Barboza, P.S. (2018) *Bison body size
  and climate change.* Ecology & Evolution 8:4564–4574. doi:10.1002/ece3.4019
  (Supplement S1–S3). Holocene subset n = 694 (0–12 ka); mass = (DstL/11.49)³.
- **Local temperature (this repo):** Holocene data-assimilation reconstruction
  (Temp12k → HadCM3 + TraCE-21ka; see repo root). Extracted from the release
  NetCDF; `scripts/02` auto-downloads it via the GitHub release.
- **Greenland temperature:** GISP2 column in Martin et al. (their predictor).
- **Paleo-bioclim (TraCE):** TraCE-21ka-derived 1-ka bioclim slices (external).

## Methods (brief)

Each specimen is assigned a temperature three ways (Greenland proxy; TraCE
paleo-bioclim; this reconstruction's DA anomaly). Body mass is modelled with
linear mixed models, `mass ~ z(temperature) + (1 | locality)`, compared by AIC
and marginal/conditional R². Species is **excluded** as a predictor (anagenetic
continuum; body size is itself species-diagnostic → circular).

**Result:** body size tracks time more than temperature. Among the three temperature
representations the observation-constrained DA field is the most informative, yet none
improves on a simple time trend, and local temperature adds essentially nothing once
time is included (LRT p = 0.53; see `results/time_vs_temperature.txt`). Precipitation
adds no detectable signal; productivity / effective moisture remains an open question.
Because the data are effectively cross-sectional — one age per locality, age and local
temperature collinear (r ≈ −0.89) — this *bounds* the thermal signal rather than
delivering a within-site temporal test. Full discussion in `../README_NOTES.md`.
