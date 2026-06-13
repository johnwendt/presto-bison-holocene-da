# Does warming shrink bison?
### Testing the thermal hypothesis among competing explanations, with a custom PReSto reconstruction
John A. F. Wendt — New Mexico State University · PReSto Workshop 2026

> *Slides as delivered at the workshop. The repository README reflects the current, broader framing — bison diminution as a question with several competing hypotheses (thermal, forage, density, hunting); see `../README_NOTES.md`.*

Figures live in `analysis/figures/`. Reproduce the models and figures with `reproduce.R` (portable; no large data needed); provenance scripts are `scripts/01`–`04`.

---

## 1 · Background & question

- North American bison body mass fell **~37% since the Last Glacial Maximum** (Martin et al. 2018, *Ecol. Evol.*).
- Proposed mechanism: **Bergmann's rule — warming → smaller**. Evidence: body mass regressed on **Greenland (GISP2) ice-core temperature**; extrapolated **+4 °C → −46% body mass**.
- The catch: a single hemispheric proxy varies **only in time**, so *temperature* and *time* are perfectly confounded — any trait declining through the Holocene will correlate with it.
- **Question:** does spatially-explicit, observation-constrained paleotemperature actually support a thermal driver of body size?

> *Speaker note:* This 46% number gets cited as a climate-warming warning for large grazers. I tested whether the temperature signal is real or an artifact of the single-proxy design.

---

## 2 · A custom PReSto Holocene reconstruction
**figure: `recon_map_6ka.png`** (mean ΔT at 6 ka)

Built in the **PReSto Custom Engine** — Holocene paleoclimate **data assimilation** (offline DA; Erb et al. 2022):

- **Target variable:** annual-mean surface temperature anomaly (`tas`)
- **Proxies assimilated:** Temp12k (Kaufman et al. 2020), North America — **423 records / 173 datasets**; all seasonalities; ≤500-yr resolution; ≥100-yr length
- **Model prior:** HadCM3 + TraCE-21ka transient simulations (multi-model); **104-member ensemble**
- **Domain / time:** 0–12,000 cal BP, **100-yr steps (120 bins)**, anomalies vs 3–5 ka, ~2.8° grid
- **Product:** gridded ensemble-mean temperature **+ per-cell uncertainty**, fully reproducible (config + repo + release-archived NetCDF)

> *Speaker note:* DA blends proxies with a climate-model prior via a Kalman update — it yields a complete, uncertainty-quantified field that varies in **space and time**, which a single ice core cannot.

---

## 3 · Data & statistical approach

- **Body size:** Martin et al. (2018) — calcaneal tuber length (DstL) in **849 specimens**; Holocene subset **n = 694** (0–12 ka). Body mass M = (DstL / 11.49)³.
- **Three temperatures per specimen:**
  1. **Greenland GISP2** (Martin's predictor) — time only
  2. **Paleo-bioclim MAT** (TraCE 1-ka slices, extracted at locality) — local, free-running model
  3. **PReSto DA anomaly** (extracted at lat/lon, nearest 100-yr bin) — local, observation-constrained
- **Models:** linear mixed models, `mass ~ z(predictor) + (1 | locality)` (lme4); compared by **AIC** and **marginal/conditional R²** (MuMIn). Predictors z-scored.
- **Species excluded as a predictor:** *antiquus/occidentalis/bison* are an anagenetic continuum and body size is itself a species-diagnostic character → conditioning on species is circular.

> *Speaker note:* The random intercept handles multiple specimens per site. Marginal R² = variance explained by the climate predictor; conditional R² adds the locality effect.

---

## 4 · Result 1 — where vs. when
**figure: `fig_3panel_temp.png`**

- **Local temperature — where a bison actually lived — is flat.** cor(mass, latitude) = **−0.04**; cor(mass, local MAT) = **+0.04** → **no spatial Bergmann**: bison from colder *places* are not bigger.
- Only the *time-tracking* temperatures slope (Greenland, PReSto), and the **age colour-gradient runs along the x-axis**.
- cor(mass, **age**) = **+0.66** → the body-size variation is **temporal**, not spatial.

> *Speaker note:* Body size doesn't care where it was warm. It only "responds" to temperatures that happen to track the calendar.

---

## 5 · Result 2 — model comparison & decomposition
**figure: `fig_model_comparison.png`**

| Predictor | ΔAIC | marginal R² |
|---|---|---|
| **age (time)** | **0** | 0.35 |
| age + paleo precip | 1 *(n.s.)* | 0.35 |
| PReSto temp (DA) | 16 | 0.12 |
| Greenland (GISP2) | 23 | 0.04 |
| Paleo-bioclim MAT (TraCE) | 28 | ~0.00 |

*(ML mixed models; AIC valid across fixed effects.)*

- **Time wins.** No temperature source beats age. Constraint ladder works — **DA (PReSto) > free model (TraCE) > single proxy (Greenland)** — but all lose to time.
- **Within/between decomposition** of PReSto temp: between-locality (spatial) effect t = −4.3 **but confounded with age (r = −0.89)**; within-locality (temporal) effect **null (t = 0.6)**.
- Provisional paleo-precipitation does **not** improve fit (likelihood-ratio test *p* = 0.32) — no detectable moisture signal here.

> *Speaker note:* PReSto's apparent temperature signal is entirely between time-specific localities, and that axis is age. Within a site, temperature change predicts nothing.

---

## 6 · Conclusions & next steps

- **Holocene bison body size tracks time more than temperature.** Among three temperature representations the observation-constrained DA field is the most informative, yet none improves on a simple time trend; local temperature adds ~0 once time is included.
- A simple warming-driven extrapolation is **not reproduced** when temperature is resolved in space and time — shifting attention to non-thermal drivers (forage productivity & quality, density; cf. Hill et al. 2008).
- Provisional paleo-precipitation showed **no** signal (n.s.); whether **productivity / effective moisture** drives the trend is an **open question** — raw annual precip is a poor proxy for forage productivity.
- **Next:** observation-constrained PReSto **hydroclimate** + a productivity-appropriate metric to test moisture directly.
- **PReSto value:** only a spatially + temporally resolved, observation-constrained reconstruction can separate temperature from time — and the whole analysis was an afternoon, end-to-end (custom DA + R).

> *Takeaway:* When you give each bison the temperature it actually experienced, "warming shrinks bison" disappears — size tracks **when**, not **how warm**.
