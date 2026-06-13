# =============================================================================
# 03_bioclim_compare.R  (PALEO bioclim version)
# Three ways to assign temperature to each bison specimen:
#
#   Greenland (GISP2)         single global proxy        -> TIME only
#   Paleo-BioClim bio1 (TraCE) free-running climate model -> SPACE+TIME, unconstrained
#   PReSto DA tas             data assimilation          -> SPACE+TIME, obs-constrained
#
# Paleo-bioclim = TraCE-derived 1-ka time-slice bricks in
#   MSU/Thesis/Bison/MaxEntPoints/<NN_MM>/bioclim<NN_MM>.grd  (layers bio1..bio19).
# Each specimen is matched to its 1-ka slice by age, then bio1/bio5/bio12/bio15
# are extracted at its lon/lat. Precip (bio12) here is PALEO (TraCE), used as a
# provisional hydroclimate axis until PReSto exposes DA precipitation.
# NOTE: bioclim precip units are TraCE-native (relative); all predictors z-scored.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(ggplot2); library(tidyr); library(lme4)
  if (requireNamespace("lmerTest", quietly = TRUE)) library(lmerTest)
  has_mumin <- requireNamespace("MuMIn", quietly = TRUE)
})

proj   <- "."
outdir <- file.path(proj, "analysis", "outputs")
mxdir  <- "${PALEOBIOCLIM:-./paleo-bioclim}"

sp <- read.csv(file.path(outdir, "bison_holocene_with_presto.csv"))

# ---- match each specimen to its 1-ka paleo slice, extract paleo-bioclim -------
sp$slice  <- pmax(1L, ceiling(sp$age_bp / 1000))            # 0-999->1 (01_00), 9000-9999->10 (10_09)
sp$folder <- sprintf("%02d_%02d", sp$slice, sp$slice - 1L)
want <- c("bio1","bio5","bio12","bio15")
sp[want] <- NA_real_
for (fo in sort(unique(sp$folder))) {
  grd <- file.path(mxdir, fo, sprintf("bioclim%s.grd", fo))
  if (!file.exists(grd)) { message("missing slice: ", fo); next }
  b   <- rast(grd)[[want]]
  idx <- which(sp$folder == fo)
  vals <- terra::extract(b, vect(sp[idx, c("lon","lat")], geom = c("lon","lat"), crs = "EPSG:4326"))[, -1]
  sp[idx, want] <- vals
}
sp <- sp |> rename(bioclim_mat = bio1, bioclim_tmax = bio5,
                   bioclim_precip = bio12, bioclim_pseas = bio15)
sp$locality <- factor(sp$locality)
n_na <- sum(is.na(sp$bioclim_mat))
if (n_na) cat(sprintf("Note: %d specimens missing paleo-bioclim (off-grid) — dropped.\n", n_na))
sp <- sp |> filter(!is.na(bioclim_mat))
cat(sprintf("Paleo-bioclim extracted for %d specimens (bio1 %.1f..%.1f C)\n",
            nrow(sp), min(sp$bioclim_mat), max(sp$bioclim_mat)))
write.csv(sp, file.path(outdir, "bison_holocene_with_bioclim.csv"), row.names = FALSE)

# ---- 3-panel temperature comparison -----------------------------------------
theme_set(theme_bw(base_size = 12))
long <- sp |>
  mutate(age_ka = age_bp/1000) |>
  select(mass_kg, age_ka,
         `Greenland GISP2`           = gisp_temp,
         `Paleo-bioclim MAT (TraCE)` = bioclim_mat,
         `PReSto DA anomaly`         = presto_tas) |>
  pivot_longer(-c(mass_kg, age_ka), names_to = "source", values_to = "temp_C") |>
  mutate(source = factor(source, levels = c(
    "Greenland GISP2","Paleo-bioclim MAT (TraCE)","PReSto DA anomaly"))) |>
  group_by(source) |> mutate(temp_z = as.numeric(scale(temp_C))) |> ungroup()  # standardize for comparability
p3 <- ggplot(long, aes(temp_z, mass_kg)) +
  geom_point(aes(colour = age_ka), alpha = .55, size = 1) +
  geom_smooth(method = "lm", colour = "black", linewidth = .7) +
  scale_colour_viridis_c("Age (ka)", direction = -1) +
  facet_wrap(~ source) +
  labs(x = "Temperature (standardized, SD units)", y = "Body mass (kg)",
       title = "Bison body mass vs. three temperature estimates",
       subtitle = "n = 694 (0–12 ka). Points coloured by specimen age; black line = OLS fit.")
ggsave(file.path(outdir, "fig_3panel_temp.png"), p3, width = 11, height = 4.4, dpi = 150)

# ---- models ------------------------------------------------------------------
# Fit with ML (REML = FALSE): REML log-likelihoods are NOT comparable across
# models with different fixed effects, so AIC across these models requires ML.
z  <- function(x) as.numeric(scale(x))
ml <- function(f) lmer(f, data = sp, REML = FALSE)
m_age     <- ml(mass_kg ~ z(age_bp)                         + (1|locality))
m_gisp    <- ml(mass_kg ~ z(gisp_temp)                      + (1|locality))
m_bioT    <- ml(mass_kg ~ z(bioclim_mat)                    + (1|locality))
m_presto  <- ml(mass_kg ~ z(presto_tas)                     + (1|locality))
m_inter   <- ml(mass_kg ~ z(bioclim_tmax)*z(bioclim_precip) + (1|locality))
m_pp      <- ml(mass_kg ~ z(presto_tas) + z(bioclim_precip) + (1|locality))
m_pxp     <- ml(mass_kg ~ z(presto_tas) * z(bioclim_precip) + (1|locality))
m_agepp   <- ml(mass_kg ~ z(age_bp)     + z(bioclim_precip) + (1|locality))

mods <- list(`age (time)`=m_age, `GISP2 Greenland`=m_gisp,
             `Paleo-BioClim temp (TraCE)`=m_bioT, `PReSto temp (DA)`=m_presto,
             `paleo Tmax x precip`=m_inter, `PReSto temp + paleo precip`=m_pp,
             `PReSto temp x paleo precip`=m_pxp, `age + paleo precip`=m_agepp)
tab <- data.frame(model = names(mods),
                  AIC = sapply(mods, AIC),
                  R2m = if (has_mumin) sapply(mods, \(m) MuMIn::r.squaredGLMM(m)[1]) else NA) |>
  arrange(AIC); tab$dAIC <- tab$AIC - min(tab$AIC)

sink(file.path(outdir, "model_comparison.txt"))
cat("Bison body mass — model comparison (n =", nrow(sp), "; 0-12 ka). Paleo-bioclim = TraCE.\n")
cat("Precip is PALEO TraCE (provisional until PReSto DA precip). All predictors z-scored.\n")
cat("===========================================================================\n")
print(format(tab, digits = 4), row.names = FALSE)
cat("\n--- single-predictor temperature coefficients (per SD) ---\n")
for (nm in c("GISP2 Greenland","Paleo-BioClim temp (TraCE)","PReSto temp (DA)")) {
  cf <- summary(mods[[nm]])$coefficients; cat(sprintf("  %-26s slope=%7.1f  t=%5.2f\n", nm, cf[2,1], cf[2,"t value"]))
}
cat("\n--- does paleo precip improve fit? (likelihood-ratio tests) ---\n")
print(anova(m_age,    m_agepp))   # precip on top of time
print(anova(m_presto, m_pp))      # precip on top of PReSto temp
cat("\n--- paleo Tmax x precip (z-scored) ---\n"); print(round(summary(m_inter)$coefficients,2))
sink()
cat("\nModel comparison:\n"); print(format(tab, digits = 4), row.names = FALSE)

# ---- model-comparison figure -------------------------------------------------
tab2 <- tab; tab2$model <- factor(tab2$model, levels = rev(tab2$model))
pa <- ggplot(tab2, aes(dAIC, model, fill = dAIC)) +
  geom_col(width=.72) + geom_text(aes(label = sprintf("R²m=%.2f", R2m)), hjust=-0.12, size=3.4) +
  scale_fill_viridis_c(guide="none") + xlim(0, max(tab2$dAIC)*1.28) +
  labs(x="ΔAIC (0 = best)", y=NULL,
       title="Predictors of Holocene bison body mass — model comparison",
       subtitle="ML mixed models (random intercept by locality); n = 694; predictors z-scored; labels = marginal R².")
ggsave(file.path(outdir, "fig_model_comparison.png"), pa, width = 9.6, height = 4.4, dpi = 150)
cat("\nWrote fig_3panel_temp.png, fig_model_comparison.png, model_comparison.txt\n")
