# =============================================================================
# 02_extract_refit.R
# PReSto workshop — bison body size vs. climate (step 2)
#
# Purpose: sample the PReSto Holocene-DA temperature field at each bison
# specimen's (lat, lon, age), then refit body-mass models to ask whether
# spatially-explicit, observation-constrained LOCAL temperature predicts body
# mass — and contrast it with Martin's single-proxy GREENLAND (GISP2) predictor.
#
# Run this AFTER the PReSto reconstruction finishes and you've downloaded the
# NetCDF from the generated repo (recons/ ... .nc).
#
# The NetCDF variable/coord names are auto-detected; if detection fails, set the
# names by hand in the CONFIG block below.
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4); library(dplyr); library(ggplot2); library(lme4); library(tidyr)
  if (requireNamespace("lmerTest", quietly = TRUE)) library(lmerTest)  # p-values if available
  has_mumin <- requireNamespace("MuMIn", quietly = TRUE)
})

# ---- CONFIG ------------------------------------------------------------------
proj   <- "."
outdir <- file.path(proj, "analysis", "outputs")
nc_path <- Sys.getenv("PRESTO_NC", "")               # <- set to the downloaded .nc
if (nc_path == "") {
  cands <- list.files(proj, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)
  cands <- cands[grepl("recon|holocene|tas", cands, ignore.case = TRUE)]
  if (length(cands)) nc_path <- cands[1]
}
stopifnot("Set nc_path / PRESTO_NC to the reconstruction NetCDF" = file.exists(nc_path))
cat("Using NetCDF:", nc_path, "\n")

# Optional manual overrides (leave NULL to auto-detect):
VAR_MEAN <- NULL   # e.g. "recon_tas_mean"  (gridded ensemble-mean tas)
VAR_ENS  <- NULL   # e.g. "recon_tas_ens"   (selected ensemble members; for uncertainty)

# ---- open & describe ---------------------------------------------------------
nc <- nc_open(nc_path)
cat("\n--- NetCDF variables ---\n")
for (v in names(nc$var)) cat(sprintf("  %-22s [%s]\n", v,
     paste(sapply(nc$var[[v]]$dim, function(d) d$name), collapse = ",")))
cat("--- dims ---\n")
for (d in names(nc$dim)) cat(sprintf("  %-12s n=%d\n", d, nc$dim[[d]]$len))

pick <- function(cands, pool) { hit <- pool[tolower(pool) %in% tolower(cands)]
  if (length(hit)) hit[1] else pool[grep(paste(cands, collapse="|"), pool, ignore.case=TRUE)][1] }
if (is.null(VAR_MEAN)) VAR_MEAN <- pick(c("recon_tas_mean","tas_mean","tas"), names(nc$var))
if (is.null(VAR_ENS )) VAR_ENS  <- pick(c("recon_tas_ens","tas_ens"),         names(nc$var))
cat(sprintf("\nUsing mean var: %s | ensemble var: %s\n", VAR_MEAN,
            ifelse(is.na(VAR_ENS),"<none>",VAR_ENS)))

# coordinate vectors
dimname <- function(cands) pick(cands, names(nc$dim))
d_age <- dimname(c("age","ages","time")); d_lat <- dimname(c("lat","latitude"))
d_lon <- dimname(c("lon","longitude"))
ages <- as.numeric(nc$dim[[d_age]]$vals)   # yr BP
lats <- as.numeric(nc$dim[[d_lat]]$vals)
lons <- as.numeric(nc$dim[[d_lon]]$vals)
cat(sprintf("ages: %.0f..%.0f (n=%d) | lat %.1f..%.1f | lon %.1f..%.1f\n",
            min(ages),max(ages),length(ages),min(lats),max(lats),min(lons),max(lons)))

# read mean field and align dim order -> [age, lat, lon]
arr  <- ncvar_get(nc, VAR_MEAN)
dord <- sapply(nc$var[[VAR_MEAN]]$dim, function(d) d$name)
perm <- match(c(d_age, d_lat, d_lon), dord)
field <- aperm(arr, perm)                  # [age, lat, lon]

# ensemble spread (sd across members) if available -> [age, lat, lon]
sdfield <- NULL
if (!is.na(VAR_ENS)) {
  ens  <- ncvar_get(nc, VAR_ENS)
  edn  <- sapply(nc$var[[VAR_ENS]]$dim, function(d) d$name)
  d_en <- setdiff(edn, c(d_age,d_lat,d_lon))[1]
  ensp <- aperm(ens, match(c(d_age,d_lat,d_lon,d_en), edn))
  sdfield <- apply(ensp, c(1,2,3), sd, na.rm = TRUE)
}
nc_close(nc)

# longitude convention: convert specimen lon to match the grid
lon_to_grid <- function(x) if (min(lons) >= 0) ((x %% 360)) else x
nearest <- function(x, grid) which.min(abs(grid - x))

# ---- specimens ---------------------------------------------------------------
sp <- read.csv(file.path(outdir, "bison_holocene_specimens.csv"))
sp$presto_tas <- NA_real_; sp$presto_sd <- NA_real_
sp$cell_lat <- NA_real_; sp$cell_lon <- NA_real_; sp$cell_age <- NA_real_
for (i in seq_len(nrow(sp))) {
  ia <- nearest(sp$age_bp[i], ages)
  it <- nearest(sp$lat[i],    lats)
  io <- nearest(lon_to_grid(sp$lon[i]), lons)
  sp$presto_tas[i] <- field[ia, it, io]
  if (!is.null(sdfield)) sp$presto_sd[i] <- sdfield[ia, it, io]
  sp$cell_age[i] <- ages[ia]; sp$cell_lat[i] <- lats[it]; sp$cell_lon[i] <- lons[io]
}
n_drop <- sum(is.na(sp$presto_tas))
if (n_drop) cat(sprintf("Note: %d specimens fell on NaN cells (ocean/edge) — dropped.\n", n_drop))
sp <- sp |> filter(!is.na(presto_tas))
sp$locality <- factor(sp$locality)
cat(sprintf("Extracted PReSto temperature for %d specimens.\n", nrow(sp)))
write.csv(sp, file.path(outdir, "bison_holocene_with_presto.csv"), row.names = FALSE)

# ---- models: LOCAL (PReSto) vs GREENLAND (GISP2), same data ------------------
sink(file.path(outdir, "bodysize_climate_models.txt"))
cat("PReSto bison body-size models  (n =", nrow(sp), "specimens, 0-12 ka)\n")
cat("=====================================================================\n")

# REML = FALSE (ML): the two models have different fixed-effect predictors, so their
# likelihoods (and hence AIC) are only comparable when fit by ML, not REML.
m_local <- lmer(mass_kg ~ scale(presto_tas) + (1 | locality), data = sp, REML = FALSE)
m_gisp  <- lmer(mass_kg ~ scale(gisp_temp)  + (1 | locality), data = sp, REML = FALSE)

cat("\n--- (A) LOCAL observation-constrained temperature (PReSto Holocene DA) ---\n")
print(summary(m_local)$coefficients)
cat("\n--- (B) GREENLAND single-proxy temperature (Martin's GISP2 predictor) ---\n")
print(summary(m_gisp)$coefficients)

# NOTE: species is intentionally excluded as a predictor (anagenetic continuum;
# body size is itself a species-diagnostic character -> circular). The contrast of
# interest is purely environmental: does spatially-explicit LOCAL temperature
# explain Holocene body size better than the single GREENLAND proxy?
cat("\n--- model comparison (AIC; lower = better) ---\n")
print(AIC(m_gisp, m_local))
if (has_mumin) {
  cat("\n--- R^2 (marginal / conditional) ---\n")
  cat("LOCAL :", round(MuMIn::r.squaredGLMM(m_local), 3), "\n")
  cat("GISP2 :", round(MuMIn::r.squaredGLMM(m_gisp ), 3), "\n")
}
sink()
cat("Wrote model summary -> outputs/bodysize_climate_models.txt\n")

# ---- the money figure: mass vs LOCAL temp  vs  mass vs GREENLAND temp ---------
theme_set(theme_bw(base_size = 12))
long <- sp |>
  select(mass_kg, species, `Local (PReSto DA)` = presto_tas, `Greenland (GISP2)` = gisp_temp) |>
  pivot_longer(c(`Local (PReSto DA)`, `Greenland (GISP2)`),
               names_to = "predictor", values_to = "temp_C")
p_cmp <- ggplot(long, aes(temp_C, mass_kg)) +
  geom_point(aes(colour = species), alpha = .5, size = 1) +
  geom_smooth(method = "lm", colour = "black") +
  facet_wrap(~ predictor, scales = "free_x") +
  labs(x = "Temperature (°C)", y = "Body mass (kg)",
       title = "Bison body mass vs. local DA temperature vs. single Greenland proxy",
       subtitle = "Holocene (0-12 ka). If LOCAL is flat while GREENLAND slopes, the classic signal is time-conflation.")
ggsave(file.path(outdir, "fig_local_vs_greenland.png"), p_cmp, width = 9, height = 4.5, dpi = 150)

# extraction-context figure: where each specimen landed in age
p_chk <- ggplot(sp, aes(age_bp, presto_tas)) +
  geom_errorbar(aes(ymin = presto_tas - presto_sd, ymax = presto_tas + presto_sd),
                colour = "grey70", na.rm = TRUE) +
  geom_point(aes(colour = lat), size = 1) + scale_x_reverse() +
  labs(x = "Age (yr cal BP)", y = "PReSto local T anomaly (°C)", colour = "Lat",
       title = "Extracted local temperature per specimen (±ensemble SD)")
ggsave(file.path(outdir, "fig_extracted_timeseries.png"), p_chk, width = 7.5, height = 4.5, dpi = 150)

cat("\nDone. Figures + tables in analysis/outputs/.\n")
