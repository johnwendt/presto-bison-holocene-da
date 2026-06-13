# =============================================================================
# 01_prep_bodysize.R
# PReSto workshop — bison body size vs. climate
#
# Purpose: read Martin et al. 2018 (Ecology & Evolution) supplement, build a
# clean per-specimen body-size dataset, reproduce Martin's body-mass ~ Greenland
# (GISP2) temperature result, and write the Holocene (0-12 ka) subset that will
# be joined to the PReSto Holocene-DA temperature field in 02_extract_refit.R.
#
# Runs WITHOUT the reconstruction — do this while the DA job is in the queue.
#
# Data note: the column labelled 'Cal Age (ka)' is mislabelled — its values are
# in YEARS cal BP (identical to 'IntCal13 (CAL BP)'). We use CAL BP throughout.
# Body mass follows Martin: M(kg) = (DstL/11.49)^3  (matches their sheet-3 stats).
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(ggplot2); library(lme4)
  if (requireNamespace("lmerTest", quietly = TRUE)) library(lmerTest)  # p-values if available
})

# ---- paths -------------------------------------------------------------------
proj   <- "."
xlsx   <- file.path(proj, "ece34019-sup-0001-datas1-s3.xlsx")
outdir <- file.path(proj, "analysis", "outputs")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- read & clean ------------------------------------------------------------
raw <- read_excel(xlsx, sheet = "1_Database")

dat <- raw |>
  transmute(
    locality  = trimws(.data[["Collection/Population Locality"]]),
    specimen  = as.character(SpecimenNumber),
    age_bp    = as.numeric(.data[["IntCal13 (CAL BP)"]]),   # years cal BP (authoritative)
    age_lo    = as.numeric(.data[["Age bound (lower)"]]),
    age_hi    = as.numeric(.data[["Age bound (upper)"]]),
    lat       = as.numeric(Lat),
    lon       = as.numeric(Long),                            # already negative (W)
    gisp_temp = as.numeric(.data[["GISP2 Temp (IntCal13)"]]),# Martin's predictor
    genus     = as.character(Genus),
    species   = as.character(Species),
    dstl      = as.numeric(DstL)
  ) |>
  filter(!is.na(dstl), !is.na(lat), !is.na(lon)) |>
  mutate(
    mass_kg = (dstl / 11.49)^3,                              # Martin allometry
    species = factor(species),
    holo    = age_bp >= 0 & age_bp <= 12000                  # PReSto-DA window
  )

cat(sprintf("Specimens with DstL+coords: %d  | localities: %d\n",
            nrow(dat), dplyr::n_distinct(dat$locality)))
cat(sprintf("Holocene (0-12 ka): %d  | Pleistocene (>12 ka): %d\n",
            sum(dat$holo), sum(!dat$holo)))

# ---- sanity check: mass formula vs. published locality means -----------------
# (eyeball a few against sheet '3_Locality Stats' "Mass Derived from DstL" Mean)
chk <- dat |> group_by(locality) |>
  summarise(n = n(), dstl_mean = mean(dstl), mass_mean = mean(mass_kg), .groups = "drop") |>
  arrange(desc(n))
cat("\nMass-formula check (compare mass_mean to sheet 3):\n"); print(head(chk, 6))

# ---- reproduce Martin: body mass ~ GISP2 (Greenland) temperature, FULL data --
m_martin <- lm(mass_kg ~ gisp_temp, data = dat)
cat("\n=== Martin-style FULL-DATA regression: mass ~ GISP2 temp ===\n")
print(summary(m_martin)$coefficients)
cat(sprintf("r = %.3f, n = %d\n",
            cor(dat$gisp_temp, dat$mass_kg, use = "complete.obs"), nrow(dat)))

# ---- Holocene-only baseline with the SAME (GISP2) predictor ------------------
# This is the "old way" we will contrast against PReSto local temperature.
holo <- dat |> filter(holo)
m_gisp_holo <- lmer(mass_kg ~ gisp_temp + (1 | locality), data = holo)
cat("\n=== Holocene-only (n=", nrow(holo), "): mass ~ GISP2 temp + (1|locality) ===\n", sep = "")
print(summary(m_gisp_holo)$coefficients)

# NOTE: species is deliberately NOT used as a predictor. The antiquus/occidentalis/
# bison forms are an anagenetic continuum and body size is itself a diagnostic
# character for them — conditioning on species would partial out the variation we
# are trying to explain (circular). Species is retained only as a descriptive
# colour in the EDA plots below.

# ---- EDA figures -------------------------------------------------------------
theme_set(theme_bw(base_size = 12))

p_age <- ggplot(dat, aes(age_bp/1000, mass_kg, colour = species)) +
  geom_point(alpha = .5, size = 1) +
  geom_vline(xintercept = 12, linetype = 2) +
  scale_x_reverse() +
  labs(x = "Age (ka cal BP)", y = "Body mass (kg)",
       title = "Bison body mass through time (dashed = 12 ka PReSto-DA limit)")

p_gisp <- ggplot(dat, aes(gisp_temp, mass_kg, colour = species)) +
  geom_point(alpha = .5, size = 1) + geom_smooth(method = "lm", colour = "black") +
  labs(x = "GISP2 (Greenland) temperature (°C)", y = "Body mass (kg)",
       title = "Martin's predictor: a single hemispheric proxy")

ggsave(file.path(outdir, "eda_mass_vs_age.png"),  p_age,  width = 7, height = 4.5, dpi = 150)
ggsave(file.path(outdir, "eda_mass_vs_gisp.png"), p_gisp, width = 6, height = 4.5, dpi = 150)

# ---- write Holocene subset for the extraction step ---------------------------
holo_out <- holo |>
  transmute(locality, specimen, age_bp, age_lo, age_hi,
            lat, lon, mass_kg, dstl, gisp_temp, genus, species)
write.csv(holo_out, file.path(outdir, "bison_holocene_specimens.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s  (n=%d specimens, 0-12 ka)\n",
            file.path(outdir, "bison_holocene_specimens.csv"), nrow(holo_out)))

saveRDS(list(all = dat, holo = holo, m_martin = m_martin), file.path(outdir, "prep.rds"))
cat("Done. Next: run 02_extract_refit.R once the PReSto NetCDF arrives.\n")
