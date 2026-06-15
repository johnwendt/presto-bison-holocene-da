# =============================================================================
# reproduce.R  — portable, self-contained entry point
#
# Reproduces the body-size analysis (models + figures) from the committed,
# pre-joined dataset `data/bison_holocene_with_bioclim.csv`. No large rasters
# or the reconstruction NetCDF are required — those are only needed to rebuild
# the joined data from scratch (see scripts/01–04 and README.md for provenance).
#
# Run from the analysis/ directory:   Rscript reproduce.R
# Needs: dplyr, ggplot2, tidyr, lme4, patchwork (+ optional MuMIn, lmerTest)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr); library(lme4); library(patchwork)
  if (requireNamespace("lmerTest", quietly = TRUE)) library(lmerTest)
  has_mumin <- requireNamespace("MuMIn", quietly = TRUE)
})
dir.create("figures", showWarnings = FALSE); dir.create("results", showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))

d <- read.csv("data/bison_holocene_with_bioclim.csv") |> mutate(age_ka = age_bp/1000)
d$locality <- factor(d$locality)
z <- function(x) as.numeric(scale(x))
cat(sprintf("n = %d specimens (0–12 ka); %d localities\n", nrow(d), nlevels(d$locality)))

# ---- where vs when -----------------------------------------------------------
cat(sprintf("cor(mass, age)=%+.2f  cor(mass, latitude)=%+.2f  cor(mass, local MAT)=%+.2f\n",
            cor(d$mass_kg,d$age_bp), cor(d$mass_kg,d$lat), cor(d$mass_kg,d$bioclim_mat)))

# ---- models ------------------------------------------------------------------
mods <- list(
  `age (time)`                 = lmer(mass_kg ~ z(age_bp)                          + (1|locality), d, REML = FALSE),
  `GISP2 Greenland`            = lmer(mass_kg ~ z(gisp_temp)                       + (1|locality), d, REML = FALSE),
  `Paleo-BioClim temp (TraCE)` = lmer(mass_kg ~ z(bioclim_mat)                     + (1|locality), d, REML = FALSE),
  `PReSto temp (DA)`           = lmer(mass_kg ~ z(presto_tas)                      + (1|locality), d, REML = FALSE),
  `paleo Tmax x precip`        = lmer(mass_kg ~ z(bioclim_tmax)*z(bioclim_precip)  + (1|locality), d, REML = FALSE),
  `PReSto temp + paleo precip` = lmer(mass_kg ~ z(presto_tas) + z(bioclim_precip)  + (1|locality), d, REML = FALSE),
  `PReSto temp x paleo precip` = lmer(mass_kg ~ z(presto_tas) * z(bioclim_precip)  + (1|locality), d, REML = FALSE),
  `age + paleo precip`         = lmer(mass_kg ~ z(age_bp)     + z(bioclim_precip)  + (1|locality), d, REML = FALSE),
  `age + PReSto temp`          = lmer(mass_kg ~ z(age_bp)     + z(presto_tas)      + (1|locality), d, REML = FALSE),
  `age + GISP Greenland`       = lmer(mass_kg ~ z(age_bp)     + z(gisp_temp)       + (1|locality), d, REML = FALSE))
# NB: models fit with ML (REML = FALSE) so AIC is valid across different fixed effects.
tab <- data.frame(model = names(mods), AIC = sapply(mods, AIC),
                  R2m = if (has_mumin) sapply(mods, \(m) MuMIn::r.squaredGLMM(m)[1]) else NA) |>
  arrange(AIC); tab$dAIC <- tab$AIC - min(tab$AIC)
capture.output(print(format(tab, digits = 4), row.names = FALSE),
               file = "results/model_comparison.txt")
print(format(tab, digits = 4), row.names = FALSE)

# ---- does local temperature add anything once time is in the model? ----------
# The central claim ("temperature adds ~0 beyond the time trend") is tested here,
# not just inferred from the single-predictor table above.
m_age  <- mods[["age (time)"]]
m_ageP <- mods[["age + PReSto temp"]]
sink("results/time_vs_temperature.txt")
cat(sprintf("Does local temperature add anything once the time trend is included?\n"))
cat(sprintf("(ML / REML = FALSE; n = %d specimens, %d localities.)\n\n", nrow(d), nlevels(d$locality)))
cat("-- age + PReSto local temperature --\n"); print(round(summary(m_ageP)$coefficients, 3))
cat(sprintf("\nLRT  age  vs  age + PReSto temp:  dAIC = %+.2f,  p = %.3f\n",
            AIC(m_ageP) - AIC(m_age), anova(m_age, m_ageP)[2, "Pr(>Chisq)"]))
cat(sprintf("PReSto temp alone: t = %.2f (negative); controlling for time: t = %.2f (sign reverses, n.s.).\n",
            summary(mods[["PReSto temp (DA)"]])$coefficients["z(presto_tas)", "t value"],
            summary(m_ageP)$coefficients["z(presto_tas)", "t value"]))
cat("\n-- the design is essentially cross-sectional --\n")
wb <- d |> group_by(locality) |> summarise(age_sd = sd(age_bp), tas_sd = sd(presto_tas), .groups = "drop")
cat(sprintf("age varies within a locality in %d of %d localities;\n", sum(wb$age_sd > 0, na.rm = TRUE), nrow(wb)))
cat(sprintf("PReSto temperature varies within a locality in %d of %d localities.\n", sum(wb$tas_sd > 0, na.rm = TRUE), nrow(wb)))
cat(sprintf("specimen-level cor(age, PReSto temp) = %+.2f  -> age and local temperature are strongly collinear,\n", cor(d$age_bp, d$presto_tas)))
cat("so this comparison reflects differences between time-specific localities, not a within-site temporal test.\n")
sink()
cat("Wrote results/time_vs_temperature.txt (time-vs-temperature confound test)\n")

# ---- Fig 1: three temperatures (standardized) --------------------------------
long <- d |>
  select(mass_kg, age_ka, `Greenland GISP2`=gisp_temp,
         `Paleo-bioclim MAT (TraCE)`=bioclim_mat, `PReSto DA anomaly`=presto_tas) |>
  pivot_longer(-c(mass_kg, age_ka), names_to="source", values_to="t") |>
  mutate(source=factor(source, levels=c("Greenland GISP2","Paleo-bioclim MAT (TraCE)","PReSto DA anomaly"))) |>
  group_by(source) |> mutate(tz=as.numeric(scale(t))) |> ungroup()
ggsave("figures/fig_3panel_temp.png",
  ggplot(long, aes(tz, mass_kg)) +
    geom_point(aes(colour=age_ka), alpha=.55, size=1) + geom_smooth(method="lm", colour="black", linewidth=.7) +
    scale_colour_viridis_c("Age (ka)", direction=-1) + facet_wrap(~source) +
    labs(x="Temperature (standardized, SD units)", y="Body mass (kg)",
         title="Bison body mass vs. three temperature estimates"),
  width=11, height=4.4, dpi=150)

# ---- Fig 2: model comparison -------------------------------------------------
t2 <- tab; t2$model <- factor(t2$model, levels=rev(t2$model))
ggsave("figures/fig_model_comparison.png",
  ggplot(t2, aes(dAIC, model, fill=dAIC)) +
    geom_col(width=.72) + geom_text(aes(label=sprintf("R²m=%.2f", R2m)), hjust=-0.12, size=3.4) +
    scale_fill_viridis_c(guide="none") + xlim(0, max(t2$dAIC)*1.28) +
    labs(x="ΔAIC (0 = best)", y=NULL, title="Predictors of Holocene bison body mass"),
  width=9.6, height=4.4, dpi=150)

# ---- Fig 3: time-conflation (age axis reversed: 0 ka on right) ----------------
base <- function(df,x,xlab,rev=FALSE){ p <- ggplot(df, aes(.data[[x]], mass_kg)) +
  geom_point(aes(colour=age_ka), alpha=.55, size=1) + geom_smooth(method="lm", colour="black", linewidth=.7) +
  scale_colour_viridis_c("Age (ka)", direction=-1) + labs(x=xlab, y="Body mass (kg)")
  if (rev) p <- p + scale_x_reverse(); p }
fig <- (base(d,"age_ka","Age (ka cal BP)",TRUE) | (base(d,"presto_tas","PReSto DA temp anomaly (°C)")+labs(y=NULL)) |
        (base(d,"gisp_temp","Greenland GISP2 temp (°C)")+labs(y=NULL))) +
  plot_layout(guides="collect") +
  plot_annotation(title="Bison body mass vs. age and two temperature estimates") &
  theme(legend.position="right")
ggsave("figures/fig_timeconflation.png", fig, width=11, height=4.2, dpi=150)

# ---- Fig 4: reconstructed local temperature through time (with ensemble SD) ---
loc <- d |> group_by(locality) |>
  summarise(age_bp=first(age_bp), presto_tas=first(presto_tas),
            presto_sd=first(presto_sd), lat=first(lat), .groups="drop")
ggsave("figures/fig_extracted_temp.png",
  ggplot(loc, aes(age_bp, presto_tas)) +
    geom_errorbar(aes(ymin=presto_tas-presto_sd, ymax=presto_tas+presto_sd),
                  colour="grey75", width=0, na.rm=TRUE) +
    geom_point(aes(colour=lat), size=2) + scale_x_reverse() +
    scale_colour_viridis_c("Latitude (°N)") +
    labs(x="Age (yr cal BP)", y="Reconstructed local temperature anomaly (°C)",
         title="Reconstructed local temperature through time (±ensemble SD)"),
  width=8, height=4.6, dpi=150)

# ---- Fig 5: body-mass R² vs. each predictor's shared variance with age --------
r2age <- function(x) summary(lm(x ~ d$age_bp))$r.squared
conf <- data.frame(
  predictor = c("Age (time)","PReSto DA anomaly","Greenland GISP2","BioClim MAT (TraCE)"),
  shared    = c(1, r2age(d$presto_tas), r2age(d$gisp_temp), r2age(d$bioclim_mat)),
  r2m       = c(tab$R2m[tab$model=="age (time)"],          tab$R2m[tab$model=="PReSto temp (DA)"],
                tab$R2m[tab$model=="GISP2 Greenland"],      tab$R2m[tab$model=="Paleo-BioClim temp (TraCE)"]),
  hj        = c(1.12, -0.12, -0.12, -0.12))
ggsave("figures/fig_age_confound.png",
  ggplot(conf, aes(shared, r2m)) +
    geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", colour="grey70", linewidth=0.5) +
    geom_point(size=3, colour="#2C7FB8") +
    geom_text(aes(label=predictor, hjust=hj), vjust=-0.7, size=3.4) +
    expand_limits(x=c(-0.05,1.18), y=c(0,0.40)) +
    labs(x="Shared variance with age,  R²(predictor ~ age)",
         y="Body-mass variance explained  (marginal R²)",
         title="Body-mass variance explained vs. shared variance with age"),
  width=7.4, height=5, dpi=150)

cat("Done. Figures in figures/, tables in results/.\n")
