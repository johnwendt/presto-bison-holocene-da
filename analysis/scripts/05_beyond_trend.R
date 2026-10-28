# =============================================================================
# 05_beyond_trend.R
# Tests behind the Presto paper's workshop paragraph (McKay et al., Sect. 3.5.3):
#   1. does either temperature series add to a model that already has age?
#   2. latitude added to age (the DA field is an anomaly field, so latitude is
#      the only spatial test of a Bergmann cline here)
#   3. spread of DA R2m and of the age-adjusted test across ensemble members
#   4. lagged thermal response: exponentially weighted memory of older local
#      temperature (memory truncated at 12 ka, where the reconstruction starts)
#
# Run from the analysis/ directory:   Rscript scripts/05_beyond_trend.R
# Needs the reconstruction NetCDF (~320 MB). Set PRESTO_NC to a local copy, or
# the script downloads it from this repository's reconstruction release.
# Needs: dplyr, lme4, MuMIn, ncdf4
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(lme4); library(MuMIn); library(ncdf4) })

nc_name <- "holocene_recon_2026-06-04_17-31-28.730862_bison-body-size-temp.nc"
nc_url  <- paste0("https://github.com/johnwendt/presto-bison-holocene-da/releases/download/",
                  "recon-26968202830/", nc_name)
nc_path <- Sys.getenv("PRESTO_NC", file.path("data", nc_name))
if (!file.exists(nc_path)) {
  options(timeout = 3600)
  message("Downloading reconstruction NetCDF to ", nc_path)
  download.file(nc_url, nc_path, mode = "wb")
}

sp <- read.csv("data/bison_holocene_with_bioclim.csv")
sp$locality <- factor(sp$locality)
z   <- function(x) as.numeric(scale(x))
ml  <- function(f, d = sp) lmer(f, data = d, REML = FALSE)
lrt <- function(m0, m1) anova(m0, m1)[2, "Pr(>Chisq)"]

dir.create("results", showWarnings = FALSE)
sink("results/beyond_trend.txt")
cat(sprintf("n = %d specimens, %d localities (one age per locality). ML fits.\n\n",
            nrow(sp), nlevels(sp$locality)))

m_age <- ml(mass_kg ~ z(age_bp) + (1 | locality))
cat(sprintf("age only: R2m = %.2f\n", r.squaredGLMM(m_age)[1]))

cat("\n1. temperature added to age (likelihood-ratio p):\n")
for (v in c("presto_tas", "gisp_temp")) {
  m <- ml(as.formula(sprintf("mass_kg ~ z(age_bp) + z(%s) + (1 | locality)", v)))
  cat(sprintf("   %-11s coef/SD = %6.1f  p = %.2f\n", v, fixef(m)[3], lrt(m_age, m)))
}

m_lat <- ml(mass_kg ~ z(age_bp) + z(lat) + (1 | locality))
cat(sprintf("\n2. latitude added to age: p = %.2f; cor(lat, mass) = %.2f\n",
            lrt(m_age, m_lat), cor(sp$lat, sp$mass_kg)))

nc <- nc_open(nc_path)
ens <- ncvar_get(nc, "recon_tas_ens"); mu <- ncvar_get(nc, "recon_tas_mean")
ages <- nc$dim$ages$vals; lats <- nc$dim$lat$vals; lons <- nc$dim$lon$vals
nc_close(nc)
ia <- sapply(sp$age_bp,     \(x) which.min(abs(ages - x)))
it <- sapply(sp$lat,        \(x) which.min(abs(lats - x)))
io <- sapply(sp$lon %% 360, \(x) which.min(abs(lons - x)))

r2 <- p <- numeric(dim(ens)[3])
for (e in seq_along(r2)) {
  sp$t_e <- ens[cbind(io, it, e, ia)]
  r2[e] <- r.squaredGLMM(ml(mass_kg ~ z(t_e) + (1 | locality)))[1]
  p[e]  <- lrt(m_age, ml(mass_kg ~ z(age_bp) + z(t_e) + (1 | locality)))
}
cat(sprintf("\n3. %d ensemble members: DA R2m median %.2f (95%% %.2f-%.2f); temp-given-age p min %.2f, median %.2f\n",
            length(r2), median(r2), quantile(r2, .025), quantile(r2, .975), min(p), median(p)))

ew <- function(i, tau) { k <- ia[i]:length(ages); w <- exp(-(ages[k] - ages[ia[i]]) / (tau * 1000))
  sum(w * mu[cbind(io[i], it[i], k)]) / sum(w) }
cat("\n4. lagged temperature (memory tau, kyr): dAIC vs age-only, R2m, p(lag | age), p(age | lag)\n")
for (tau in c(0.5, 1, 2, 4, 8)) {
  sp$t_lag <- sapply(seq_len(nrow(sp)), ew, tau = tau)
  m_l  <- ml(mass_kg ~ z(t_lag) + (1 | locality))
  m_al <- ml(mass_kg ~ z(age_bp) + z(t_lag) + (1 | locality))
  cat(sprintf("   tau=%3.1f  dAIC=%+5.1f  R2m=%.2f  p(lag|age)=%.2f  p(age|lag)=%.3f\n", tau,
              AIC(m_l) - AIC(m_age), r.squaredGLMM(m_l)[1], lrt(m_age, m_al), lrt(m_l, m_al)))
}
sink()
cat(readLines("results/beyond_trend.txt"), sep = "\n")
