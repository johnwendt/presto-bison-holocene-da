# =============================================================================
# make_recon_map.R — North American map of the reconstructed 6 ka temperature
# anomaly, with the body-size localities overlaid.
#
# Unlike reproduce.R, this needs the reconstruction NetCDF (it is too large to
# commit; it is archived as an asset on the `recon-26968202830` GitHub release).
# Point PRESTO_NC at it, or drop the .nc anywhere under the repo, then:
#   Rscript analysis/make_recon_map.R
# Needs: ncdf4, ggplot2, dplyr, maps.
# =============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(ggplot2); library(dplyr); library(maps) })

ncp <- Sys.getenv("PRESTO_NC", "")
if (ncp == "") {
  cand <- list.files(".", pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)
  cand <- cand[grepl("recon|holocene", cand, ignore.case = TRUE)]
  if (length(cand)) ncp <- cand[1]
}
stopifnot("Set PRESTO_NC to the reconstruction NetCDF (recon-26968202830 release asset)" = file.exists(ncp))

nc   <- nc_open(ncp)
lon  <- ncvar_get(nc, "lon"); lat <- ncvar_get(nc, "lat"); ages <- ncvar_get(nc, "ages")
tas  <- ncvar_get(nc, "recon_tas_mean")          # [lon, lat, ages]
nc_close(nc)

k  <- which.min(abs(ages - 6000))                # nearest 100-yr bin to 6 ka
df <- expand.grid(lon = lon, lat = lat)
df$tas <- as.vector(tas[, , k])
df$lon <- ifelse(df$lon > 180, df$lon - 360, df$lon)   # 0–360 -> -180..180
df <- df |> filter(lon >= -170, lon <= -50, lat >= 25, lat <= 75)

ll <- read.csv("analysis/data/bison_holocene_with_bioclim.csv") |>
  group_by(locality) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")

zr <- max(abs(df$tas), na.rm = TRUE)
p <- ggplot(df, aes(lon, lat, fill = tas)) +
  geom_raster() +
  scale_fill_gradient2("ΔT (°C)", low = "#2166AC", mid = "#F7F7F7",
                       high = "#B2182B", midpoint = 0, limits = c(-zr, zr)) +
  borders("world", colour = "grey30", linewidth = 0.3) +
  geom_point(data = ll, aes(lon, lat), inherit.aes = FALSE,
             shape = 21, fill = "black", colour = "white", size = 1.6, stroke = 0.3) +
  coord_quickmap(xlim = c(-170, -50), ylim = c(25, 75), expand = FALSE) +
  labs(x = "Longitude (°E)", y = "Latitude (°N)",
       title = "Reconstructed temperature anomaly at 6 ka, with body-size localities") +
  theme_bw(base_size = 12)

dir.create("analysis/figures", showWarnings = FALSE)
ggsave("analysis/figures/recon_map_6ka.png", p, width = 8, height = 5.4, dpi = 150)
cat("Wrote analysis/figures/recon_map_6ka.png (6 ka =", round(ages[k]), "yr BP bin)\n")
