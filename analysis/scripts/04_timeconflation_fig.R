# =============================================================================
# 04_timeconflation_fig.R
# Body mass vs age and two temperature estimates, as three independent panels
# (patchwork) so each gets its own x-axis title and the age panel can run
# old -> young (0 ka on the RIGHT; time advances left to right).
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
outdir <- "./analysis/outputs"
d <- read.csv(file.path(outdir, "bison_holocene_with_presto.csv")) |> mutate(age_ka = age_bp/1000)
theme_set(theme_bw(base_size = 12))

base <- function(df, xvar, xlab, reverse = FALSE) {
  p <- ggplot(df, aes(.data[[xvar]], mass_kg)) +
    geom_point(aes(colour = age_ka), alpha = .55, size = 1) +
    geom_smooth(method = "lm", colour = "black", linewidth = .7) +
    scale_colour_viridis_c("Age (ka)", direction = -1) +
    labs(x = xlab, y = "Body mass (kg)")
  if (reverse) p <- p + scale_x_reverse()   # 0 on the right; time L->R
  p
}
pA <- base(d, "age_ka",     "Age (ka cal BP)",                 reverse = TRUE)
pB <- base(d, "presto_tas", "PReSto DA temp anomaly (°C)") + labs(y = NULL)
pC <- base(d, "gisp_temp",  "Greenland GISP2 temp (°C)")   + labs(y = NULL)

fig <- (pA | pB | pC) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Bison body mass vs. age and two temperature estimates",
    subtitle = "n = 694 (0–12 ka). Points coloured by specimen age; black line = OLS fit.") &
  theme(legend.position = "right")

ggsave(file.path(outdir, "fig_timeconflation.png"), fig, width = 11, height = 4.2, dpi = 150)
cat("wrote fig_timeconflation.png (age axis reversed: 0 ka on right; per-panel x titles)\n")
