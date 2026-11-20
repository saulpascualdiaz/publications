# Fig. 3 — Data-driven partitioning of cortical connections by length.
#   a  SIFT2-weighted connection-length distributions (ADNI-4 CN per participant and pooled; BTC pooled),
#      with the three bands defined by the weighted tertiles of the ADNI-4 distribution
#   b  per-participant lower and upper weighted-tertile cuts, ADNI-4 vs. BTC
#   c  pooled distributions per unit log-length, with the kernel-density mode and shoulder
#   d  ADNI-4 pooled distribution weighted by SIFT2 strength or by connection count
source("scripts/s00_setup.R")

v        <- read_values("distance_bands/values.csv")
pooled   <- read_data("distance_bands/density_pooled.csv")
adni_sub <- read_data("distance_bands/adni_density_per_subject.csv")
adni_cut <- read_data("distance_bands/adni_cuts_per_subject.csv")
btc_cut  <- read_data("distance_bands/btc_cuts_per_subject.csv")
kde      <- read_data("distance_bands/kde_modes.csv")

cuts      <- v[c("cut_lo_mm", "cut_hi_mm")]
btc_cuts  <- v[c("btc_cut_lo_mm", "btc_cut_hi_mm")]
ADNI_LAB  <- sprintf("ADNI-4 CN, n = %d", nrow(adni_cut))
BTC_LAB   <- sprintf("BTC controls, n = %d", nrow(btc_cut))
BTC_COL   <- "grey30"; CUT_COL <- "grey20"
XMAX      <- 90
adni <- filter(pooled, cohort == "ADNI"); btc <- filter(pooled, cohort == "BTC"); btc10 <- filter(pooled, cohort == "BTC_min10mm")
linear <- function(d) filter(d, bin_center_mm <= XMAX + 1)

# ---- a: weighted length density and the three bands
ymax_a <- max(adni_sub$w_density, linear(btc)$w_density, linear(btc10)$w_density); ytop_a <- ymax_a * 1.22
bands <- tibble(band = BAND_LEVELS, xmin = c(0, cuts), xmax = c(cuts, XMAX),
                share = v[paste0("share_", BAND_LEVELS)], label = sprintf("%s\n%.2f", band, share))
pA <- ggplot(mapping = aes(bin_center_mm, w_density)) +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = band), alpha = .35, inherit.aes = FALSE) +
  scale_fill_manual(values = BAND_COLS) +
  geom_line(data = adni_sub, aes(group = subject), colour = "grey72", linewidth = .18, alpha = .8) +
  geom_line(data = linear(adni), colour = BLUE, linewidth = .6) +
  geom_line(data = linear(btc), colour = BTC_COL, linewidth = .5, linetype = "42") +
  geom_line(data = linear(btc10), colour = BTC_COL, linewidth = .4, linetype = "13") +
  geom_vline(xintercept = cuts, colour = CUT_COL, linetype = "22", linewidth = .35) +
  geom_text(data = bands, aes((xmin + xmax) / 2, ytop_a, label = label), size = KEY, colour = "grey25", vjust = 1, lineheight = .85) +
  annotate("text", x = cuts + 0.7, y = ymax_a * 1.04, label = sprintf("%s mm", fmt(cuts)), hjust = 0, size = KEY, colour = CUT_COL) +
  annotate("text", x = XMAX - 0.5, y = ymax_a * c(1.04, 0.95, 0.86), hjust = 1, size = KEY, colour = c(BLUE, BTC_COL, BTC_COL),
           label = c(ADNI_LAB, paste(BTC_LAB, "(dashed)"), "BTC, streamlines ≥ 10 mm (ADNI minimum) (dotted)")) +
  annotate("segment", x = btc_cuts, xend = btc_cuts, y = 0, yend = 0.03 * ytop_a, colour = BTC_COL, linewidth = .6) +
  annotate("text", x = btc_cuts + 0.6, y = 0.045 * ytop_a, label = fmt(btc_cuts), size = SMALL, colour = BTC_COL, hjust = 0, vjust = 0) +
  scale_x_continuous(breaks = seq(0, XMAX, 15), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(xlim = c(0, XMAX), ylim = c(0, ytop_a)) +
  labs(x = "Euclidean connection length (mm)", y = "SIFT2 weight fraction per mm") + theme_fig()

# ---- b: per-participant tertile cuts
cut_long <- bind_rows(mutate(adni_cut, cohort = ADNI_LAB), mutate(btc_cut, cohort = BTC_LAB)) %>%
  pivot_longer(c(cut_lo_mm, cut_hi_mm), names_to = "cut", values_to = "mm") %>%
  mutate(cohort = factor(cohort, c(ADNI_LAB, BTC_LAB)), cut = factor(cut, c("cut_lo_mm", "cut_hi_mm")))
CUT_FILL <- c(cut_lo_mm = BAND_COLS[["short"]], cut_hi_mm = BAND_COLS[["medium"]])
yr <- range(cut_long$mm); ylim_b <- c(yr[1] - 0.08 * diff(yr), yr[2] + 0.14 * diff(yr))
key_y <- ylim_b[2] - c(0.03, 0.11) * diff(ylim_b)
pB <- ggplot(cut_long, aes(cohort, mm)) +
  geom_hline(yintercept = cuts, linetype = "13", colour = "grey45", linewidth = .35) +
  geom_point(aes(fill = cut, group = cut), shape = 21, colour = "grey35", stroke = .15, size = .75, alpha = .85,
             position = position_jitterdodge(jitter.width = .28, dodge.width = .6, seed = 1)) +
  geom_boxplot(aes(group = interaction(cohort, cut)), alpha = 0, colour = "grey15", width = .42, linewidth = .3,
               coef = 0, outlier.shape = NA, position = position_dodge(width = .6)) +
  scale_fill_manual(values = CUT_FILL) +
  annotate("text", x = 2.62, y = cuts + 0.02 * diff(ylim_b), label = fmt(cuts), hjust = 1, vjust = 0, size = SMALL, colour = "grey40") +
  annotate("point", x = 1.62, y = key_y, shape = 22, size = 1.9, fill = CUT_FILL, colour = "grey35", stroke = .2) +
  annotate("text", x = 1.72, y = key_y, label = c("lower tertile cut", "upper tertile cut"), hjust = 0, size = KEY, colour = "grey30") +
  scale_x_discrete(labels = function(l) sub(", n = ", "\nn = ", l), expand = expansion(add = c(.55, .7))) +
  coord_cartesian(ylim = ylim_b) +
  labs(x = NULL, y = "weighted tertile cut (mm)") + theme_fig() + theme(axis.text.x = element_text(lineheight = .85))

# ---- c: density per unit log-length; a local maximum of low prominence is a shoulder, not a mode
XLIM <- c(3, 120)
loglen <- filter(pooled, cohort != "BTC_min10mm", bin_center_mm >= XLIM[1] - 0.5, bin_center_mm <= XLIM[2] + 0.5)
bw <- min(kde$bandwidth_log)
modes <- kde %>% filter(kind == "mode", bandwidth_log == bw) %>%
  mutate(label = ifelse(prominence_frac >= 0.5, sprintf("mode %s mm", fmt(length_mm)),
                        sprintf("shoulder %s mm\nprominence %s %%", fmt(length_mm), fmt(100 * prominence_frac, 2))))
ymax_c <- max(loglen$w_density_per_loglength); ytop_c <- ymax_c * 1.36
pC <- ggplot(mapping = aes(bin_center_mm, w_density_per_loglength)) +
  annotate("segment", x = cuts, xend = cuts, y = -Inf, yend = ymax_c * 1.02, colour = CUT_COL, linetype = "22", linewidth = .35) +
  geom_line(data = filter(loglen, cohort == "ADNI"), colour = BLUE, linewidth = .6) +
  geom_line(data = filter(loglen, cohort == "BTC"), colour = BTC_COL, linewidth = .5, linetype = "42") +
  annotate("text", x = cuts * 1.04, y = 0.02 * ytop_c, label = fmt(cuts), hjust = 0, vjust = 0, size = SMALL, colour = CUT_COL) +
  geom_point(data = modes, aes(length_mm, ymax_c * 1.08), shape = 6, size = 1.6, colour = BLUE, stroke = .5) +
  geom_text(data = modes, aes(length_mm, ymax_c * 1.14, label = label), size = SMALL, colour = BLUE, vjust = 0, lineheight = .85) +
  annotate("text", x = XLIM[1] * 1.04, y = ytop_c * c(0.99, 0.92), hjust = 0, vjust = 1, size = KEY, colour = c(BLUE, BTC_COL),
           label = c("ADNI-4 CN pooled", "BTC controls pooled (dashed)")) +
  annotate("point", x = XLIM[1] * 1.07, y = ytop_c * 0.828, shape = 6, size = 1.4, colour = BLUE, stroke = .5) +
  annotate("text", x = XLIM[1] * 1.15, y = ytop_c * 0.85, hjust = 0, vjust = 1, size = KEY, colour = BLUE, lineheight = .85,
           label = sprintf("KDE mode / shoulder\n(log bw %s)", bw)) +
  scale_x_log10(breaks = c(3, 5, 10, 20, 50, 100), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(xlim = XLIM, ylim = c(0, ytop_c)) +
  labs(x = "Euclidean connection length (mm, log scale)", y = "weight density per unit log-length") + theme_fig()

# ---- d: SIFT2-weighted vs. connection-count density
d <- linear(adni); med <- v[c("median_weighted_mm", "median_edge_count_mm")]
ymax_d <- max(d$w_density, d$e_density); ytop_d <- ymax_d * 1.1
pD <- ggplot(d, aes(bin_center_mm)) +
  geom_line(aes(y = e_density), colour = "grey40", linewidth = .5) +
  geom_line(aes(y = w_density), colour = BLUE, linewidth = .6) +
  annotate("text", x = d$bin_center_mm[which.max(d$w_density)] + 2.5, y = ymax_d * 1.03, label = "SIFT2-weighted", hjust = 0, size = KEY, colour = BLUE) +
  annotate("text", x = 40, y = approx(d$bin_center_mm, d$e_density, 40)$y + 0.08 * ymax_d, label = "edge count (unweighted)", hjust = 0, vjust = 0, size = KEY, colour = "grey40") +
  annotate("segment", x = med, xend = med, y = 0, yend = 0.03 * ytop_d, colour = c(BLUE, "grey40"), linewidth = .6) +
  annotate("text", x = med + c(-0.8, 0.8), y = 0.065 * ytop_d, label = sprintf("median %s", fmt(med)), hjust = c(1, 0), vjust = 0, size = SMALL, colour = c(BLUE, "grey40")) +
  scale_x_continuous(breaks = seq(0, XMAX, 15), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(xlim = c(0, XMAX), ylim = c(0, ytop_d)) +
  labs(x = "Euclidean connection length (mm)", y = "fraction per mm") + theme_fig()

top <- (pA | pB) + plot_layout(widths = c(2.05, 1))
bot <- (pC | pD) + plot_layout(widths = c(1.1, 1))
fig <- tag_panels(top / bot + plot_layout(heights = c(1, 0.95)))
save_figure(fig, "fig3_distance_bands", 180, 105)
