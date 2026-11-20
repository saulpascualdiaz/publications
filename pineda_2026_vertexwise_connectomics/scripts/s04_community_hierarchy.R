# Fig. 4 — A nested community hierarchy whose modules become shorter-range and more compact
# with increasing resolution (Leiden modularity, seven resolutions, ADNI-4 CN).
#   a  number of modules per participant
#   b  within-module connection length, module radius of gyration and within-module connectivity fraction
#   c  containment of each level's modules within the preceding level
#   d  hemispheric purity of the modules
#   e  fraction of native vertices assigned to a module vs. cortical coverage
source("scripts/s00_setup.R")

levels_  <- read_data("hierarchy/level_summary.csv")            # cohort mean and bootstrap 95 % CI per level
subj     <- read_data("hierarchy/subject_level_metrics.csv")    # one row per participant and level
subj_sum <- read_data("hierarchy/subject_summary.csv")          # one row per participant
TERT     <- c("low", "mid", "high")                             # cortical-coverage tertiles
tertiles <- read_data("hierarchy/coverage_tertile_curves.csv") %>% mutate(tertile = factor(tertile, TERT))
contain  <- read_data("hierarchy/containment_transitions.csv")
hemi     <- read_data("hierarchy/hemisphere_summary.csv")

PURITY_MIN <- 95                              # a module is hemisphere-pure when >= 95 % of it lies in one hemisphere
GAMMAS  <- sort(unique(subj$resolution))
N_LAB   <- sprintf("n = %d", n_distinct(subj$subject))
XLIM_G  <- c(0.044, 5.7)
x_gamma <- scale_x_log10(breaks = GAMMAS, labels = format(GAMMAS, drop0trailing = TRUE, trim = TRUE))
th <- theme_fig()
subject_lines <- function(y) geom_line(data = subj, aes(resolution, .data[[y]], group = subject), colour = "grey74", linewidth = .15, alpha = .55)
level_of <- function(m) filter(levels_, metric == m)

# ---- a: module count per level
la <- level_of("n_modules"); ymax_a <- max(subj$n_modules)
pA <- ggplot(la, aes(resolution)) + subject_lines("n_modules") +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, colour = BLUE, linewidth = .4) +
  geom_line(aes(y = mean), colour = BLUE, linewidth = .7) + geom_point(aes(y = mean), colour = BLUE, size = 1.3) +
  annotate("text", x = XLIM_G[1] * 1.08, y = ymax_a * c(1.22, 0.98), hjust = 0, size = KEY, colour = c(BLUE, "grey40"),
           label = c(sprintf("mean %s → %s, 95 %% CI", fmt(min(la$mean)), fmt(max(la$mean))), sprintf("per participant, %s", N_LAB))) +
  x_gamma + scale_y_log10(breaks = c(2, 5, 10, 20, 50, 100)) +
  coord_cartesian(xlim = XLIM_G, ylim = c(1.8, ymax_a * 1.35)) +
  labs(x = "Leiden resolution γ", y = "modules per participant") + th

# ---- b: three level curves, with the means of the cortical-coverage tertiles
TERT_COLS <- setNames(BLUE_RAMP(5)[3:5], TERT)
TERT_LT   <- c(low = "solid", mid = "31", high = "11")
tert_lab  <- tertiles %>% distinct(tertile, cov_min, cov_max) %>% arrange(tertile) %>%
  with(sprintf("%s %.0f–%.0f %%", tertile, cov_min, cov_max))
level_panel <- function(m, ylab, pad, key = FALSE) {
  ylim <- range(subj[[m]]) + c(-pad, pad)
  p <- ggplot(level_of(m), aes(resolution)) + subject_lines(m) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, colour = BLUE, linewidth = .4) +
    geom_line(aes(y = mean), colour = BLUE, linewidth = .8) +
    geom_line(data = filter(tertiles, metric == m), aes(y = mean, colour = tertile, linetype = tertile), linewidth = .45) +
    scale_colour_manual(values = TERT_COLS) + scale_linetype_manual(values = TERT_LT) +
    x_gamma + coord_cartesian(xlim = XLIM_G, ylim = ylim) + labs(x = "Leiden resolution γ", y = ylab) + th
  if (!key) return(p)
  yk <- ylim[2] - diff(ylim) * c(0.03, 0.10, 0.17, 0.24, 0.31)
  p + annotate("segment", x = 0.23, xend = 0.36, y = yk[1], yend = yk[1], colour = BLUE, linewidth = .8) +
    annotate("text", x = 0.40, y = yk[1:2], label = c("mean, 95 % CI", "coverage tertiles:"), hjust = 0, size = KEY, colour = c(BLUE, "grey30")) +
    annotate("segment", x = 0.23, xend = 0.36, y = yk[3:5], yend = yk[3:5], colour = TERT_COLS, linetype = TERT_LT, linewidth = .45) +
    annotate("text", x = 0.40, y = yk[3:5], label = tert_lab, hjust = 0, size = KEY, colour = TERT_COLS)
}
pB1 <- level_panel("mean_within_len_mm", "SIFT2-weighted within-module\nconnection length (mm)", 0.5)
pB2 <- level_panel("mean_compactness_mm", "module radius of gyration (mm)", 1, key = TRUE)
pB3 <- level_panel("frac_within_edges", "within-module connectivity fraction", 0.02)

# ---- c: containment in the previous level, per transition
steps <- contain %>% filter(level_order > 0) %>%
  mutate(label = sprintf("γ %s\n→ %s", format(gamma_from, drop0trailing = TRUE, trim = TRUE), format(gamma_to, drop0trailing = TRUE, trim = TRUE)))
subj_c <- filter(subj, level_order >= 1)
mean_all  <- contain$mean[contain$transition == "mean_over_6_transitions"]
mean_rest <- contain$mean[contain$transition == "mean_over_5_transitions_excl_L0toL1"]
YLIM_C <- c(min(subj_c$containment_prev) - 0.03, 1.03)
pC <- ggplot(steps, aes(level_order)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_point(data = subj_c, aes(y = containment_prev), colour = "grey74", size = .45, alpha = .7,
             position = position_jitter(width = .22, height = 0, seed = 1)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, colour = BLUE, linewidth = .5) +
  geom_line(aes(y = mean), colour = BLUE, linewidth = .5) + geom_point(aes(y = mean), colour = BLUE, size = 1.5) +
  annotate("text", x = 0.58, y = 1.016, label = "hemispheric split", hjust = 0, size = KEY, colour = "grey30") +
  annotate("text", x = 6.45, y = YLIM_C[1] + 0.012 + c(0.048, 0.032, 0.016, 0), hjust = 1, vjust = 0, size = KEY,
           colour = c("grey30", "grey30", BLUE, "grey40"),
           label = c(sprintf("all six transitions, mean %s", fmt(mean_all, 3)), sprintf("excluding the first, mean %s", fmt(mean_rest, 3)),
                     "mean, 95 % CI", sprintf("per participant, %s", N_LAB))) +
  scale_x_continuous(breaks = 1:6, labels = steps$label, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.75, 1, 0.05)) +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = YLIM_C) +
  labs(x = "transition between consecutive levels", y = "containment in the previous level") + th +
  theme(axis.text.x = element_text(lineheight = .85), panel.grid.major.x = element_blank())

# ---- d: hemispheric purity per level
purity <- hemi %>% select(resolution, mean = size_weighted_purity_mean, median = size_weighted_purity_median, min = size_weighted_purity_min) %>%
  pivot_longer(-resolution, names_to = "stat", values_to = "purity")
PUR_LT <- c(mean = "solid", median = "31", min = "11")
counts0 <- table(subj_sum$level0_n_modules)
level0_lab <- sprintf("level 0 (γ %s), modules per participant:\n%s\ntwo ≥ %d %%-pure hemispheres in %d of %d", GAMMAS[1],
                      paste(sprintf("%s in %d", names(counts0), counts0), collapse = ", "), PURITY_MIN,
                      sum(subj_sum$level0_modules_are_hemispheres), nrow(subj_sum))
yk_d <- 0.805 - c(0, 0.032, 0.064, 0.096, 0.128, 0.158)
pD <- ggplot(mapping = aes(resolution)) +
  geom_line(data = purity, aes(y = purity, linetype = stat), colour = BLUE, linewidth = .55) +
  geom_point(data = filter(purity, stat == "mean"), aes(y = purity), colour = BLUE, size = 1.3) +
  geom_line(data = hemi, aes(y = frac_assigned_vtx_in_pure95_mean), colour = "grey30", linewidth = .55) +
  geom_point(data = hemi, aes(y = frac_assigned_vtx_in_pure95_mean), shape = 21, fill = "white", colour = "grey30", size = 1.3, stroke = .5) +
  scale_linetype_manual(values = PUR_LT) +
  annotate("segment", x = 0.19, xend = 0.28, y = yk_d[2:4], yend = yk_d[2:4], colour = BLUE, linetype = PUR_LT, linewidth = .55) +
  annotate("point", x = 0.235, y = yk_d[2], colour = BLUE, size = 1.3) +
  annotate("point", x = 0.235, y = yk_d[5], shape = 21, fill = "white", colour = "grey30", size = 1.3, stroke = .5) +
  annotate("text", x = 0.31, y = yk_d, hjust = 0, size = KEY, colour = rep(c(BLUE, "grey30"), c(4, 2)),
           label = c("hemispheric purity of the modules", "size-weighted, mean", "median", "min",
                     "assigned vertices in", sprintf("≥ %d %%-pure modules, mean", PURITY_MIN))) +
  annotate("text", x = XLIM_G[2] * 0.97, y = 0.508, hjust = 1, vjust = 0, size = SMALL, colour = "grey30", lineheight = .95, label = level0_lab) +
  x_gamma + scale_y_continuous(breaks = seq(0.5, 1, 0.1)) +
  coord_cartesian(xlim = XLIM_G, ylim = c(0.5, 1)) +
  labs(x = "Leiden resolution γ", y = "hemispheric purity / fraction") + th

# ---- e: assigned vertices vs. cortical coverage
rho_e  <- cor(subj_sum$coverage_pct, subj_sum$frac_assigned, method = "spearman")
ols    <- coef(lm(frac_assigned ~ coverage_pct, subj_sum))
ratio  <- subj_sum$frac_assigned_over_coverage
finest <- filter(subj, level_order == max(level_order))
rho_n  <- cor.test(finest$coverage_pct, finest$n_modules, method = "spearman", exact = FALSE)
sci <- function(p) {                                    # 9.7e-09 -> "9.7 × 10⁻⁹"
  e <- floor(log10(p)); sup <- c("⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹")
  sprintf("%.1f × 10%s%s", p / 10^e, if (e < 0) "⁻" else "", paste(sup[as.integer(strsplit(as.character(abs(e)), "")[[1]]) + 1], collapse = ""))
}
XLIM_E <- range(subj_sum$coverage_pct) + c(-3, 3); YLIM_E <- range(subj_sum$frac_assigned) + c(-0.04, 0.05)
rows_e <- diff(YLIM_E) * c(0.02, 0.085, 0.15, 0.215)
pE <- ggplot(subj_sum, aes(coverage_pct, frac_assigned)) +
  annotate("line", x = XLIM_E, y = ols[1] + ols[2] * XLIM_E, colour = "grey45", linetype = "22", linewidth = .4) +
  geom_point(shape = 21, fill = PALE_BLUE, colour = "grey35", stroke = .25, size = 1.1, alpha = .9) +
  annotate("text", x = XLIM_E[1] + 0.5, y = YLIM_E[2] - rows_e, hjust = 0, size = KEY, colour = "grey30",
           label = c(sprintf("Spearman ρ = %s, %s", fmt(rho_e, 3), N_LAB), sprintf("OLS %s per 10 %% (dashed)", fmt(10 * ols[2], 3)),
                     "assigned ÷ coverage:", sprintf("mean %s (%s–%s)", fmt(mean(ratio), 2), fmt(min(ratio), 2), fmt(max(ratio), 2)))) +
  annotate("text", x = XLIM_E[2] - 0.5, y = YLIM_E[1] + rev(rows_e), hjust = 1, size = KEY, colour = "grey30",
           label = c("coverage vs", sprintf("module count at γ %s:", max(GAMMAS)), sprintf("Spearman ρ = %s", fmt(rho_n$estimate, 2)),
                     sprintf("p = %s", sci(rho_n$p.value)))) +
  scale_x_continuous(breaks = seq(40, 100, 10)) + scale_y_continuous(breaks = seq(0.4, 0.9, 0.1)) +
  coord_cartesian(xlim = XLIM_E, ylim = YLIM_E) +
  labs(x = "coverage: % of native vertices\nwith ≥ 1 streamline endpoint", y = "fraction of native vertices\nassigned to a module") + th

top <- (pA | pB1 | pB2 | pB3) + plot_layout(widths = c(0.9, 1.05, 1.05, 1.05))
bot <- (pC | pD | pE) + plot_layout(widths = c(1.1, 1.2, 0.85))
fig <- tag_panels(top / bot, c("a", "b", "", "", "c", "d", "e"))
save_figure(fig, "fig4_community_hierarchy", 180, 110)
