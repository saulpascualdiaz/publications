# Fig. 6 — An individually referenced short-to-long balance captures cognitive and affective
# variation across the Alzheimer's continuum (short-band normative deviation z, ADNI-4).
#   a  group contrasts (patients vs. CN; AD vs. CN): t within clusters significant at pFWE < .05
#   b  pooled vertex-wise models of episodic memory, CDR-Sum of Boxes and GDS
#   c  cortex-wide partial r per outcome under the pooled, diagnosis-adjusted, within-patient and within-CN models
#   d  spatial correlation of each outcome's map with the patients-vs-CN map (spin test)
#   e  split-half reliability of the deviation maps by diagnostic group
source("scripts/s00_setup.R")
source("scripts/s00_surface.R")

v           <- read_values("clinical/values.csv")
results     <- read_data("clinical/results.csv")
contrasts   <- read_data("clinical/contrasts.csv")
spin        <- read_data("clinical/spin.csv")
reliability <- read_data("clinical/reliability.csv")
maps        <- read_data("clinical/vertex_maps.csv") %>% split(.$map)

th <- theme_fig() + theme(panel.grid.major.y = element_blank())
OUTCOMES <- c(PHC_MEM = "episodic memory", R_mem = "memory reserve", CDRSB = "CDR-SB", GDS = "GDS", PSS = "PSS")
MODELS   <- c(m1_pooled = "pooled", m2_pooled_dx = "diagnosis-adjusted", m3_patients = "within patients", m4_cn = "within CN")
zero_line <- geom_vline(xintercept = 0, colour = "grey55", linewidth = .35)
minus     <- label_number(style_negative = "minus", accuracy = 0.1)

# ---- a, b: cluster-corrected maps
MAPS_A <- c(patients_vs_CN = "patients (MCI + AD) vs CN", AD_vs_CN = "AD vs CN")
MAPS_B <- c(PHC_MEM = "episodic memory (PHC_MEM), pooled", CDRSB = "CDR-SB, pooled", GDS = "GDS, pooled")
map_label <- function(d, labels) {
  d <- d[match(names(labels), d$outcome), ]
  sprintf("%s, n = %d: %d cluster%s at pFWE < .05, top pFWE = %s", labels, d$n, d$n_clusters_fwe05, ifelse(d$n_clusters_fwe05 == 1, "", "s"), fmt_p(d$pfwe_cluster))
}
rng   <- t_range(maps, v["t_cluster_forming"])
tiles <- lapply(maps[c(names(MAPS_A), names(MAPS_B))], render_map, t_min = v["t_cluster_forming"], t_max = rng$max)
size  <- tile_size(tiles)
pA <- brain_panel(tiles[names(MAPS_A)], map_label(contrasts, MAPS_A), size)
pB <- brain_panel(tiles[names(MAPS_B)], map_label(filter(results, model == "m1_pooled"), MAPS_B), size)
pBar <- colour_bar(v["t_cluster_forming"], rng$max, rng$two_sided,
                   "t within clusters at pFWE < .05 (grey: n.s.)\ncovariates: age, sex, education, site, coverage")

# ---- c: cortex-wide partial r under four models
c_ <- results %>% filter(outcome %in% names(OUTCOMES), model %in% names(MODELS)) %>%
  mutate(y = length(OUTCOMES) - match(outcome, names(OUTCOMES)) + 1 + (2.5 - match(model, names(MODELS))) * 0.19,
         col = MODEL_COLS[model], fill = ifelse(pfwe_cluster < ALPHA, col, "white"))
n_range <- c_ %>% group_by(model) %>% summarise(n = ifelse(min(n) == max(n), as.character(min(n)), paste0(min(n), "–", max(n))))
pC <- ggplot(c_) + zero_line +
  geom_segment(aes(x = global_r_ci_lo, xend = global_r_ci_hi, y = y, yend = y, colour = col), linewidth = .4) +
  geom_point(aes(global_partial_r, y, colour = col, fill = fill), shape = 21, size = 1.5, stroke = .5) +
  scale_colour_identity() + scale_fill_identity() +
  annotate("text", x = c(-0.40, 0.03, -0.40, 0.03), y = c(6.5, 6.5, 6.17, 6.17), hjust = 0, vjust = 1, size = SMALL, colour = MODEL_COLS,
           label = sprintf("%s, n %s", MODELS, n_range$n[match(names(MODELS), n_range$model)])) +
  annotate("text", x = -0.40, y = 5.84, label = "filled: cluster pFWE < .05", hjust = 0, vjust = 1, size = SMALL, colour = "grey30") +
  scale_y_continuous(breaks = length(OUTCOMES):1, labels = OUTCOMES) + scale_x_continuous(breaks = seq(-0.4, 0.4, 0.2), labels = minus) +
  coord_cartesian(xlim = c(-0.40, 0.50), ylim = c(0.5, 6.55)) +
  labs(x = "cortex-wide partial r (95 % CI)", y = NULL) + th

# ---- d: spin correlation with the patients-vs-CN map
d_ <- spin %>% filter(outcome %in% names(OUTCOMES), model %in% c("m1_pooled", "m2_pooled_dx")) %>%
  mutate(row = length(OUTCOMES) - match(outcome, names(OUTCOMES)) + 1, pooled = model == "m1_pooled", y = row + ifelse(pooled, 0.17, -0.17),
         col = MODEL_COLS[model], fill = ifelse(p_spin < ALPHA, col, "white"))
null_d <- max(d_$null_p95_abs)
pD <- ggplot(d_) +
  annotate("rect", xmin = -null_d, xmax = null_d, ymin = -Inf, ymax = Inf, fill = "grey93") + zero_line +
  geom_point(aes(r, y, colour = col, fill = fill), shape = 21, size = 1.5, stroke = .5) +
  geom_text(aes(ifelse(pooled, 1.10, 1.40), row, label = fmt_p(p_spin), colour = col), hjust = 0, size = SMALL) +
  scale_colour_identity() + scale_fill_identity() +
  annotate("text", x = 1.10, y = 5.78, label = "p spin", hjust = 0, vjust = 1, size = SMALL, colour = "grey30") +
  annotate("text", x = -0.98, y = 5.78, label = "pooled (blue), diagnosis-adjusted (black)", hjust = 0, vjust = 1, size = KEY, colour = "grey30") +
  annotate("text", x = -0.98, y = 0.45, label = "filled: p spin < .05", hjust = 0, vjust = 0, size = SMALL, colour = "grey30") +
  annotate("text", x = 0.26, y = 0.45, label = sprintf("spin null, 95 %% |r| ≤ %s", fmt_r(null_d)), hjust = 0, vjust = 0, size = SMALL, colour = "grey45") +
  scale_y_continuous(breaks = length(OUTCOMES):1, labels = OUTCOMES) + scale_x_continuous(breaks = seq(-1, 1, 0.5), labels = minus) +
  coord_cartesian(xlim = c(-1, 1.62), ylim = c(0.4, 5.85)) +
  labs(x = sprintf("spatial r with the patients-vs-CN map (%s rotations)", format(v["n_rotations"], big.mark = ",")), y = NULL) + th

# ---- e: split-half reliability of the deviation maps by group
GROUPS <- c("CN", "MCI", "AD")
e_ <- reliability %>%
  pivot_longer(c(within_r_mean, between_r_mean), names_to = "which", values_to = "r") %>%
  mutate(within = which == "within_r_mean", sd = ifelse(within, within_r_sd, between_r_sd), smoothed = variant != "unsmoothed",
         y = length(GROUPS) - match(group, GROUPS) + 1 + ifelse(smoothed, -0.17, 0.17), col = ifelse(within, RED, BLUE), fill = ifelse(smoothed, "white", col))
n_group <- reliability %>% distinct(group, n) %>% arrange(match(group, GROUPS))
pE <- ggplot(e_) + zero_line +
  geom_segment(aes(x = r - sd, xend = r + sd, y = y, yend = y, colour = col), linewidth = .35) +
  geom_point(aes(r, y, colour = col, fill = fill), shape = 21, size = 1.5, stroke = .5) +
  scale_colour_identity() + scale_fill_identity() +
  annotate("text", x = -0.19, y = c(4.08, 3.74), hjust = 0, vjust = 1, size = KEY, colour = "grey30",
           label = c("within-participant (red), between-participant (blue)", sprintf("filled: unsmoothed, open: smoothed (k = %d)", v["k_smooth"]))) +
  annotate("text", x = -0.19, y = 0.42, hjust = 0, vjust = 0, size = KEY, colour = "grey30",
           label = sprintf("identification %d %% (both directions, every group)", as.integer(100 * min(reliability$ident_AtoB, reliability$ident_BtoA)))) +
  scale_y_continuous(breaks = length(GROUPS):1, labels = sprintf("%s (n = %d)", n_group$group, n_group$n)) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = label_number(style_negative = "minus", accuracy = 0.01)) +
  coord_cartesian(xlim = c(-0.2, 1.05), ylim = c(0.35, 4.15)) +
  labs(x = "split-half similarity of the deviation maps, r (mean ± SD)", y = NULL) + th

# ---- assemble: the left column holds the fixed-aspect brain panels, so the figure height follows from them
W <- 180; LEFT <- 0.56; H_BAR <- 0.05
hA <- attr(pA, "aspect"); hB <- attr(pB, "aspect")
bar   <- (plot_spacer() | pBar | plot_spacer()) + plot_layout(widths = c(0.2, 0.6, 0.2))
left  <- (pA / pB / bar) + plot_layout(heights = c(hA, hB, H_BAR))
right <- (pC / pD / pE) + plot_layout(heights = c(1.4, 1.15, 0.9))
fig <- tag_panels((left | right) + plot_layout(widths = c(LEFT, 1 - LEFT)), c("a", "b", "", "c", "d", "e"))
save_figure(fig, "fig6_clinical_balance", W, (W * LEFT - 7.5) * (hA + hB + H_BAR) + 15)
