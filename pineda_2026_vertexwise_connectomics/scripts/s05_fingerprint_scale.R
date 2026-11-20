# Fig. 5 — The individual fingerprint across spatial scales and atlases.
#   a  Idiff of the band-strength maps (matched split halves) from the vertex scale through
#      k-NN smoothing kernels to atlas parcels, with within- and between-participant similarity (ADNI-4 CN)
#   b  the same for the short band in the BTC controls, before and after removing cortical geometry
#   c  sensitivity to misalignment (BTC): similarity when maps are displaced along the surface
#   d  community-partition fingerprint: within- vs. between-participant AMI per Leiden resolution
#   e  share of SIFT2 weight on intra-parcel connections, per atlas and distance band
#   f  Idiff of different representations of the same atlas
source("scripts/s00_setup.R")

v        <- read_values("fingerprint/values.csv")
band     <- read_data("fingerprint/band_identifiability.csv")
kernels  <- read_data("fingerprint/knn_kernels.csv")
btc      <- read_data("fingerprint/btc_identifiability.csv")
shift    <- read_data("fingerprint/btc_misalignment.csv")
ami      <- read_data("fingerprint/partition_ami.csv")
weight   <- read_data("fingerprint/atlas_weight_share.csv")
atlas_id <- read_data("fingerprint/atlas_identifiability.csv")

th <- theme_fig() + theme(axis.text.x = element_text(lineheight = .85))
ADNI_LAB <- sprintf("ADNI-4 CN, n = %d", unique(band$n_subjects))
BTC_LAB  <- sprintf("BTC controls, n = %d", unique(btc$n))
ATLAS    <- c(glasser = "Glasser-360", dk = "DK-68")
divider  <- function(x) geom_vline(xintercept = x, linetype = "13", colour = "grey60", linewidth = .3)   # kernels | parcellations

# ---- a: Idiff across scale, with the within- and between-participant similarity underneath
RES_A <- c("vertex", paste0("knn", kernels$k), "glasser360", "dk68")
LAB_A <- c(sprintf("vertex\n%.0f mm²", v["vertex_area_mm2"]), sprintf("kNN %d\n%.0f mm²", kernels$k, kernels$area_measured_mm2), "Glasser\n360", "DK\n68")
a <- band %>% mutate(x = match(resolution, RES_A), band = factor(band, BAND_LEVELS))
XLIM_A <- c(0.45, 8.5)
band_scales <- list(scale_colour_manual(values = BAND_COLS), scale_fill_manual(values = BAND_COLS),
                    scale_x_continuous(breaks = 1:8, labels = LAB_A, expand = c(0, 0)))
pA1 <- ggplot(a, aes(x, idiff, colour = band)) + divider(6.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, linewidth = .4) +
  geom_line(linewidth = .5) + geom_point(aes(fill = band), shape = 21, colour = "grey35", stroke = .25, size = 1.4) +
  annotate("point", x = 0.7, y = v["idiff_random_360_vertices"], shape = 21, fill = "white", colour = "grey20", size = 1.7, stroke = .5) +
  annotate("text", x = 0.6, y = 0.99, label = "open: random 360 vertices", hjust = 0, size = KEY, colour = "grey20") +
  annotate("text", x = 0.6, y = 0.32, hjust = 0, vjust = 0, size = KEY, colour = "grey30",
           label = sprintf("identification %d %% at every scale", round(100 * min(a$acc_row, a$acc_col)))) +
  annotate("text", x = 8.45, y = 0.99, label = "parcellations", hjust = 1, size = KEY, colour = "grey40") +
  band_scales + scale_y_continuous(breaks = seq(0.3, 1, 0.1)) +
  coord_cartesian(xlim = XLIM_A, ylim = c(0.29, 1.02)) +
  labs(x = NULL, y = "identifiability (Idiff)") + th + theme(axis.text.x = element_blank())
similarity <- a %>% pivot_longer(c(within_mean, between_mean), names_to = "which", values_to = "r")
key_a <- c(0.90, 0.76, 0.62)
pA2 <- ggplot(similarity, aes(x, r, colour = band, linetype = which)) + divider(6.5) +
  geom_line(linewidth = .45) + geom_point(aes(fill = band), shape = 21, colour = "grey35", stroke = .25, size = 1.1) +
  scale_linetype_manual(values = c(within_mean = "13", between_mean = "solid")) +
  annotate("point", x = 0.85, y = key_a, shape = 22, size = 1.9, fill = BAND_COLS, colour = "grey35", stroke = .2) +
  annotate("text", x = 1.05, y = key_a, label = BAND_LEVELS, hjust = 0, size = KEY, colour = "grey30") +
  annotate("text", x = 8.45, y = c(0.88, 0.76), hjust = 1, size = KEY, colour = "grey30",
           label = c("within-participant (dotted)", "between-participant (solid)")) +
  annotate("text", x = 2.3, y = 0.62, label = ADNI_LAB, hjust = 0, size = KEY, colour = "grey30") +
  band_scales + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(xlim = XLIM_A, ylim = c(0, 1.04)) +
  labs(x = "spatial scale (k-NN kernel area; parcellation)", y = "similarity (r)") + th

# ---- b: BTC controls, raw vs. geometry removed (curvature, sulcal depth, thickness, coverage regressed out)
RES_B <- c(vertex = "vertex", knn10 = "kNN\n10", knn30 = "kNN\n30", knn100 = "kNN\n100", knn150 = "kNN\n150", glasser360 = "Glasser\n360", dk68 = "DK\n68")
b <- btc %>% filter(resolution %in% names(RES_B)) %>%
  mutate(model = factor(model, c("raw", "resid_L")), x = match(resolution, names(RES_B)) + ifelse(model == "raw", -0.12, 0.12)) %>% arrange(model, x)
SHORT <- BAND_COLS[["short"]]
key_b <- c(0.985, 0.925)
pB <- ggplot(b, aes(x, idiff, linetype = model)) + divider(5.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, colour = SHORT, linewidth = .4) +
  geom_line(colour = SHORT, linewidth = .5) + geom_point(aes(fill = model), shape = 21, colour = "grey35", stroke = .3, size = 1.4) +
  scale_linetype_manual(values = c(raw = "solid", resid_L = "22")) + scale_fill_manual(values = c(raw = SHORT, resid_L = "white")) +
  annotate("segment", x = 3.6, xend = 4.1, y = key_b, yend = key_b, colour = SHORT, linewidth = .5, linetype = c("solid", "22")) +
  annotate("point", x = 3.85, y = key_b, shape = 21, fill = c(SHORT, "white"), colour = "grey35", stroke = .3, size = 1.4) +
  annotate("text", x = 4.25, y = key_b, label = c("raw", "geometry removed"), hjust = 0, size = KEY, colour = "grey30") +
  annotate("text", x = 0.6, y = 0.32, hjust = 0, vjust = 0, size = KEY, colour = "grey30", lineheight = .9,
           label = sprintf("%s\nidentification %d/%d at every scale", BTC_LAB, min(b$n_correct_row), unique(b$n))) +
  scale_x_continuous(breaks = 1:7, labels = RES_B, expand = c(0, 0)) + scale_y_continuous(breaks = seq(0.3, 1, 0.1)) +
  coord_cartesian(xlim = c(0.5, 7.5), ylim = c(0.29, 1.02)) +
  labs(x = "spatial scale (short band)", y = "identifiability (Idiff)") + th

# ---- c: misalignment: own map displaced (self-similarity) vs. others' maps displaced (between-participant)
shift <- shift %>% mutate(resolution = factor(resolution, c("vertex", "glasser360")))
own    <- shift %>% filter(protocol == "own") %>% rename(r = self_r_mean)
others <- shift %>% filter(protocol == "others") %>% rename(r = other_r_mean)
far    <- own %>% filter(fsavg_white_mm == max(fsavg_white_mm)) %>%
  mutate(label = sprintf("identification %d %%", round(100 * acc_row_mean)), up = resolution == "vertex", y = r + ifelse(up, 0.09, -0.065))
pC <- ggplot(mapping = aes(fsavg_white_mm, r, linetype = resolution)) +
  geom_line(data = others, colour = "grey50", linewidth = .35) + geom_line(data = own, colour = RED, linewidth = .5) +
  geom_point(data = own, colour = RED, size = 1.2) + geom_point(data = others, colour = "grey50", size = 0.9) +
  geom_text(data = far, aes(fsavg_white_mm - 0.1, y, label = label, vjust = as.numeric(!up)), hjust = 1, size = KEY, colour = RED) +
  scale_linetype_manual(values = c(vertex = "solid", glasser360 = "22")) +
  annotate("text", x = max(shift$fsavg_white_mm), y = c(1.21, 1.15, 1.09, 1.03), hjust = 1, size = KEY, colour = c(RED, "grey50", "grey30", "grey30"),
           label = c("self-similarity, own map displaced", "between-participant, others displaced", "solid: vertex, dashed: Glasser-360", BTC_LAB)) +
  scale_x_continuous(breaks = unique(shift$fsavg_white_mm), labels = function(x) fmt(x), expand = expansion(add = 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) + coord_cartesian(ylim = c(0, 1.24)) +
  labs(x = "displacement on the white surface (mm)", y = "fingerprint similarity (r)") + th

# ---- d: partition fingerprint (adjusted mutual information per Leiden level)
acc_d <- floor(100 * range(pmin(ami$acc_row, ami$acc_col)))
ci_d  <- range(ami$within_ci_lo, ami$within_ci_hi, ami$between_AB_subjboot_lo, ami$between_AB_subjboot_hi)
# y range: 17.5 % headroom for two key lines, and >= 25 % free under the level-0 points for their label
y_lo <- ci_d[1] - 0.03
for (i in 1:4) y_lo <- min(ci_d[1] - 0.03, min(ami$between_AB_subjboot_lo[1], ami$within_ci_lo[1]) - 0.25 * (ci_d[2] - y_lo) / 0.825)
y_lo <- floor(10 * y_lo) / 10
YLIM_D <- c(y_lo, y_lo + (ci_d[2] - y_lo) / 0.825); STEP_D <- 0.07 * diff(YLIM_D); X_D <- max(ami$level_order) + 0.3
pD <- ggplot(ami, aes(level_order)) +
  annotate("rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_errorbar(aes(ymin = between_AB_subjboot_lo, ymax = between_AB_subjboot_hi), width = 0, colour = BLUE, linewidth = .4) +
  geom_line(aes(y = between_AB_mean), colour = BLUE, linewidth = .5) + geom_point(aes(y = between_AB_mean), colour = BLUE, size = 1.4) +
  geom_errorbar(aes(ymin = within_ci_lo, ymax = within_ci_hi), width = 0, colour = RED, linewidth = .4) +
  geom_line(aes(y = within_mean), colour = RED, linewidth = .5) + geom_point(aes(y = within_mean), colour = RED, size = 1.4) +
  annotate("text", x = 0, y = YLIM_D[1] + 0.1 * STEP_D, hjust = .5, vjust = 0, size = SMALL, colour = "grey30", lineheight = .9,
           label = sprintf("hemispheric\nsplit\n(2 modules\nin %d/%d)", v["level0_n_two_modules"], unique(ami$n_subjects))) +
  annotate("text", x = X_D, y = YLIM_D[2] - c(0.5, 1.5) * STEP_D, hjust = 1, size = KEY, colour = c(RED, "grey30"),
           label = c("within-participant (own split halves)", sprintf("identification %d–%d %% across levels", acc_d[1], acc_d[2]))) +
  annotate("text", x = X_D, y = YLIM_D[1] + c(1.5, 0.5) * STEP_D, hjust = 1, size = KEY, colour = c(BLUE, "grey30"),
           label = c("between-participant (all pairs, matched)", sprintf("gap within − between %s–%s", fmt(min(ami$gap_AB), 2), fmt(max(ami$gap_AB), 2)))) +
  scale_x_continuous(breaks = ami$level_order, labels = sprintf("%g\n%.1f", ami$resolution, ami$mean_n_modules), expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + coord_cartesian(ylim = YLIM_D) +
  labs(x = "Leiden γ (upper), mean modules (lower)", y = "partition similarity (AMI)") + th

# ---- e: share of the SIFT2 weight on connections with both endpoints in the same parcel
BAND_E <- c("all", BAND_LEVELS); COL_E <- c(all = "grey78", BAND_COLS)
e <- weight %>% filter(cls == "intra", band %in% BAND_E) %>%
  mutate(band = factor(band, BAND_E), x = match(atlas, names(ATLAS)) + (as.integer(band) - 2.5) * 0.2)
unlabelled <- weight %>% filter(cls == "unlabelled", band == "all") %>% arrange(match(atlas, names(ATLAS)))
key_e <- c(0.88, 0.82, 0.76, 0.70)
pE <- ggplot(e, aes(x, mean)) +
  geom_col(aes(fill = band), width = 0.17, colour = "grey35", linewidth = .25) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.07, colour = "grey30", linewidth = .35) +
  scale_fill_manual(values = COL_E) +
  annotate("point", x = 0.58, y = key_e, shape = 22, size = 1.9, fill = COL_E, colour = "grey35", stroke = .2) +
  annotate("text", x = 0.66, y = key_e, label = c("all connections", paste(BAND_LEVELS, "band")), hjust = 0, size = KEY, colour = "grey30") +
  annotate("text", x = 0.58, y = 0.63, label = paste0(ADNI_LAB, ", mean ± SD"), hjust = 0, size = KEY, colour = "grey30") +
  annotate("text", x = 2.47, y = 0.88, hjust = 1, vjust = 1, size = KEY, colour = "grey30", lineheight = .9,
           label = sprintf("unlabelled endpoints (no parcel):\n%s %% of all weight, Glasser-360\n%s %%, DK-68", fmt(100 * unlabelled$mean[1]), fmt(100 * unlabelled$mean[2]))) +
  scale_x_continuous(breaks = 1:2, labels = ATLAS, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), labels = label_percent(accuracy = 1, suffix = " %")) +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(0, 0.92)) +
  labs(x = "parcellation", y = "intra-parcel share of SIFT2 weight") + th

# ---- f: the same atlas, four representations
READS <- tribble(
  ~family,                      ~band,   ~intra,  ~label,                                       ~fill,                ~shape, ~offset,
  "edges_log1p",                "all",   "cross", "parcel connectome, edges, all connections",  "grey55",             21,     -0.3,
  "parcel_strength",            "all",   "cross", "parcel connectome, node strength, all",      "white",              23,     -0.1,
  "edges_log1p",                "short", "cross", "parcel connectome, edges, short band",       BAND_COLS[["short"]], 21,      0.1,
  "vertex_strength_parcelmean", "short", "na",    "vertex-wise map, parcel means, short band",  "#5E86C0",            24,      0.3)
f <- atlas_id %>% inner_join(READS, by = c("family", "band", "intra")) %>% filter(atlas %in% names(ATLAS)) %>%
  mutate(x = match(atlas, names(ATLAS)) + offset, label = factor(label, READS$label))
YTOP_F <- max(f$idiff_ci_hi) + 0.30; key_f <- YTOP_F - 0.025 - (0:3) * 0.055
stopifnot(all(pmin(f$acc_row, f$acc_col) == 1))          # the note below states 100 % identification
pF <- ggplot(f, aes(x, idiff)) +
  geom_errorbar(aes(ymin = idiff_ci_lo, ymax = idiff_ci_hi), width = 0, colour = "grey35", linewidth = .4) +
  geom_point(aes(fill = label, shape = label), colour = "grey35", stroke = .3, size = 1.8) +
  scale_fill_manual(values = setNames(READS$fill, READS$label)) + scale_shape_manual(values = setNames(READS$shape, READS$label)) +
  annotate("point", x = 1.38, y = key_f, shape = READS$shape, fill = READS$fill, colour = "grey35", stroke = .3, size = 1.8) +
  annotate("text", x = 1.46, y = key_f, label = READS$label, hjust = 0, size = KEY, colour = "grey30") +
  annotate("text", x = 0.55, y = key_f[1], hjust = 0, vjust = 1, size = KEY, colour = "grey30", lineheight = .95,
           label = paste0(ADNI_LAB, "\nmatched halves\nidentification 100 %\nfor every reading")) +
  scale_x_continuous(breaks = 1:2, labels = ATLAS, expand = c(0, 0)) + scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(0, YTOP_F)) +
  labs(x = "parcellation", y = "identifiability (Idiff)") + th

row1 <- ((pA1 / pA2 + plot_layout(heights = c(1.55, 1))) | pB) + plot_layout(widths = c(1.5, 1))
fig <- tag_panels(row1 / (pC | pD) / (pE | pF) + plot_layout(heights = c(1.25, 1, 1)), c("a", "", "b", "c", "d", "e", "f"))
save_figure(fig, "fig5_fingerprint_scale", 180, 170)
