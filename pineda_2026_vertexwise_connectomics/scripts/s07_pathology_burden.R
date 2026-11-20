# Fig. 7 — The individually referenced short-to-long balance tracks Alzheimer's pathology burden
# beyond diagnosis (composite = first principal component of five markers, complete cases).
#   a  composite maps, pooled and diagnosis-adjusted: t within clusters significant at pFWE < .05
#   b  cortex-wide partial r for the composite and each marker, on the common sample
#   c  leave-one-site-out analysis of the pooled composite association
#   d  composite association under different covariate sets, and with plasma within two years of imaging
#   e  composite construction: marker loadings and component scores by diagnostic group
source("scripts/s00_setup.R")
source("scripts/s00_surface.R")

v        <- read_values("pathology/values.csv")
results  <- read_data("pathology/results.csv")
pairs    <- read_data("pathology/marker_comparison.csv")
loso     <- read_data("pathology/leave_one_site_out.csv", col_types = cols(dropped_site = col_character()))
loadings <- read_data("pathology/pca_loadings.csv")
scores   <- read_data("pathology/composite_scores.csv")
maps     <- read_data("pathology/vertex_maps.csv") %>% split(.$map)

th <- theme_fig() + theme(panel.grid.major.y = element_blank())
MARKERS   <- c(GFAP = "GFAP (plasma)", NfL = "NfL (plasma)", pT217 = "p-tau217 (plasma)", tauPET = "tau-PET", amyloid = "amyloid-PET")
zero_line <- geom_vline(xintercept = 0, colour = "grey55", linewidth = .35)
minus     <- label_number(style_negative = "minus", accuracy = 0.1)
composite <- filter(results, analysis == "composite")
pooled    <- filter(composite, model == "m1_pooled")
N         <- pooled$n
forest <- function(d) list(                                            # CI segment + point, filled when significant
  geom_segment(data = d, aes(x = global_r_ci_lo, xend = global_r_ci_hi, y = y, yend = y, colour = col), linewidth = .4),
  geom_point(data = d, aes(global_partial_r, y, colour = col, fill = ifelse(significant, col, "white")), shape = 21, size = 1.5, stroke = .5),
  scale_colour_identity(), scale_fill_identity())

# ---- a: composite maps
MAPS <- c(pooled = "composite (PC1), pooled (m1)", diagnosis_adjusted = "diagnosis-adjusted (m2)")
rng   <- t_range(maps, v["t_cluster_forming"])
tiles <- lapply(maps[names(MAPS)], render_map, t_min = v["t_cluster_forming"], t_max = rng$max)
pA <- brain_panel(tiles, sprintf("%s, n = %d: %d cluster%s at pFWE < .05, top pFWE = %s", MAPS, composite$n, composite$n_clusters_fwe05,
                                 ifelse(composite$n_clusters_fwe05 == 1, "", "s"), fmt_p(composite$pfwe_cluster)), tile_size(tiles))
pBar <- colour_bar(v["t_cluster_forming"], rng$max, rng$two_sided,
                   "t within clusters at pFWE < .05 (grey: n.s.)\ncovariates: age, sex, education, site, coverage")

# ---- b: composite and markers on the common sample (markers: Holm-corrected over the five, within model)
ROWS_B <- c(composite = "composite (PC1)", MARKERS)
b_ <- results %>% filter(analysis %in% c("composite", "marker")) %>%
  mutate(y = length(ROWS_B) - match(outcome, names(ROWS_B)) + 1 + ifelse(model == "m1_pooled", 0.17, -0.17),
         col = MODEL_COLS[model], significant = coalesce(holm_p_cluster, pfwe_cluster) < ALPHA)
steiger <- pairs %>% group_by(model) %>% summarise(p = min(steiger_p), n = n())
pB <- ggplot() + zero_line + forest(b_) +
  annotate("text", x = -0.15, y = c(8.2, 7.72, 7.24), hjust = 0, vjust = 1, size = SMALL, colour = "grey30",
           label = c(sprintf("pooled, n = %d (blue); diagnosis-adjusted (black)", N), "filled: cluster pFWE < .05", "(markers: Holm over five, within model)")) +
  annotate("text", x = -0.15, y = 0.15, hjust = 0, vjust = 1, size = SMALL, colour = "grey30", lineheight = .95,
           label = sprintf("all %d pairwise marker differences n.s.\n(Steiger p ≥ %s pooled, ≥ %s diagnosis-adjusted)", steiger$n[1],
                           fmt_r(steiger$p[steiger$model == "m1_pooled"]), fmt_r(steiger$p[steiger$model == "m2_pooled_dx"]))) +
  scale_y_continuous(breaks = length(ROWS_B):1, labels = ROWS_B) + scale_x_continuous(breaks = seq(0, 0.4, 0.2), labels = minus) +
  coord_cartesian(xlim = c(-0.15, 0.50), ylim = c(-0.8, 8.3)) +
  labs(x = "cortex-wide partial r (95 % CI), common sample", y = NULL) + th

# ---- c: leave one site out
c_ <- loso %>% arrange(desc(n_dropped), dropped_site) %>% mutate(x = row_number(), significant = pfwe_cluster < ALPHA)
pC <- ggplot(c_, aes(x)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = pooled$global_r_ci_lo, ymax = pooled$global_r_ci_hi, fill = "grey93") +
  geom_hline(yintercept = pooled$global_partial_r, colour = "grey45", linewidth = .4, linetype = "22") +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = .35) +
  geom_segment(aes(xend = x, y = global_r_ci_lo, yend = global_r_ci_hi), colour = BLUE, linewidth = .4) +
  geom_point(aes(y = global_partial_r, fill = ifelse(significant, BLUE, "white")), shape = 21, colour = BLUE, size = 1.5, stroke = .5) +
  scale_fill_identity() +
  annotate("text", x = 0.6, y = c(0.595, 0.54), hjust = 0, vjust = 1, size = KEY, colour = "grey30",
           label = c(sprintf("full sample r = %s [%s, %s], n = %d", fmt_r(pooled$global_partial_r), fmt_r(pooled$global_r_ci_lo), fmt_r(pooled$global_r_ci_hi), N),
                     sprintf("filled: top-cluster pFWE < .05 (%d/%d sites); %s permutations per drop", sum(c_$significant), nrow(c_),
                             format(v["n_permutations_loso"], big.mark = ",")))) +
  scale_x_continuous(breaks = c_$x, labels = sprintf("%s\n%d", c_$dropped_site, c_$n_dropped), expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.2), labels = minus) + coord_cartesian(ylim = c(-0.05, 0.60)) +
  labs(x = "site dropped (upper), participants removed (lower); pooled model", y = "cortex-wide partial r (95 % CI)") +
  theme_fig() + theme(axis.text.x = element_text(size = 5.5, lineheight = .85), panel.grid.major.x = element_blank())

# ---- d: dependence of the composite association on the covariates
ROWS_D <- c(cont_nocov = "no covariates", cont_cov_nosite = "age, sex, education, coverage", cont_siteonly = "site only",
            m1_pooled = "full (m1): + site", m2_pooled_dx = "diagnosis-adjusted (m2)", plasma_within_2y = "|gap| ≤ 2 y, full (m1)")
d_ <- results %>% filter(analysis != "marker") %>%
  mutate(key = ifelse(analysis == "plasma_within_2y", analysis, model), y = length(ROWS_D) - match(key, names(ROWS_D)) + 1,
         col = ifelse(key == "m2_pooled_dx", MODEL_COLS[["m2_pooled_dx"]], BLUE), significant = pfwe_cluster < ALPHA)
pD <- ggplot() + zero_line + forest(d_) +
  geom_text(data = d_, aes(0.56, y, label = sprintf("n = %d", n)), hjust = 0, size = SMALL, colour = "grey30") +
  annotate("text", x = -0.29, y = 0.55, label = "filled: cluster pFWE < .05", hjust = 0, vjust = 0, size = KEY, colour = "grey30") +
  scale_y_continuous(breaks = length(ROWS_D):1, labels = ROWS_D) + scale_x_continuous(breaks = seq(-0.2, 0.4, 0.2), labels = minus) +
  coord_cartesian(xlim = c(-0.30, 0.72), ylim = c(0.5, 6.5)) +
  labs(x = "cortex-wide partial r (95 % CI), composite", y = NULL) + th

# ---- e: composite construction
e1 <- loadings %>% mutate(y = rev(row_number()))
pE1 <- ggplot(e1, aes(loading, y)) +
  geom_col(fill = PALE_BLUE, colour = "grey35", linewidth = .25, width = .62, orientation = "y") +
  geom_text(aes(loading + 0.012, label = fmt(loading, 2)), hjust = 0, size = KEY, colour = "grey25") +
  annotate("text", x = 0.005, y = 6.5, hjust = 0, vjust = 1, size = SMALL, colour = "grey30", lineheight = .9,
           label = sprintf("PC1: %.1f %% of the variance\n(five z-scored markers, n = %d)", 100 * v["pc1_variance_explained"], nrow(scores))) +
  scale_y_continuous(breaks = e1$y, labels = MARKERS[e1$marker]) + scale_x_continuous(breaks = seq(0, 0.6, 0.2), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 0.72), ylim = c(0.4, 6.55)) +
  labs(x = "PC1 loading", y = NULL) + th
e2 <- scores %>% mutate(group = factor(group, c("CN", "MCI", "AD")))
group_mean <- e2 %>% group_by(group) %>% summarise(mean = mean(composite), n = n())
yr_e <- range(e2$composite); ytop_e <- yr_e[2] + 0.3 * diff(yr_e)
pE2 <- ggplot(e2, aes(group, composite)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .3) +
  geom_point(aes(fill = group), shape = 21, colour = "grey35", stroke = .15, size = .8, alpha = .8, position = position_jitter(width = .22, seed = 1)) +
  geom_boxplot(alpha = 0, colour = "grey15", width = .5, linewidth = .3, coef = 0, outlier.shape = NA) +
  geom_text(data = group_mean, aes(y = ytop_e, label = sprintf("mean\n%s", fmt_r(mean))), vjust = 1, size = SMALL, colour = "grey30", lineheight = .9) +
  scale_fill_manual(values = GROUP_COLS) + scale_x_discrete(labels = sprintf("%s\nn = %d", group_mean$group, group_mean$n)) +
  coord_cartesian(ylim = c(yr_e[1] - 0.05 * diff(yr_e), ytop_e)) +
  labs(x = NULL, y = "composite (PC1 score)") + theme_fig() + theme(axis.text.x = element_text(lineheight = .85))

# ---- assemble: row heights in mm; the first row follows the aspect of the brain panel
W <- 180; LEFT <- 0.56
heights <- c((W * LEFT - 7.5) * attr(pA, "aspect"), 11, 34, 36)
top <- (pA | pB) + plot_layout(widths = c(LEFT, 1 - LEFT))
bar <- (plot_spacer() | pBar | plot_spacer()) + plot_layout(widths = c(LEFT * 0.2, LEFT * 0.6, 1 - LEFT * 0.8))
bot <- (pD | pE1 | pE2) + plot_layout(widths = c(0.44, 0.28, 0.28))
fig <- tag_panels(free(top) / free(bar) / free(pC) / free(bot) + plot_layout(heights = heights), c("a", "b", "", "c", "d", "e", ""))
save_figure(fig, "fig7_pathology_burden", W, sum(heights) + 30)
