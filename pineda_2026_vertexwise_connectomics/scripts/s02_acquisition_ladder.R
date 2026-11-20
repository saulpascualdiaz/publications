# Fig. 2 — Identifiability of the short-range fingerprint under streamline subsampling and
# acquisition degradation (BTC controls).
#   a  streamline count: Idiff against the full 5 M tractogram, and between two subsamples
#   b  shell sampling: Idiff after removing shells and halving the directions
#   c  similarity to the participant's own reference vs. to other participants' references
#   d  participant-by-participant similarity matrices per acquisition
source("scripts/s00_setup.R")

LEVELS <- c(L0_3shell = "3 shells\n96 dir", L1_2shell = "2 shells\n46 dir",
            L2_1shell = "1 shell\n30 dir", L3_1shell_lowdir = "1 shell\n15 dir")

str_sum <- read_data("acquisition_ladder/streamlines_summary.csv")
str_sub <- read_data("acquisition_ladder/streamlines_per_subject.csv") %>%
  group_by(subject, n_streamlines) %>% summarise(d = mean(d), .groups = "drop")      # mean over the subsamples
str_rt  <- read_data("acquisition_ladder/streamlines_test_retest.csv")
acq_sum <- read_data("acquisition_ladder/acquisition_summary.csv") %>% mutate(x = match(rung, names(LEVELS)))
acq_sub <- read_data("acquisition_ladder/acquisition_per_subject.csv") %>% mutate(x = match(level, names(LEVELS)))
sim     <- read_data("acquisition_ladder/acquisition_similarity.csv") %>%
  mutate(level = factor(level, names(LEVELS), gsub("\n", ", ", LEVELS)))

th   <- theme_fig() + theme(axis.text.x = element_text(size = 5.5, lineheight = .85))
YLIM <- c(-0.02, 0.9)

pA <- ggplot(mapping = aes(n_streamlines)) +
  geom_line(data = str_sub, aes(y = d, group = subject), colour = "grey70", linewidth = .25) +
  geom_line(data = str_rt, aes(y = tr_idiff), colour = BLUE, linewidth = .35, linetype = "22") +
  geom_point(data = str_rt, aes(y = tr_idiff), colour = BLUE, fill = "white", shape = 21, size = 1.3, stroke = .5) +
  geom_errorbar(data = str_sum, aes(ymin = ci_lo, ymax = ci_hi), width = 0, colour = BLUE, linewidth = .4) +
  geom_line(data = str_sum, aes(y = idiff), colour = BLUE, linewidth = .5) +
  geom_point(data = str_sum, aes(y = idiff), colour = BLUE, size = 1.5) +
  annotate("text", x = 1.15e5, y = 0.895, hjust = 0, vjust = 1, size = KEY, colour = "grey30",
           label = "filled: versus own 5 M reference\nopen: between two subsamples") +
  scale_x_log10(breaks = c(1e5, 2.5e5, 5e5, 1e6, 2e6, 5e6), labels = c("100k", "250k", "500k", "1M", "2M", "5M")) +
  coord_cartesian(ylim = YLIM) +
  labs(x = "streamlines (subsampled from 5 M)", y = "identifiability (Idiff)") + theme_fig()

pB <- ggplot(mapping = aes(x)) +
  geom_line(data = acq_sub, aes(y = d, group = subject), colour = "grey70", linewidth = .25) +
  geom_errorbar(data = acq_sum, aes(ymin = ci_lo, ymax = ci_hi), width = 0, colour = BLUE, linewidth = .4) +
  geom_line(data = acq_sum, aes(y = idiff), colour = BLUE, linewidth = .5) +
  geom_point(data = acq_sum, aes(y = idiff), colour = BLUE, size = 1.5) +
  annotate("point", x = 1, y = 0.775, shape = 25, size = 1.6, fill = "grey30", colour = "grey30") +
  annotate("text", x = 0.6, y = 0.895, hjust = 0, vjust = 1, size = KEY, colour = "grey30",
           label = "ADNI-4 protocol: 3 shells, 114 dir,\nb ≤ 2000 s mm⁻² (not measured)") +
  scale_x_continuous(breaks = 1:4, expand = expansion(add = .45),
                     labels = sprintf("%s\n(%d/%d)", LEVELS, acq_sum$n_correct_row, acq_sum$n)[order(acq_sum$x)]) +
  coord_cartesian(ylim = YLIM) +
  labs(x = "BTC acquisition, reprocessed", y = "identifiability (Idiff)") + th

own_other <- acq_sub %>%
  pivot_longer(c(self_r, other_r_mean), names_to = "reference", values_to = "r") %>%
  mutate(reference = factor(reference, c("self_r", "other_r_mean")))
own_other_mean <- own_other %>% group_by(x, reference) %>% summarise(m = mean(r), sd = sd(r), .groups = "drop")
pC <- ggplot(mapping = aes(x, colour = reference)) +
  geom_point(data = own_other, aes(y = r), position = position_jitter(width = .07, height = 0, seed = 1), size = .7, alpha = .55) +
  geom_errorbar(data = own_other_mean, aes(ymin = m - sd, ymax = m + sd), width = 0, linewidth = .4) +
  geom_line(data = own_other_mean, aes(y = m), linewidth = .5) +
  geom_point(data = own_other_mean, aes(y = m), size = 1.5) +
  annotate("text", x = 4.4, y = c(1.13, 1.055), hjust = 1, vjust = 1, size = KEY, colour = c(RED, BLUE),
           label = c("similarity to own 3-shell reference", "similarity to others’ references")) +
  scale_colour_manual(values = c(self_r = RED, other_r_mean = BLUE)) +
  scale_x_continuous(breaks = 1:4, labels = LEVELS, expand = expansion(add = .45)) +
  coord_cartesian(ylim = c(0, 1.14)) +
  labs(x = "BTC acquisition, reprocessed", y = "fingerprint similarity (r)") + th

n <- max(sim$degraded)
pD <- ggplot(sim, aes(reference, degraded, fill = r)) +
  geom_raster() + facet_wrap(~level, nrow = 1) +
  scale_fill_gradientn(colours = BLUE_RAMP(100), limits = c(0, 1), breaks = c(0, .5, 1), oob = squish, name = "r") +
  scale_x_continuous(breaks = c(1, ceiling(n / 2), n), expand = c(0, 0)) +
  scale_y_reverse(breaks = c(1, ceiling(n / 2), n), expand = c(0, 0)) +
  coord_equal() + labs(x = "reference participant (3 shells)", y = "degraded participant") +
  theme_fig() + theme(legend.position = "right", legend.key.height = unit(3.5, "mm"), legend.key.width = unit(2, "mm"),
                      legend.title = element_text(size = 6), legend.text = element_text(size = 5.5),
                      strip.text = element_text(size = 6.5), panel.grid = element_blank(), panel.spacing.x = unit(3, "mm"))

fig <- tag_panels((pA | pB | pC) / pD + plot_layout(heights = c(1, 0.78)))
save_figure(fig, "fig2_acquisition_ladder", 180, 118)
