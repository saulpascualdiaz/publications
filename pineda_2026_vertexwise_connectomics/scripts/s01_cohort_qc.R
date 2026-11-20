# Fig. 1 — Cohort characteristics and connectome quality-control criteria.
#   a–c  age, education and sex by diagnostic group (ADNI-4, analysed sample)
#   d    cortical coverage by group, with the 40 % quality-control threshold
#   e    age and sex of the BTC controls
#   f    coverage vs. endpoint distance for every enrolled participant, with both thresholds
source("scripts/s00_setup.R")

adni <- read_data("cohort/adni_participants.csv") %>%
  mutate(group = factor(group, levels = names(GROUP_COLS)), sex = factor(sex, levels = names(SEX_COLS)))
btc  <- read_data("cohort/btc_participants.csv") %>% mutate(sex = factor(sex, levels = names(SEX_COLS)))
analysed <- adni %>% filter(analysed) %>% mutate(group = droplevels(group))

COVERAGE_MIN <- 40    # % of white-surface vertices reached by a streamline endpoint
DISTANCE_MAX <- 5     # median endpoint-to-surface distance, mm

th <- theme_fig() + theme(axis.title = element_text(size = 7), legend.text = element_text(size = 6),
                          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, -2, 0))

# "CN\n(n = 208)" axis labels, counting the participants with data for `var`
group_n <- function(var) {
  n <- analysed %>% filter(!is.na(.data[[var]])) %>% count(group)
  setNames(sprintf("%s\n(n = %d)", n$group, n$n), n$group)
}
violin <- function(var, ylab) {
  ggplot(filter(analysed, !is.na(.data[[var]])), aes(group, .data[[var]], fill = group)) +
    geom_violin(alpha = .35, colour = NA, width = .9) +
    geom_boxplot(width = .16, outlier.size = .4, linewidth = .3) +
    scale_fill_manual(values = GROUP_COLS) + scale_x_discrete(labels = group_n(var)) +
    labs(x = NULL, y = ylab) + th
}

pA <- violin("age", "age (years)")
pB <- violin("education", "education (years)")

pC <- analysed %>% filter(!is.na(sex)) %>% count(group, sex) %>%
  group_by(group) %>% mutate(p = n / sum(n)) %>%
  ggplot(aes(group, p, fill = sex)) +
  geom_col(width = .7) +
  geom_text(aes(label = n), position = position_stack(vjust = .5), size = 2.1, colour = "grey25") +
  scale_fill_manual(values = SEX_COLS, name = NULL) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, .02))) +
  scale_x_discrete(labels = group_n("sex")) +
  labs(x = NULL, y = "proportion") + th + theme(legend.position = "top")

pD <- violin("coverage_pct", "cortical coverage (%)") +
  geom_hline(yintercept = COVERAGE_MIN, linetype = "dashed", linewidth = .3, colour = "grey40") +
  annotate("text", x = 3.5, y = COVERAGE_MIN + 3.5, label = sprintf("QC gate (%d %%)", COVERAGE_MIN), size = 2, colour = "grey40", hjust = 1)

pE <- ggplot(btc, aes(sprintf("BTC controls\n(n = %d)", nrow(btc)), age)) +
  geom_boxplot(width = .35, outlier.shape = NA, linewidth = .3, fill = PALE_TEAL, alpha = .5) +
  geom_jitter(aes(colour = sex), width = .08, height = 0, size = 1.3) +
  scale_colour_manual(values = SEX_COLS, name = NULL) +
  labs(x = NULL, y = "age (years)") + th + theme(legend.position = "top")

excluded <- filter(adni, qc_fail)
pF <- ggplot(adni, aes(coverage_pct, ep_d_mm)) +
  geom_hline(yintercept = DISTANCE_MAX, linetype = "dashed", linewidth = .3, colour = "grey40") +
  geom_vline(xintercept = COVERAGE_MIN, linetype = "dashed", linewidth = .3, colour = "grey40") +
  geom_point(data = filter(adni, !qc_fail), aes(colour = group), size = .9, alpha = .7) +
  geom_point(data = excluded, shape = 4, size = 2, stroke = .6) +
  geom_text(data = excluded, aes(label = subject), size = 1.8, hjust = -0.15, vjust = -0.4) +
  scale_colour_manual(values = GROUP_COLS, name = NULL) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
  labs(x = "cortical coverage (%)", y = "endpoint distance (mm)") + th + theme(legend.position = "top")

fig <- tag_panels((pA | pB | pC) / (pD | pE | pF))
save_figure(fig, "fig1_cohort_qc", 180, 105)
