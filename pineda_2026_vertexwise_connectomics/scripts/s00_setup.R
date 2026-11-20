# s00_setup.R — shared packages, paths, palette, theme and helpers.
# Sourced at the top of every figure script; run the scripts from the repository root.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(scales)
})

# ---- paths (override with the environment variables DATA_DIR and FIG_DIR) -------------------
DATA_DIR <- Sys.getenv("DATA_DIR", "data")
FIG_DIR  <- Sys.getenv("FIG_DIR", "figures")

# ---- palette -------------------------------------------------------------------------------
GROUP_COLS  <- c(CN = "#009E73", MCI = "#E69F00", AD = "#D55E00", `non-AD dem` = "#999999")   # Okabe-Ito
SEX_COLS    <- c(Female = "#E8B0D0", Male = "#9BBBE0")
BAND_LEVELS <- c("short", "medium", "long")
BAND_COLS   <- c(short = "#7FA5D6", medium = "#AEC4E8", long = "#D8E3F3")
BLUE        <- "#3E6CA6"     # main line colour; between-participant similarity
RED         <- "#C0504D"     # within-participant similarity
PALE_BLUE   <- "#AEC4E8"     # pale fill
PALE_TEAL   <- "#9FD8CB"     # BTC fill
BLUE_RAMP   <- colorRampPalette(c("#EFF4FB", "#AEC4E8", "#5E86C0"))
MODEL_COLS  <- c(m1_pooled = BLUE, m2_pooled_dx = "grey10", m3_patients = "#D55E00", m4_cn = "#009E73")

# ---- theme ---------------------------------------------------------------------------------
KEY   <- 2      # in-panel key text (mm)
SMALL <- 1.8    # smallest in-panel text (mm)

theme_fig <- function() {
  theme_minimal(base_size = 7) +
    theme(axis.text = element_text(size = 6), panel.grid.minor = element_blank(),
          legend.position = "none", plot.margin = margin(2, 4, 2, 2))
}

# add bold lower-case panel letters; `tags` is "a" (automatic) or a vector with "" for untagged panels
tag_panels <- function(fig, tags = "a") {
  fig + plot_annotation(tag_levels = if (length(tags) == 1) tags else list(tags)) &
    theme(plot.tag = element_text(face = "bold", size = 8))
}

# ---- input / output ------------------------------------------------------------------------
read_data <- function(file, ...) {
  path <- file.path(DATA_DIR, file)
  if (!file.exists(path)) stop("Missing ", path, " (see the Data section of README.md)", call. = FALSE)
  read_csv(path, show_col_types = FALSE, ...)
}

# name,value table -> named numeric vector
read_values <- function(file) { v <- read_data(file); setNames(v$value, v$name) }

# PNG (300 dpi) + PDF in figures/
save_figure <- function(plot, name, width_mm, height_mm) {
  dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = w, height = h, dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = w, height = h, bg = "white", device = cairo_pdf)
  message("saved ", file.path(FIG_DIR, name), ".png / .pdf")
}

# ---- number formatting ---------------------------------------------------------------------
fmt   <- function(x, digits = 1) formatC(x, format = "f", digits = digits)
fmt_r <- function(x) paste0(ifelse(x < 0, "−", ""), sub("^0", "", fmt(abs(x), 2)))                # .28, −.05
fmt_p <- function(p) ifelse(p < 0.01, sub("0+$", "", sub("^0", "", fmt(p, 4))), sub("^0", "", fmt(p, 3)))
