# s00_surface.R — cluster-corrected vertex maps on the inflated fsaverage surface (Figs. 6 and 7).
# Vertex maps are stored on fsaverage5; each fsaverage vertex takes the value of its nearest
# fsaverage5 vertex on the sphere. Views are rendered off screen with fsbrain / rgl.
#
# Needs FreeSurfer's fsaverage: set SUBJECTS_DIR (or FREESURFER_HOME), or run once
#   fsbrain::download_fsaverage(accept_freesurfer_license = TRUE)

suppressPackageStartupMessages({ library(fsbrain); library(rgl); library(magick); library(freesurferformats) })
options(rgl.useNULL = TRUE)

SUBJECTS_DIR <- fsbrain::fsaverage.path()
N_FSA5    <- 10242L        # vertices per hemisphere, fsaverage5 (= the first 10,242 vertices of fsaverage)
VIEWS     <- c(sd_lateral_lh = "left, lateral", sd_medial_lh = "left, medial", sd_medial_rh = "right, medial", sd_lateral_rh = "right, lateral")
GREY_NS   <- "#DADADA"     # cortex outside the significant clusters
GREY_WALL <- "#B5B5B5"     # medial wall
RAMP_POS  <- colorRampPalette(c("#E8A6A1", "#C0504D", "#8B2E2B"))(100)
RAMP_NEG  <- colorRampPalette(c("#AEC4E8", "#5E86C0", "#2F4F8A"))(100)
ALPHA     <- 0.05          # cluster-level family-wise error threshold

# index of the nearest fsaverage5 vertex for every fsaverage vertex of one hemisphere
nearest_fsa5 <- function(hemi, chunk = 2048L) {
  s <- read.fs.surface(file.path(SUBJECTS_DIR, "fsaverage", "surf", paste0(hemi, ".sphere")))$vertices
  s <- s / sqrt(rowSums(s^2)); s5 <- t(s[seq_len(N_FSA5), ])
  unlist(lapply(split(seq_len(nrow(s)), ceiling(seq_len(nrow(s)) / chunk)),
                function(i) max.col(s[i, , drop = FALSE] %*% s5, ties.method = "first")), use.names = FALSE)
}
FSA5_INDEX <- list(lh = nearest_fsa5("lh"), rh = nearest_fsa5("rh"))

# vertices displayed: cortex inside a cluster significant at pFWE < ALPHA
shown <- function(map) map$cortex == 1 & !is.na(map$cluster_pfwe) & map$cluster_pfwe < ALPHA

# colour range shared by a set of maps: largest |t| displayed (rounded up), and whether any displayed t is negative
t_range <- function(maps, t_min) {
  t <- unlist(lapply(maps, function(m) m$t[shown(m)]))
  list(max = ceiling(10 * max(abs(t), t_min + 0.1)) / 10, two_sided = any(t < 0))
}

# one colour per fsaverage vertex of a hemisphere: |t| from t_min to t_max on the red (t > 0) or blue (t < 0) ramp
vertex_colours <- function(map, hemi, t_min, t_max) {
  m <- map[map$hemi == hemi, ]; m <- m[order(m$hemi_vertex), ]
  col <- ifelse(m$cortex == 1, GREY_NS, GREY_WALL); s <- shown(m)
  k <- 1L + as.integer(round(99 * pmin(pmax((abs(m$t[s]) - t_min) / (t_max - t_min), 0), 1)))
  col[s] <- ifelse(m$t[s] >= 0, RAMP_POS[k], RAMP_NEG[k])
  col[FSA5_INDEX[[hemi]]]
}

# the four views of one map as trimmed images (the snapshot is retried when it comes back blank)
render_map <- function(map, t_min, t_max, px = 1000L) {
  lh <- vertex_colours(map, "lh", t_min, t_max); rh <- vertex_colours(map, "rh", t_min, t_max)
  lapply(setNames(names(VIEWS), names(VIEWS)), function(view) {
    vis.color.on.subject(SUBJECTS_DIR, "fsaverage", lh, rh, surface = "inflated", views = view, rgloptions = list(windowRect = c(0L, 0L, px, px)))
    on.exit(try(close3d(), silent = TRUE))
    file <- tempfile(fileext = ".png")
    for (attempt in 1:6) {
      try(snapshot3d(file, webshot = TRUE), silent = TRUE)
      img <- tryCatch(image_trim(image_read(file)), error = function(e) NULL)
      if (!is.null(img) && all(unlist(image_info(img)[c("width", "height")]) > 80)) return(img)
    }
    stop("blank snapshot for view ", view)
  })
}

# largest tile over several maps, so that the brain scale is identical across panels
tile_size <- function(tiles) {
  info <- do.call(rbind, lapply(unlist(tiles, recursive = FALSE), image_info))
  c(w = max(info$width), h = max(info$height))
}

# ggplot panel with one row of views per map: view names on top, a text label above each row.
# The height-to-width ratio of the panel is returned as attribute "aspect".
brain_panel <- function(tiles, labels, size, width_px = 2400L) {
  tw <- as.integer(size[["w"]]); th <- as.integer(size[["h"]])
  gap <- as.integer(round(0.03 * tw)); lab <- as.integer(round(0.19 * th)); hdr <- as.integer(round(0.10 * th))
  rows <- lapply(tiles, function(views) {
    views <- lapply(views, image_extent, geometry = sprintf("%dx%d", tw, th), gravity = "center", color = "white")
    row <- views[[1]]
    for (v in views[-1]) row <- image_append(c(row, image_blank(gap, th, "white"), v))
    image_append(c(image_blank(image_info(row)$width, lab, "white"), row), stack = TRUE)
  })
  block <- image_blank(image_info(rows[[1]])$width, hdr, "white")
  for (r in rows) block <- image_append(c(block, r), stack = TRUE)
  k <- width_px / image_info(block)$width                                    # rescale to a fixed width
  block <- image_resize(block, sprintf("%dx%d!", width_px, as.integer(round(image_info(block)$height * k))))
  W <- image_info(block)$width; H <- image_info(block)$height
  x_view <- (seq_along(VIEWS) - 1) * (tw + gap) * k + tw * k / 2
  y_row  <- H - hdr * k - (seq_along(rows) - 1) * (lab + th) * k - 0.5 * lab * k
  p <- ggplot() + annotation_raster(as.raster(block), 0, W, 0, H) +
    annotate("text", x = x_view, y = H - 0.5 * hdr * k, label = VIEWS, size = SMALL, colour = "grey45") +
    annotate("text", x = 0, y = y_row, label = labels, hjust = 0, size = KEY, colour = "grey20", lineheight = .9) +
    coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(1, 2, 0, 2))
  structure(p, aspect = H / W)
}

# thresholded colour bar (blue | grey | red when two_sided, grey | red otherwise)
colour_bar <- function(t_min, t_max, two_sided, title) {
  seg <- function(lo, hi, cols) { e <- seq(lo, hi, length.out = length(cols) + 1); data.frame(xmin = head(e, -1), xmax = e[-1], fill = cols) }
  if (two_sided) {
    d <- rbind(seg(-t_max, -t_min, rev(RAMP_NEG)), data.frame(xmin = -t_min, xmax = t_min, fill = GREY_NS), seg(t_min, t_max, RAMP_POS))
    breaks <- c(-t_max, -t_min, t_min, t_max)
  } else {
    d <- rbind(data.frame(xmin = t_min - 0.6 * (t_max - t_min), xmax = t_min, fill = GREY_NS), seg(t_min, t_max, RAMP_POS))
    breaks <- c(t_min, t_max)
  }
  ggplot(d) + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = fill)) + scale_fill_identity() +
    scale_x_continuous(breaks = breaks, labels = sub("-", "−", fmt(breaks)), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(d$xmin[1], t_max), clip = "off") + labs(x = title, y = NULL) + theme_fig() +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_line(colour = "grey40", linewidth = .3),
          axis.ticks.length.x = unit(1.2, "pt"), axis.title.x = element_text(size = 6, margin = margin(t = 1)), plot.margin = margin(0, 2, 1, 2))
}
