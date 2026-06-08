# ============================================================
# lumen_style.R
#
# Reusable plotting style for Lumen Computational Biology figures.
# Source this at the top of any analysis script:
#
#     source("lumen_style.R")
#
# Provides:
#   - lumen_palette        : named Lumen brand colors
#   - lumen_module_colors  : the hdWGCNA module color map
#   - theme_lumen()        : ggplot2 theme for bars / heatmaps / scatter
#   - lumen_network()      : ggraph network plot with the house style
#   - lumen_gradient_fill() / lumen_gradient_color() : accent gradients
#
# Dependencies: ggplot2 (always), igraph + ggraph (only if you
# call lumen_network).
# ============================================================

# ------------------------------------------------------------
# Palette
# ------------------------------------------------------------
lumen_palette <- c(
  bg     = "#faf7f2",   # warm paper background
  ink    = "#2a2622",   # near-black text
  sub    = "#5c544b",   # muted subtitle grey-brown
  accent = "#b15835",   # terracotta (primary accent)
  sage   = "#7d8c6e",   # sage green
  plum   = "#7c5c6b",   # muted plum
  ochre  = "#b88a3e"    # ochre
)

# hdWGCNA module colors (consistent across the whole project)
lumen_module_colors <- c(
  yellow    = "#e8b800",
  brown     = "#964b00",
  turquoise = "#4cb5b0",
  blue      = "#4472c4"
)

# A discrete categorical sequence drawn from the brand palette,
# for when you need several distinguishable colors.
lumen_discrete <- unname(lumen_palette[c("accent", "sage", "plum",
                                         "ochre", "ink")])

# ------------------------------------------------------------
# ggplot2 theme
# ------------------------------------------------------------
theme_lumen <- function(base_size = 11, grid = TRUE) {
  th <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1,
                                            color = lumen_palette["ink"]),
      plot.subtitle = ggplot2::element_text(size = base_size - 1,
                                            color = lumen_palette["sub"]),
      plot.caption  = ggplot2::element_text(size = base_size - 3,
                                            color = lumen_palette["sub"]),
      axis.text     = ggplot2::element_text(color = lumen_palette["ink"],
                                            size = base_size - 2),
      axis.title    = ggplot2::element_text(color = lumen_palette["sub"],
                                            size = base_size - 1),
      strip.text    = ggplot2::element_text(face = "bold", size = base_size,
                                            color = lumen_palette["ink"]),
      legend.title  = ggplot2::element_text(size = base_size - 2,
                                            color = lumen_palette["sub"]),
      legend.text   = ggplot2::element_text(size = base_size - 2,
                                            color = lumen_palette["ink"]),
      panel.grid.minor = ggplot2::element_blank()
    )
  if (!grid) {
    th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  th
}

# Convenience gradient scales (white -> accent)
lumen_gradient_fill <- function(name = NULL, ...) {
  ggplot2::scale_fill_gradient(
    low = lumen_palette["bg"], high = lumen_palette["accent"],
    name = name, ...
  )
}
lumen_gradient_color <- function(name = NULL, ...) {
  ggplot2::scale_color_gradient(
    low = lumen_palette["bg"], high = lumen_palette["accent"],
    name = name, ...
  )
}

# ------------------------------------------------------------
# ggraph network plot in the house style
# ------------------------------------------------------------
# Draws a directed regulatory network with:
#   - TF nodes as squares, target nodes as circles
#   - edges colored by a single accent color (or per-module)
#   - edge width mapped to a weight column
#   - repelled labels, TFs bold
#
# Arguments:
#   edges      data.frame with columns: from, to, and a weight col
#   weight_col name of the numeric weight column in `edges`
#   tf_nodes   character vector of node names that are TFs (squares)
#   accent     edge/node color (default Lumen accent; pass a module
#              color like lumen_module_colors["yellow"] for themed nets)
#   title      optional plot title
#   layout     ggraph layout (default "fr")
#   max_targets_per_tf  optional cap on targets per TF for readability
#
# Returns a ggraph/ggplot object.
# ------------------------------------------------------------
lumen_network <- function(edges,
                          weight_col = "weight",
                          tf_nodes = NULL,
                          accent = unname(lumen_palette["accent"]),
                          title = NULL,
                          layout = "fr",
                          max_targets_per_tf = NULL,
                          seed = 42) {
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    stop("lumen_network requires the 'igraph' and 'ggraph' packages.")
  }
  set.seed(seed)

  e <- edges
  names(e)[names(e) == weight_col] <- "weight"

  # Optional thinning for readability
  if (!is.null(max_targets_per_tf)) {
    e <- e |>
      dplyr::group_by(from) |>
      dplyr::slice_max(weight, n = max_targets_per_tf, with_ties = FALSE) |>
      dplyr::ungroup()
  }

  g <- igraph::graph_from_data_frame(
    e[, c("from", "to", "weight")], directed = TRUE
  )

  if (is.null(tf_nodes)) tf_nodes <- unique(e$from)
  igraph::V(g)$is_tf <- igraph::V(g)$name %in% tf_nodes

  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight),
      color = accent, alpha = 0.4,
      arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"),
                             type = "closed"),
      end_cap = ggraph::circle(3, "mm")
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = ifelse(is_tf, 6, 3),
                   shape = ifelse(is_tf, "TF", "target")),
      color = accent
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = name,
                   fontface = ifelse(is_tf, "bold", "plain")),
      size = 3, repel = TRUE, color = unname(lumen_palette["ink"])
    ) +
    ggraph::scale_edge_width(range = c(0.3, 2), guide = "none") +
    ggplot2::scale_size_identity() +
    ggplot2::scale_shape_manual(values = c("TF" = 15, "target" = 16),
                                guide = "none") +
    ggplot2::theme_void()

  if (!is.null(title)) {
    p <- p + ggplot2::labs(title = title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12,
                                           color = accent, hjust = 0.5)
      )
  }
  p
}

message("Lumen style loaded: theme_lumen(), lumen_network(), ",
        "lumen_palette, lumen_module_colors")
