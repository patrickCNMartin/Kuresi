
#' Plot competition scores spatially
#'
#' Creates a ggplot visualization of competition scores with optional spatial
#' image background and score binning.
#'
#' @param score Data frame with columns for x, y coordinates and score values.
#' @param cex Numeric text size multiplier (default: 10).
#' @param cex_pt Numeric point size (default: 1).
#' @param alpha Numeric transparency value between 0 and 1 (default: 0.65).
#' @param bins Integer number of bins for discrete coloring, or NULL
#'   for continuous (default: NULL).
#' @param img Optional imager object for background image.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot geom_point geom_raster
#'   scale_colour_gradientn scale_color_manual scale_fill_identity
#'   theme_classic theme element_text guides guide_legend labs
#'   scale_x_continuous scale_y_continuous aes
#' @importFrom ggnewscale new_scale
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @importFrom grDevices rgb
#' @importFrom methods is
score_plot <- function(score,
                       cex = 10,
                       cex_pt = 1,
                       alpha = 0.65,
                       bins = NULL,
                       img = NULL) {
  targets <- colnames(score)
  colnames(score) <- c("x", "y", "score")
  g <- ggplot()
  if (!is.null(img) && is(img, "imager")) {
    img <- as.data.frame(img, wide = "c") %>%
        mutate(rgb_val = rgb(c.1, c.2, c.3))
      img$x <- rev(img$x)
    g <- g + geom_raster(
      data = img,
      aes(x = x, y = y, fill = rgb_val)
    ) +
      scale_fill_identity() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      new_scale("fill")
  }
  
  if (is.null(bins)) {
    g <- g +
    geom_point(
      data = score,
      aes(x, y, col = as.numeric(score)),
      size = cex_pt,
      alpha = alpha)+
      scale_colour_gradientn(colours = c(
        "#850101",
        "#cd8878",
        "#f1f1b1",
        "#9CAAC4",
        "#1F3B70"
      ))
  } else if (!is.null(bins) && bins == 0) {
    score$score <- as.factor(score$score)
    cols <- create_palette(score)
    g <- g +
      geom_point(
        data = score,
        aes(x, y, col = score),
        size = cex_pt,
        alpha = alpha) +
      scale_color_manual(values = cols)
  } else if (!is.null(bins) && bins > 1) {
    score <- bin_score(score, bins)
    score$score <- as.factor(score$score)
    cols <- create_palette(score)
    g <- g +
      geom_point(
        data = score,
        aes(x, y, col = score),
        size = cex_pt,
        alpha = alpha) +
      scale_color_manual(values = cols)
  } else {
    stop("Excuse me? What are you giving to bin? ",
      "Because it's not a null or numeric")
  }
  g <- g +
    theme_classic() +
    theme(
      legend.text = element_text(size = cex * 1.2),
      axis.text = element_text(size = cex * 1.2),
      axis.title = element_text(size = cex * 1.2),
      plot.title = element_text(size = cex * 1.5),
      legend.title = element_text(size = cex * 1.2)
    ) +
    guides(colour = guide_legend(
      override.aes = list(size = cex * 0.3)
    )) +
    labs(
      colour = "Score", title = targets[3L],
      x = targets[1L], y = targets[2L]
    )
  return(g)

}


#' Bin score values into discrete categories
#'
#' Converts continuous score values into discrete bins and ranks them.
#'
#' @param score Data frame with a "score" column containing numeric values.
#' @param bins Integer number of bins to create.
#' @return Data frame with binned and ranked score values.
bin_score <- function(score, bins) {
  score$score <- cut(score$score, breaks = bins)
  ranked <- order(unique(score$score), decreasing = TRUE)
  locs <- match(score$score, unique(score$score)[ranked])
  score$score <- locs
  return(score)
}




#' Create color palette for scores
#'
#' Generates a color palette based on the number of unique score levels.
#' Uses a predefined color scheme with colors randomly assigned to each
#' territory. Note that as the territory plot returns a ggplot object,
#' you can easily override the color scheme.
#'
#' @param score Data frame with a "score" column (should be a factor).
#' @return Character vector of color values.
#' @importFrom grDevices colorRampPalette
create_palette <- function(score) {
  ter_col <- length(levels(score$score))
  base_colours <- rev(c(
    "#850101",
    "#cd8878",
    "#f1f1b1",
    "#9CAAC4",
    "#1F3B70"
  ))
  if (ter_col < length(base_colours)) {
    ter_pal <- colorRampPalette(base_colours[seq(1, ter_col)])
  } else {
    ter_pal <- colorRampPalette(base_colours)
  }
  ter_col <- ter_pal(ter_col)
  return(ter_col)
}


#' Create alpha values for territory highlighting
#'
#' Creates transparency (alpha) values for territories, with optional
#' highlighting of specific territories.
#'
#' @param territories Data frame with a "trial" column containing
#'   territory assignments.
#' @param highlight Numeric or character vector describing which
#'   territories should be highlighted.
#' @param alpha Numeric transparency factor (default alpha for
#'   non-highlighted territories is alpha * 0.25).
#' @details If highlight is NULL, will return the same alpha values
#'   for all territories.
#' @return Numeric vector of alpha values matching territories.
create_alpha <- function(territories, highlight, alpha) {
  if (!is.null(highlight)) {
    ter_col <- rep(alpha * 0.25,
      length(levels(territories$trial)))
    loc <- as.character(levels(territories$trial)) %in%
      highlight
    ter_col[loc] <- alpha
  } else {
    ter_col <- rep(alpha, length(levels(territories$trial)))
  }
  return(ter_col[as.integer(territories$trial)])
}


# adjust_cooridnates <- function(trial, vesalius_assay) {
#   orig_coord <- vesalius_assay@meta$orig_coord
#   #-------------------------------------------------------------------------#
#   # First let's split barcodes between adjusted and not
#   #-------------------------------------------------------------------------#
#   adj_barcodes <- grep("_et_", trial$barcodes, value = TRUE)
#   non_adj_barcodes <- trial$barcodes[!trial$barcodes %in% adj_barcodes]
#   #-------------------------------------------------------------------------#
#   # next get original coordinates and match them in trial for non adjusted
#   #-------------------------------------------------------------------------#
#   in_trial <- match(orig_coord$barcodes, non_adj_barcodes) %>%
#     na.exclude()
#   in_orig <- match(non_adj_barcodes, orig_coord$barcodes) %>%
#     na.exclude()

#   trial$x[in_trial] <- orig_coord$x_orig[in_orig]
#   trial$y[in_trial] <- orig_coord$y_orig[in_orig]
#   #-------------------------------------------------------------------------#
#   # unpack adjusted
#   #-------------------------------------------------------------------------#
#   adj_barcodes_sp <- split(adj_barcodes, "_et_")
#   for (i in seq_along(adj_barcodes_sp)) {
#     x <- orig_coord[
#       orig_coord$barcodes %in% adj_barcodes_sp[[i]],
#       "x_orig"
#     ]
#     y <- orig_coord[
#       orig_coord$barcodes %in% adj_barcodes_sp[[i]],
#       "y_orig"
#     ]
#     trial$x[trial$barcodes == adj_barcodes[i]] <- median(x)
#     trial$y[trial$barcodes == adj_barcodes[i]] <- median(y)
#   }
#   #-------------------------------------------------------------------------#
#   # adjuste using scale - note!!! This is dodgy as fuck!
#   # mainly because my original use of scale was intended to be used this way
#   # Using this like this since it is faster for ad hoc analysis
#   #-------------------------------------------------------------------------#
#   scale <- vesalius_assay@meta$scale$scale
#   trial$x <- trial$x * scale
#   trial$y <- trial$y * scale
#   return(trial)
# }


#' Visualize competition score matrix
#'
#' Creates a heatmap visualization of a competition score matrix, showing
#' scores between seed and query groups as a color-coded tile plot.
#'
#' @param scores A matrix or data frame with competition scores, where rows
#'   represent query groups and columns represent seed groups.
#' @param limits Optional numeric vector of length 2 specifying the color
#'   scale limits (default: NULL, uses data range).
#' @return A ggplot object showing the score matrix as a heatmap.
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradientn
#'   theme_minimal labs
#' @importFrom RColorBrewer brewer.pal
#' @export
view_scores <- function(scores, limits = NULL) {
  scores <- as.data.frame(scores) %>%
    # Move row names into a column called "Query"
    rownames_to_column(var = "Query") %>%
    pivot_longer(
      cols = -Query, # Reshape all other columns
      names_to = "Seed",
      values_to = "Score"
    )
  scores$Seed <- as.factor(scores$Seed)
  scores$Query <- as.factor(scores$Query)
  g <- ggplot(scores, aes(x = Seed, y = Query, fill = Score)) +
    geom_tile(color = "white", width = 1.2) +
    scale_fill_gradientn(
      colors = rev(brewer.pal(11, "Spectral")),
      limits = limits
    ) + # Adjust color scale
    theme_minimal() +
    labs(
      fill = "Score",
      x = "Seed Territories",
      y = "Query Territories")
  return(g)
}

#' Visualize ELO ratings
#'
#' Creates a bar plot showing ELO ratings for seed and query groups.
#'
#' @param elo_scores List containing "elo_seed" and "elo_query" numeric vectors.
#' @return A ggplot object showing ELO ratings as bar plots.
#' @importFrom ggplot2 ggplot geom_bar scale_fill_manual theme_minimal
#'   labs coord_cartesian theme element_text aes
#' @export
view_elo <- function(elo_scores) {
  seed <- elo_scores$elo_seed
  seed_names <- paste0("Seed_", names(seed))
  query <- elo_scores$elo_query
  query_names <- paste0("Query_", names(query))
  h <- data.frame(
    "Territory" = c(seed_names, query_names),
    "ELO" = c(seed, query),
    "Group" = c(rep("Seed", length(seed)), rep("Query", length(query)))
  )
  ord <- order(h$ELO, decreasing = TRUE)
  ord <- h$Territory[ord]
  h$Territory <- as.factor(h$Territory)
  h$Territory <- factor(h$Territory, levels = ord)
  h$ELO <- as.numeric(h$ELO)
  h$Group <- as.factor(h$Group)
  g <- ggplot(h, aes(x = Territory, y = ELO, fill = Group)) +
    geom_bar(stat = "identity") +
    coord_cartesian(ylim = c(min(h$ELO) - 10, max(h$ELO))) +
    scale_fill_manual(
      values = c("Seed" = "#284259", "Query" = "#f2d99d")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(g)
}

#' Create color palette to visualize seed and query scores
#'
#' Generates color palettes for visualizing seed and query groups,
#' assigning colors based on the union of groups present in both seed
#' and query.
#'
#' @param palette Character vector of base colors to use for palette
#'   generation.
#' @param seed Character or numeric vector of seed group identifiers.
#' @param query Character or numeric vector of query group identifiers.
#' @return A list with two character vectors: "seed" and "query"
#'   containing color assignments.
#' @importFrom grDevices colorRampPalette
#' @export
generate_palette <- function(palette, seed, query) {
  pal <- colorRampPalette(palette)
  groups <- sort(union(seed, query), decreasing = FALSE)
  n_groups <- length(unique(groups))
  if ("Not Selected" %in% groups) {
    cols <- pal(n_groups - 1)
    cols <- c(cols, "grey")
  } else {
    cols <- pal(n_groups)
  }
  seed_colors <- match(groups, unique(seed))
  seed_colors <- cols[!is.na(seed_colors)]
  query_colors <- match(groups, unique(query))
  query_colors <- cols[!is.na(query_colors)]
  return(list("seed" = seed_colors, "query" = query_colors))
}
