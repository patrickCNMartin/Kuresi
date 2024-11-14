#' @export
#' @importFrom ggplot2 ggplot geom_point aes facet_wrap theme_classic
#' @importFrom ggplot2 scale_color_manual theme element_text scale_colour_gradientn
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 guides guide_legend labs geom_polygon
#' @importFrom ggplot2 scale_fill_manual scale_x_continuous scale_y_continuous
#' @importFrom imager imrotate
#' @importFrom ggnewscale new_scale
score_plot <- function(vesalius_assay,
    trial = "last",
    cex = 10,
    cex_pt = 1,
    alpha = 0.65,
    randomise = FALSE,
    use_image = FALSE) {
    #--------------------------------------------------------------------------#
    # Dirty ggplot - this is just a quick and dirty plot to show what it look
    # like
    # At the moment - I thinking the user should make their own
    # Not a prority for custom plotting functions
    # SANITY check here and format
    #--------------------------------------------------------------------------#
    territories <- check_territory_trial(vesalius_assay, trial)
    tiles <- get_tiles(vesalius_assay)
    if (use_image) {
        img <- vesalius_assay@image
        if (length(img) == 0) {
          warning(paste0("No Image found in ",
            get_assay_names(vesalius_assay)))
          img <- NULL
        } else {
          #-------------------------------------------------------------------#
          # NOTE: not ideal - i made it a list to handle multiple images 
          # will need to think about hwo to handle this data 
          # especially when it comes to interoperable object
          #-------------------------------------------------------------------#
          img <- imager::imrotate(img[[1]], 90)
          img <- as.data.frame(img, wide = "c") %>%
            mutate(rgb_val = rgb(c.1, c.2, c.3))
          img$x <- rev(img$x)
          territories <- adjust_cooridnates(territories, vesalius_assay)
          #tiles <- adjust_cooridnates(tiles, vesalius_assay)
        }
    } else {
       img <- NULL
    }
    #--------------------------------------------------------------------------#
    # prepaing plot
    #--------------------------------------------------------------------------#
    ter_plot <- ggplot()
    if (!is.null(img)) {
      ter_plot <- ter_plot + geom_raster(data = img,
        aes(x = x, y = y, fill = rgb_val)) +
        scale_fill_identity() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        new_scale("fill")
    }

    ter_plot <- ter_plot +
        geom_point(data = territories,
          aes(x, y, col = trial),
          size = cex_pt,
          alpha = alpha)
    #--------------------------------------------------------------------------#
    # Checking if need to go spectral or cats
    #--------------------------------------------------------------------------#
    if (is(territories$trial, "numeric")) {
        ter_plot <- ter_plot +
            scale_colour_gradientn(colours = c("#850101",
            "#cd8878",
            "#f1f1b1",
            "#9CAAC4",
            "#1F3B70"))
    } else {
        cols <- create_palette(territories, randomise)
        ter_plot <- ter_plot +
            scale_color_manual(values = cols)
    }
    ter_plot <- ter_plot +
      theme_classic() +
      theme(legend.text = element_text(size = cex * 1.2),
        axis.text = element_text(size = cex * 1.2),
        axis.title = element_text(size = cex * 1.2),
        plot.title = element_text(size = cex * 1.5),
        legend.title = element_text(size = cex * 1.2)) +
      guides(colour = guide_legend(
        override.aes = list(size = cex * 0.3))) +
      labs(colour = legend, title = paste("Vesalius", trial),
        x = "X coordinates", y = "Y coordinates")
    return(ter_plot)

}


#' territories present. If required the colours will be randomly assinged
#' to each territory. Note that as the territory plot 
#' return a ggplot object, you can easily override the color scheme. 
#' @return color vector
#' @importFrom grDevices colorRampPalette
create_palette <- function(territories, randomise) {
  ter_col <- length(levels(territories$trial))
  base_colours <- rev(c("#850101",
    "#cd8878",
    "#f1f1b1",
    "#9CAAC4",
    "#1F3B70"))
  if (ter_col < length(base_colours)) {
      ter_pal <- colorRampPalette(base_colours[seq(1, ter_col)])
  } else {
      ter_pal <- colorRampPalette(base_colours)
  }
  if (randomise) {
        ter_col <- sample(ter_pal(ter_col), ter_col)
  } else {
        ter_col <- ter_pal(ter_col)
  }
  return(ter_col)
}


#' create alpha value if territories need to be highlighted
#' @param territories vesalius territories taken from a vesalius_assay
#' @param highlight numeric vector describing which territories should 
#' be highlighted
#' @param alpha tranaparent factor
#' @details If highlight is null, will return the same alpha values 
#' for all territories
#' @return vector of alpha values
create_alpha <- function(territories, highlight, alpha) {
  if (!is.null(highlight)) {
    ter_col <- rep(alpha * 0.25, length(levels(territories$trial)))
    loc <- as.character(levels(territories$trial)) %in% highlight
    ter_col[loc] <- alpha
  } else {
    ter_col <- rep(alpha, length(levels(territories$trial)))
  }
  return(ter_col[as.integer(territories$trial)])
}


adjust_cooridnates <- function(trial, vesalius_assay) {
    orig_coord <- vesalius_assay@meta$orig_coord
    #-------------------------------------------------------------------------#
    # First let's split barcodes between adjusted and not 
    #-------------------------------------------------------------------------#
    adj_barcodes <- grep("_et_", trial$barcodes, value = TRUE)
    non_adj_barcodes <- trial$barcodes[!trial$barcodes %in% adj_barcodes]
    #-------------------------------------------------------------------------#
    # next get original coordinates and match them in trial for non adjusted
    #-------------------------------------------------------------------------#
    in_trial <- match(orig_coord$barcodes, non_adj_barcodes) %>%
      na.exclude()
    in_orig <- match(non_adj_barcodes, orig_coord$barcodes) %>%
      na.exclude()
    
    trial$x[in_trial] <- orig_coord$x_orig[in_orig]
    trial$y[in_trial] <- orig_coord$y_orig[in_orig]
    #-------------------------------------------------------------------------#
    # unpack adjusted 
    #-------------------------------------------------------------------------#
    adj_barcodes_sp <- split(adj_barcodes, "_et_")
    for (i in seq_along(adj_barcodes_sp)) {
        x <- orig_coord[orig_coord$barcodes %in% adj_barcodes_sp[[i]],
          "x_orig"]
        y <- orig_coord[orig_coord$barcodes %in% adj_barcodes_sp[[i]],
          "y_orig"]
        trial$x[trial$barcodes == adj_barcodes[i]] <- median(x)
        trial$y[trial$barcodes == adj_barcodes[i]] <- median(y)
    }
    #-------------------------------------------------------------------------#
    # adjuste using scale - note!!! This is dodgy as fuck!
    # mainly because my original use of scale was intended to be used this way
    # Using this like this since it is faster for ad hoc analysis 
    #-------------------------------------------------------------------------#
    scale <- vesalius_assay@meta$scale$scale
    trial$x <- trial$x * scale
    trial$y <- trial$y * scale
    return(trial)

}


#'@importFrom tidyr pivot_longer
#'@importFrom dplyr mutate
#'@importFrom tibble rownames_to_column
#'@importFrom ggplot2 ggplot geom_tile scale_fill_gradientn theme_minimal labs
#'@importFrom RColorBrewer brewer.pal
#'@export
view_scores <- function(scores, limits = NULL){
    scores <- as.data.frame(scores) %>%
  rownames_to_column(var = "Query") %>%   # Move row names into a column called "Query"
  pivot_longer(cols = -Query,              # Reshape all other columns
               names_to = "Seed", 
               values_to = "Score")
    scores$Seed <- as.factor(scores$Seed)
    scores$Query <- as.factor(scores$Query)
    g <- ggplot(scores, aes(x = Seed, y = Query, fill = Score)) +
        geom_tile(color = "white", width = 1.2)+
        scale_fill_gradientn(colors=rev(brewer.pal(11, "Spectral")),
            limits = limits) +  # Adjust color scale
        theme_minimal() +
        labs(fill = "Score", x = "Seed Territories", y = "Query Territories")
    return(g)
}

#'@importFrom ggplot2 ggplot geom_bar scale_fill_gradientn theme_minimal labs coord_cartesian theme
#'@importFrom RColorBrewer brewer.pal
#'@export
view_elo <- function(elo_scores){
    seed <- elo_scores$elo_seed
    seed_names <- paste0("Seed_", names(seed))
    query <- elo_scores$elo_query
    query_names <- paste0("Query_", names(query))
    h <- data.frame("Territory" = c(seed_names, query_names),
        "ELO" = c(seed, query),
        "Group" = c(rep("Seed", length(seed)), rep("Query",length(query))))
    ord <- order(h$ELO, decreasing = TRUE)
    ord <- h$Territory[ord]
    h$Territory <- as.factor(h$Territory)
    h$Territory <- factor(h$Territory, levels = ord)
    h$ELO <- as.numeric(h$ELO)
    h$Group <- as.factor(h$Group)
    g <- ggplot(h, aes(x = Territory, y = ELO, fill = Group)) +
        geom_bar(stat = "identity") +
        coord_cartesian(ylim = c(min(h$ELO) - 10, max(h$ELO))) +
        scale_fill_manual(values = c("Seed" = "#284259", "Query" = "#f2d99d"))+
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    return(g)
}

#' Create color palette to vizualize seed and query scores
#' @importFrom grDevices colorRampPalette
#' @export
generate_palette <- function(palette, seed, query){
    pal <- colorRampPalette(palette)
    groups <- sort(union(seed, query), descending = FALSE)
    if ("Not Selected" %in% groups) {
        cols <- pal(length(n_groups) - 1)
        cols <- c(cols, "grey")
    } else {
        cols <- pal(n_groups)
    }
    seed_colors <- cols[match(groups, unique(seed))]
    query_colors <- cols[match(groups, unique(query))]
    return(list("seed" = seed_colors, "query" = query_colors))
}
