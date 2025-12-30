# #' Loading Visium data anc creating a Vesalius object
# #' @param input_dir character string - path to directory containing visium data
# #' @param tag character string - file tag used to define data set name
# #' @return vesalius_asssay object
# #' @importFrom Seurat Read10X_h5
# #' @importFrom Matrix readMM
# #' @importFrom rjson fromJSON
# #' @importFrom imager load.image
# #' @importFrom vesalius build_vesalius_assay
# #' @export
# load_visium <- function(input_dir, tag = NULL) {
#   #-------------------------------------------------------------------------#
#   # check if spatial folder exists
#   #-------------------------------------------------------------------------#
#   sub_folder <- list.dirs(path = input_dir, recursive = FALSE)
#   sub_folder <- ifelse(length(sub_folder) == 0, "", sub_folder)
#   if (grepl(pattern = "spatial", x = sub_folder)) {
#     coord <- list.files(
#       path = paste0(input_dir, "/spatial"),
#       pattern = "tissue_positions",
#       full.names = TRUE
#     )
#     if (length(coord) == 0) {
#       stop(paste0(
#         "No Spatial coordinates found in input directory:",
#         input_dir, "/spatial"
#       ))
#     }
#     image <- list.files(
#       path = paste0(input_dir, "/spatial"),
#       pattern = "hires_image|lowres_image",
#       full.names = TRUE
#     )
#     if (length(image) == 0) {
#       stop(paste0(
#         "No hires_image found in input directory: ",
#         input_dir, "/spatial"
#       ))
#     }
#     js <- list.files(
#       path = paste0(input_dir, "/spatial"),
#       pattern = ".json",
#       full.names = TRUE
#     )
#     if (length(js) == 0) {
#       stop(paste0(
#         "No scaling data (.json) found in input directory:",
#         input_dir, "/spatial"
#       ))
#     }
#   } else {
#     #-------------------------------------------------------------------------#
#     # assume that all files in the current directory
#     #-------------------------------------------------------------------------#
#     coord <- list.files(
#       path = input_dir,
#       pattern = "tissue_positions",
#       full.names = TRUE
#     )
#     if (length(coord) == 0) {
#       stop(paste0(
#         "No Spatial coordinates found in input directory:",
#         input_dir
#       ))
#     }
#     image <- list.files(
#       path = input_dir,
#       pattern = "hires_image|lowres_image",
#       full.names = TRUE
#     )
#     if (length(image) == 0) {
#       stop(paste0(
#         "No hires_image found in input directory:",
#         input_dir
#       ))
#     }
#     js <- list.files(
#       path = input_dir,
#       pattern = ".json",
#       full.names = TRUE
#     )
#     if (length(js) == 0) {
#       stop(paste0(
#         "No scaling data (.json) found in input directory:",
#         input_dir
#       ))
#     }
#   }
#   #-------------------------------------------------------------------------#
#   # counts should be in current directory and we also create a plot tag
#   #-------------------------------------------------------------------------#
#   counts <- list.files(input_dir,
#     pattern = "filtered_feature_bc_matrix.h5|.mtx",
#     full.names = TRUE
#   )
#   #-------------------------------------------------------------------------#
#   # filter based on tags
#   #-------------------------------------------------------------------------#
#   coord <- grep(tag, coord, value = TRUE)
#   counts <- grep(tag, counts, value = TRUE)
#   image <- grep(tag, image, value = TRUE)
#   js <- grep(tag, js, value = TRUE)
#   #-------------------------------------------------------------------------#
#   # Load data
#   #-------------------------------------------------------------------------#
#   coord <- read.csv(coord, header = FALSE)
#   coord <- coord[coord[, 2] == 1, ]
#   coord <- coord[, c(1, 5, 6)]
#   colnames(coord) <- c("barcodes", "x", "y")
#   coord$x <- as.numeric(coord$x)
#   coord$y <- as.numeric(coord$y)

#   if (grepl(pattern = ".h5", x = counts)) {
#     counts <- Seurat::Read10X_h5(counts)
#   } else {
#     features <- list.files(input_dir,
#       pattern = "features",
#       full.names = TRUE
#     )
#     features <- grep(tag, features, value = TRUE)
#     barcodes <- list.files(input_dir,
#       pattern = "barcodes",
#       full.names = TRUE
#     )
#     barcodes <- grep(tag, barcodes, value = TRUE)
#     feat <- read.table(features)[, 2]
#     bar <- read.table(barcodes)[, 1]
#     counts <- Matrix::readMM(counts)
#     counts <- as(counts, "CsparseMatrix")
#     rownames(counts) <- feat
#     counts <- counts[!duplicated(feat), ]
#     colnames(counts) <- bar
#   }
#   if (length(image) > 1) {
#     img <- imager::load.image(grep(x = image, pattern = "hires", value = TRUE))
#     scale <- rjson::fromJSON(file = js)$tissue_hires_scalef
#   } else {
#     img <- imager::load.image(image)
#     scale <- rjson::fromJSON(file = js)$tissue_lowres_scalef
#   }

#   #-------------------------------------------------------------------------#
#   # build assay
#   #-------------------------------------------------------------------------#
#   ves <- build_vesalius_assay(coord,
#     counts,
#     image = img,
#     scale = scale,
#     assay = tag,
#     verbose = FALSE
#   )
#   return(ves)
# }

#' Create Win/Lose gene list
#'
#' Returns predefined lists of "win" and "lose" genes associated with
#' competitive fitness in cancer cells.
#'
#' @return A list containing two character vectors: "win"
#'   (pro-competitive genes) and "lose" (anti-competitive genes).
#' @export
win_lose_genes <- function() {
  win <- toupper(c(
    "Myc", "Mycn", "Cacfd1", "Col17a1", "Yap1", "Taz",
    paste0("Tead", 1:4), "Scrib", paste0("Stat", 1:6), paste0("Jak", 1:3),
    "Kras", "Hras", "Egfr", "Src",
    paste0("Bmp", 1:15), "Mtor", "Ppm1d"
  ))
  lose <- toupper(c(
    "Sparc", paste0("Mapk", 8:10), "Mapk14", "Apc", "Tp53"))
  return(list("win" = win, "lose" = lose))
}

#' Create Cancer marker list
#'
#' Returns predefined lists of cancer marker genes for specific cancer types.
#'
#' @param type Character string specifying cancer type: "pancreas",
#'   "breast", or "ovary".
#' @return Character vector of cancer marker gene names.
#' @export
cancer_maker_list <- function(type) {
  genes <- switch(type,
    "pancreas" = c(
      "BRCA1", "BRCA2", "PALB2", "CDKN2A", "ATM",
      "TP53", "STK11", "MLH1", "MSH2",
      "MSH6", "PMS2", "EPCAM"
    ),
    "breast" = c(
      "BRCA1", "BRCA2", "CHEK2", "PALB2", "TM",
      "ATM", "BARD1", "BRIP1", "CHEK1",
      "CDH1", "NF1", "PTEN", "RAD51C", "RAD51D",
      "STK11", "TP53"
    ),
    "ovary" = c(
      "ATM", "BRCA1", "BRCA2", "BRIP1", "EPCAM",
      "MLH1", "MSH2", "MSH6", "NBN", "PALB2", "RAD51C",
      "RAD51D"
    )
  )
  return(genes)
}




#' Generate function based on score type
#'
#' Returns an appropriate aggregation function (mean, sum, or Matrix
#' functions) based on score type and context parameters.
#'
#' @param score_type Character string containing "mean" or "sum" to
#'   specify function type.
#' @param single Logical whether working with single values
#'   (default: TRUE).
#' @param by_territory Logical whether aggregating by territory
#'   (default: TRUE).
#' @return A function (mean, sum, Matrix::colMeans, or Matrix::colSums).
#' @importFrom Matrix rowMeans rowSums colMeans colSums
function_call <- function(
    score_type,
    single = TRUE,
    by_territory = TRUE) {
  if (single && grepl("mean", score_type)) {
    fun <- mean
  } else if (single && grepl("sum", score_type)) {
    fun <- sum
  } else if (!single && grepl("mean", score_type) && by_territory) {
    fun <- Matrix::colMeans
  } else if (!single && grepl("sum", score_type) && by_territory) {
    fun <- Matrix::colSums
  } else if (!single && grepl("mean", score_type) &&
    !by_territory) {
    fun <- Matrix::colMeans
  } else if (!single && grepl("sum", score_type) &&
    !by_territory) {
    fun <- Matrix::colSums
  } else {
    stop("Cannot generate function from score type and by territory",
      " tags")
  }
  return(fun)
}

# #' @importFrom vesalius get_tiles get_counts build_vesalius_assay
# #' @importFrom dplyr %>% filter
# #' @export
# filter_assay <- function(
#     vesalius_assay,
#     cells = NULL,
#     territories = NULL,
#     trial = "last",
#     genes = NULL) {
#   if (!is.null(territories)) {
#     ter_keep <- check_territory_trial(vesalius_assay, trial)
#     ter_keep <- ter_keep$barcodes[ter_keep$trial %in% territories]
#   } else {
#     ter_keep <- c()
#   }
#   if (!is.null(cells)) {
#     cell_keep <- get_tiles(vesalius_assay) %>% filter(origin == 1)
#     cell_keep <- cell_keep$barcodes[cell_keep$barcodes %in% cells]
#   } else {
#     cell_keep <- c()
#   }
#   if (!is.null(genes)) {
#     gene_keep <- rownames(get_counts(vesalius_assay))
#     gene_keep <- gene_keep[gene_keep %in% genes]
#   } else {
#     gene_keep <- c()
#   }
#   ### NOTE this is buggy!!! Will not hold to edge cases
#   ### Need to fix for long term use!!!!
#   keep_barcodes <- unique(c(ter_keep, cell_keep))
#   if (length(keep_barcodes) == 0) {
#     stop("No Barcodes present under current subset conditions!")
#   } else {
#     coord <- vesalius_assay@meta$orig_coord
#     coord <- coord[
#       coord$barcodes %in% keep_barcodes,
#       c("barcodes", "x_orig", "y_orig")
#     ]
#     colnames(coord) <- c("barcodes", "x", "y")
#     counts <- get_counts(vesalius_assay, "raw")
#     counts <- counts[, colnames(counts) %in% keep_barcodes]
#   }
#   if (length(gene_keep) == 0 && !is.null(genes)) {
#     stop("Requested genes are not in count matrix")
#   } else if (length(gene_keep) > 0 && !is.null(genes)) {
#     counts <- counts[rownames(counts) %in% gene_keep, ]
#   }
#   if (length(vesalius_assay@image) > 0) {
#     image <- vesalius_assay@image[[1]]
#     scale <- vesalius_assay@meta$scale$scale
#   } else {
#     image <- NULL
#     scale <- "auto"
#   }
#   vesalius_assay <- build_vesalius_assay(coord,
#     counts,
#     image = image,
#     scale = scale,
#     verbose = FALSE
#   )
#   return(vesalius_assay)
# }

# #' @importFrom Seurat AddModuleScore
# #' @importFrom vesalius get_territories
# #' @export
# module_score <- function(
#     vesalius_assay,
#     seurat_object = NULL,
#     gene_list,
#     threshold = 0,
#     filter = FALSE,
#     add_name = NULL) {
#   if (!is.null(seurat_object)) {
#     seurat_object <- AddModuleScore(seurat_object,
#       list(gene_list),
#       ctrl = length(gene_list),
#       name = "module",
#       scale = TRUE
#     )

#     score <- seurat_object@meta.data
#   } else {
#     score <- count_score(vesalius_assay, gene_list)
#   }


#   territories <- get_territories(vesalius_assay)
#   if (!is.null(add_name)) {
#     new_trial <- create_trial_tag(
#       colnames(territories),
#       add_name
#     ) %>%
#       tail(1)
#   } else {
#     new_trial <- create_trial_tag(
#       colnames(territories),
#       "Module_score"
#     ) %>%
#       tail(1)
#   }
#   score <- score[match(territories$barcodes, rownames(score)), ]
#   territories <- cbind(territories[, c("barcodes", "x", "y")], score$module1)
#   colnames(territories) <- c("barcodes", "x", "y", new_trial)
#   vesalius_assay <- update_vesalius_assay(
#     vesalius_assay,
#     territories,
#     "territories"
#   )
#   commit <- create_commit_log(
#     arg_match = as.list(match.call()),
#     default = formals(module_score)
#   )
#   vesalius_assay <- commit_log(
#     vesalius_assay,
#     commit,
#     vesalius_assay@assay
#   )

#   if (filter) {
#     keep <- rownames(score)[score$module1 > threshold]
#     vesalius_assay <- filter_assay(vesalius_assay, cells = keep)
#   }
#   return(vesalius_assay)
# }


# count_score <- function(vesalius_assay, gene_list) {
#   count <- get_counts(vesalius_assay, type = "raw")
#   count <- count[rownames(count) %in% gene_list, ]
#   if (is.null(count)) {
#     score <- count
#   } else {
#     score <- colSums(count)
#   }
#   score <- data.frame("barcodes" = colnames(count), "module1" = score)
#   rownames(score) <- colnames(count)
#   return(score)
# }

# #' @export
# compress_score <- function(
#     vesalius_assay,
#     trial = "last",
#     add_name = NULL,
#     rank = TRUE,
#     bins = 1) {
#   trial <- check_territory_trial(vesalius_assay, trial)
#   if (rank) {
#     if (bins > 1) {
#       trial$trial <- cut(trial$trial, breaks = bins)
#     }
#     ranked <- order(unique(trial$trial), decreasing = TRUE)
#     locs <- match(trial$trial, unique(trial$trial)[ranked])
#     trial$trial <- locs
#   } else {
#     if (bins > 1) {
#       trial$trial <- cut(trial$trial, breaks = bins)
#     } else {
#       warning("No valid parameters - returning assay as is")
#       return(vesalius_assay)
#     }
#   }
#   if (!is.null(add_name)) {
#     new_trial <- create_trial_tag(
#       colnames(trial),
#       add_name
#     ) %>%
#       tail(1)
#   } else if (is.null(add_name) && rank) {
#     new_trial <- create_trial_tag(
#       colnames(trial),
#       "ranked"
#     ) %>%
#       tail(1)
#   } else {
#     new_trial <- create_trial_tag(
#       colnames(trial),
#       "binned"
#     ) %>%
#       tail(1)
#   }
#   colnames(trial) <- c("barcodes", "x", "y", new_trial)
#   vesalius_assay <- update_vesalius_assay(
#     vesalius_assay,
#     trial,
#     "territories"
#   )
#   commit <- create_commit_log(
#     arg_match = as.list(match.call()),
#     default = formals(compress_score)
#   )
#   vesalius_assay <- commit_log(
#     vesalius_assay,
#     commit,
#     vesalius_assay@assay
#   )
#   return(vesalius_assay)
# }

# #' @export
# get_max_range <- function(vesalius_assay) {
#   coord <- get_tiles(vesalius_assay)
#   max_range <- max(c(max(coord$x), max(coord$y)))
#   return(max_range)
# }


# filter_territories <- function(territories, trial, min_spatial_index) {
#   tmp <- unlist(sapply(territories, function(ter, trial, min_spatial_index) {
#     buffer <- nrow(trial[trial$trial == ter, ]) < min_spatial_index | ter == "isolated"
#     return(buffer)
#   }, trial, min_spatial_index, simplify = FALSE))
#   return(territories[!tmp])
# }





# #' Integrate counts using Seurat
# #' @param matched matrix - matrix containing counts from matched/mapped assay
# #' @param reference matrix - matrix containing counts from reference assay
# #' @param method character - Seurat integration method to use
# #' @param nfeatures integer - number of features to use during integration
# #' @param dimensions interger - number of dimensions integrated latent space
# #' dimensions.
# #' @param infer logical - back infer original counts by reversing reduced
# #' dimensional space roations.
# #' @param signal character - defining which signal should be returned:
# #' variable_features, all_features or custom gene list.
# #' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
# #' @importFrom Seurat ScaleData RunPCA IntegrateLayers
# count_integration <- function(
#     query,
#     seed,
#     integration_method = "cca",
#     nfeatures = 2000,
#     dimensions = 30,
#     features = NULL) {
#   query <- Seurat::CreateSeuratObject(query)
#   seed <- Seurat::CreateSeuratObject(seed)

#   if (integration_method != "none") {
#     integrated <- list(seed, query)
#     integrated <- lapply(integrated, function(x, nfeatures) {
#       x <- Seurat::NormalizeData(x, verbose = FALSE) %>%
#         Seurat::FindVariableFeatures(verbose = FALSE, nfeatures = nfeatures)
#       return(x)
#     }, nfeatures = nfeatures)
#     features <- pad_features(features, integrated, nfeatures)
#     anchors <- suppressWarnings(Seurat::FindIntegrationAnchors(
#       object.list = integrated,
#       anchor.features = features,
#       dims = seq(1, dimensions),
#       reduction = integration_method,
#       verbose = FALSE
#     ))
#     integrated <- suppressWarnings(Seurat::IntegrateData(anchors,
#       dims = seq(1, dimensions),
#       verbose = FALSE
#     ))
#     counts <- integrated@assays$integrated
#     counts <- list("raw" = counts@counts, "data" = counts@data)
#   } else {
#     integrated <- merge(seed, query)
#     integrated <- SeuratObject::JoinLayers(integrated)
#     integrated <- Seurat::NormalizeData(integrated, verbose = FALSE) %>%
#       Seurat::FindVariableFeatures(verbose = FALSE, nfeatures = nfeatures)
#     counts <- list(
#       "raw" = Seurat::GetAssayData(integrated, layer = "counts")[features, ],
#       "data" = Seurat::GetAssayData(integrated, layer = "data")[features, ]
#     )
#   }

#   return(counts)
# }

# pad_features <- function(features, integrated, nfeatures) {
#   if (length(features) >= nfeatures) {
#     features <- features[sample(features, size = nfeatures)]
#   } else if (length(integrated) == 2) {
#     pad_seed <- Seurat::VariableFeatures(integrated[[1]])
#     pad_query <- Seurat::VariableFeatures(integrated[[2]])
#     pad <- intersect(pad_seed, pad_query)
#     pad <- pad[!pad %in% features]

#     if (length(pad) < nfeatures) {
#       pad_length <- nfeatures - length(pad)
#     }
#     if (pad_length > 0) {
#       features <- c(features, pad[seq(1, pad_length)])
#     }
#   } else {
#     pad <- Seurat::VariableFeatures(integrated)
#     pad <- pad[!pad %in% features]
#     if (length(pad) < nfeatures) {
#       pad_length <- nfeatures - length(pad)
#     }
#     if (pad_length > 0) {
#       features <- c(features, pad[seq(1, pad_length)])
#     }
#   }
#   return(features)
# }

# #' @importFrom vesalius get_coordinates
# quick_assay_build <- function(
#     seed_assay,
#     query_assay,
#     integrated,
#     seed_trial,
#     query_trial,
#     seed,
#     query,
#     min_spatial_index = 10) {
#   #-------------------------------------------------------------------------#
#   # first filter out unwanted territories in seed
#   #-------------------------------------------------------------------------#
#   seed_coord <- vesalius::get_coordinates(seed_assay)[
#     ,
#     c("barcodes", "x", "y")
#   ]
#   seed_ter <- vesalius::get_territories(seed_assay)[
#     ,
#     c("barcodes", "x", "y", seed_trial)
#   ]
#   buffer <- table(seed_ter[, seed_trial])
#   seed <- names(buffer)[buffer >= min_spatial_index &
#     names(buffer) != "isolated"]
#   seed_locs <- seed_ter$barcodes[seed_ter[, seed_trial] %in% seed]
#   seed_coord <- seed_coord[seed_coord$barcodes %in% seed_locs, ]
#   seed_ter <- seed_ter[seed_ter$barcodes %in% seed_locs, ]

#   #-------------------------------------------------------------------------#
#   # next in query
#   #-------------------------------------------------------------------------#
#   query_coord <- vesalius::get_coordinates(query_assay)[
#     ,
#     c("barcodes", "x", "y")
#   ]
#   query_ter <- vesalius::get_territories(query_assay)[
#     ,
#     c("barcodes", "x", "y", query_trial)
#   ]
#   buffer <- table(query_ter[, query_trial])
#   query <- names(buffer)[buffer >= min_spatial_index &
#     names(buffer) != "isolated"]
#   query_locs <- query_ter$barcodes[query_ter[, query_trial] %in% query]
#   query_coord <- query_coord[query_coord$barcodes %in% query_locs, ]
#   query_ter <- query_ter[query_ter$barcodes %in% query_locs, ]

#   #-------------------------------------------------------------------------#
#   # then build full cooridnates
#   #-------------------------------------------------------------------------#
#   coordinates <- rbind(seed_coord, query_coord)
#   territories <- data.frame(
#     "barcodes" = c(seed_ter$barcodes, query_ter$barcodes),
#     "x" = c(seed_ter$x, query_ter$x),
#     "y" = c(seed_ter$y, query_ter$y),
#     "sample" = c(
#       rep("reference", nrow(seed_ter)),
#       rep("matched", nrow(query_ter))
#     ),
#     "trial" = c(
#       seed_ter[, seed_trial],
#       query_ter[, query_trial]
#     )
#   )
#   colnames(territories) <- gsub(
#     "trial",
#     paste0(seed_trial, "_", query_trial),
#     colnames(territories)
#   )


#   integrated <- lapply(integrated, function(x, locs) {
#     x <- x[, colnames(x) %in% locs]
#     return(x)
#   }, locs = c(seed_locs, query_locs))
#   new_assay <- vesalius:::build_vesalius_assay(
#     coordinates = coordinates,
#     verbose = FALSE
#   )
#   new_assay <- vesalius:::update_vesalius_assay(new_assay,
#     data = integrated,
#     slot = "counts",
#     append = FALSE
#   )
#   new_assay <- vesalius:::update_vesalius_assay(new_assay,
#     data = territories,
#     slot = "territories",
#     append = FALSE
#   )
#   new_assay <- vesalius:::add_active_count_tag(new_assay,
#     norm = "data"
#   )
#   return(list("assay" = new_assay, "seed" = seed, "query" = query))
# }
