#Load libs
library(devtools)
library(rmarkdown)
library(roxygen2)
roxygenise('.')
detach(package:Kuresi, unload = TRUE)
detach(package:vesalius, unload = TRUE)
system('R CMD build --no-build-vignettes --resave-data .')
install.packages('Kuresi_0.0.1.tar.gz', repo = NULL)
library(Kuresi)
library(vesalius)
message('Successfully removed, rebuilt, re-installed and reloaded Kuresi')
#message('Calculating test unit coverage of Kuresi')
# setwd("Kuresi")
# #try(report())
# setwd('../')

# max_size <- 10000 * 1024^2
# options(future.globals.maxSize = max_size)
# win_lose <- win_lose_genes()
# win <- win_lose$win
# lose <- win_lose$lose
# pancreas <- cancer_maker_list(type = "pancreas")
# breast <- cancer_maker_list(type = "breast")
# ovary <- cancer_maker_list(type = "ovary")
# vesalius <- readRDS("data/ovarian.rds")
# vesalius_2 <- readRDS("data/ovarian.rds")
# vesalius_2@tiles$barcodes <- paste0("q_",vesalius_2@tiles$barcodes)
# rownames(vesalius_2@embeddings$PCA) <- paste0("q_",rownames(vesalius_2@embeddings$PCA))
# rownames(vesalius_2@active) <- paste0("q_",rownames(vesalius_2@active))
# vesalius_2@territories$barcodes <- paste0("q_",vesalius_2@territories$barcodes)
# colnames(vesalius_2@counts$raw) <- paste0("q_",colnames(vesalius_2@counts$raw))
# colnames(vesalius_2@counts$log_norm) <- paste0("q_",colnames(vesalius_2@counts$log_norm))
# vesalius_2@meta$orig_coord$barcodes <- paste0("q_",vesalius_2@meta$orig_coord$barcodes)


# outcomes <- compute_competition_outcomes(vesalius, win, lose, vesalius_2)

#message('Generating vignette Kuresi')
#rmarkdown::render('Kuresi/vignettes/Kuresi', 'html_document')


ter1 <- coord$barcodes[coord$cell_labels %in% c(1.0, 1.1)]
c1 <- log2(as.matrix(counts[, ter1]) + 1)
ter2 <- coord$barcodes[!coord$cell_labels %in% c(1.0, 1.1)]
c2 <- log2(as.matrix(counts[, ter2]) + 1)

test <- Kuresi:::k_wilcox(c1,c2)
test <- test[test$p_value_adj < 0.05 &
             test$fold_change > 0.5 |
             test$fold_change < -0.5,]
gene_set <- test$genes[order(test$fold_change, decreasing = TRUE)]
gene_set1 <- head(gene_set)
gene_set2 <- tail(gene_set)
rownames(coord) <- coord$barcodes
ratio <- compute_ratio_bygroup(counts,coord,"cell_labels", gene_set1,gene_set2) 
