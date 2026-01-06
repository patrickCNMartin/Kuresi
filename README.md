# Kuresi

![R-CMD-check](https://github.com/patrickCNMartin/Kuresi/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)

**Competitive fitness of Cancer clones**

Kuresi is an R package for analyzing competitive fitness between cancer clones in single cell and spatial transcriptomics data. It provides tools to compute and visualize competition scores based on differential gene expression patterns.

## Overview

Kuresi analyzes competition between cancer cell populations by comparing expression patterns of "win" genes (associated with competitive fitness) and "lose" genes (associated with cell death). The package computes ratio-based competition scores and uses differential expression analysis to assess competitive outcomes between groups.

## Key Features

- **Competition Score Computation**: Calculate ratio-based scores between win and lose gene sets
- **Differential Expression Analysis**: Use statistical tests (Wilcoxon, t-test) to identify competitive outcomes
- **ELO Rating System**: Convert competition outcomes into ELO ratings for ranking clones
- **Spatial Visualization**: Visualize competition scores on spatial transcriptomics data
- **Gene Set Management**: Access predefined win/lose gene lists and cancer marker genes

## Installation

### From GitHub

```r
devtools::install_github("patrickCNMartin/Kuresi")
```
### From tarball

First, download `Kuresi_0.0.1.tar.gz` available in the root of this repository. Move the tarball to your desired location, and make
sure you are in that directory within R, then:

```r
install.packages("Kuresi_0.0.1.tar.gz", repo = NULL)
```


## Quick Start

```r
library(Kuresi)

# Get win/lose gene sets
gene_sets <- win_lose_genes()
win_genes <- gene_sets$win
lose_genes <- gene_sets$lose

# Compute ratio scores
scores <- compute_ratio_score(
  counts = count_matrix,
  genes_1 = win_genes,
  genes_2 = lose_genes,
  method = "mean_ratio"
)

# Compute competition outcomes between groups
outcomes <- compute_competition_outcomes(
  counts = count_matrix,
  groups = group_assignments,
  gene_set1 = win_genes,
  gene_set2 = lose_genes
)

# Convert to ELO ratings
elo_ratings <- outcomes_as_elo(outcomes)

# Visualize scores
score_plot(score_data)
view_scores(score_matrix)
view_elo(elo_ratings)
```

## Main Functions

- `compute_ratio_score()`: Compute competition scores based on gene set ratios
- `compute_competition_outcomes()`: Analyze competitive outcomes using differential expression
- `outcomes_as_elo()`: Convert outcomes to ELO ratings
- `outcomes_as_cost()`: Convert outcomes to cost matrix
- `outcomes_as_score_matrix()`: Convert outcomes to score matrix
- `win_lose_genes()`: Get predefined win/lose gene lists
- `cancer_maker_list()`: Get cancer-specific marker genes
- `score_plot()`: Visualize scores spatially
- `view_scores()`: Heatmap visualization of score matrices
- `view_elo()`: Bar plot of ELO ratings


## License

GPL-3
