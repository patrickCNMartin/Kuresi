library(rix)

rix(
  r_ver = "latest-upstream",
  r_pkgs = c(
  "ggplot2",
  "patchwork",
  "ggpubr",
  "rjson",
  "Matrix",
  "hdf5r",
  "Seurat",
  "RColorBrewer",
  "knitr",
  "RUnit",
  "rmarkdown",
  "testthat",
  "roxygen2",
  "devtools",
  "rix",
  "rlang",
  "styler",
  "lintr"),
  system_pkgs = NULL,
  git_pkgs = list(
    package_name = "Vesalius",
    repo_url = "https://github.com/WonLab-CS/Vesalius",
    commit = "bb4aa3d"
  ),
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)