% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{count_integration}
\alias{count_integration}
\title{Integrate counts using Seurat}
\usage{
count_integration(
  query,
  seed,
  integration_method = "cca",
  nfeatures = 2000,
  dimensions = 30,
  features = NULL
)
}
\arguments{
\item{nfeatures}{integer - number of features to use during integration}

\item{dimensions}{interger - number of dimensions integrated latent space
dimensions.}

\item{matched}{matrix - matrix containing counts from matched/mapped assay}

\item{reference}{matrix - matrix containing counts from reference assay}

\item{method}{character - Seurat integration method to use}

\item{infer}{logical - back infer original counts by reversing reduced 
dimensional space roations.}

\item{signal}{character - defining which signal should be returned:
variable_features, all_features or custom gene list.}
}
\description{
Integrate counts using Seurat
}
