# Kuresi

![R-CMD-check](https://github.com/patrickCNMartin/Kuresi/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)

**Competitive fitness of Cancer clones**

## Overview

Kuresi is an R package for analyzing competitive fitness between cancer clones in single cell and spatial transcriptomics data. It computes and visualizes competition scores based on differential gene expression patterns, comparing "win" genes (associated with competitive fitness) against "lose" genes (associated with cell death). Outcomes can be represented as ELO ratings, cost matrices, or score matrices.

## Installation

### From GitHub

```r
devtools::install_github("patrickCNMartin/Kuresi")
```

### From tarball

Download `Kuresi_0.0.1.tar.gz` from the root of this repository, move it to your desired location, ensure you are in that directory within R, then:

```r
install.packages("Kuresi_0.0.1.tar.gz", repo = NULL)
```

## Development

This repository includes a Nix flake that provides a fully reproducible development environment with pinned R and package versions. This is the recommended way to contribute.

### Prerequisites

- [Nix](https://nixos.org/download/) with flakes enabled

### Getting started

```bash
git clone https://github.com/patrickCNMartin/Kuresi.git
cd Kuresi
nix develop
```

This drops you into a shell with R and all dependencies available at pinned versions. From there, use the standard R development workflow:

```r
# Install dev dependencies and load the package
devtools::load_all()

# Run tests
devtools::test()

# Build documentation
devtools::document()
```

### Enabling flakes

If you haven't enabled Nix flakes yet, add the following to `~/.config/nix/nix.conf`:

```
experimental-features = nix-flakes
```

## License

GPL-3