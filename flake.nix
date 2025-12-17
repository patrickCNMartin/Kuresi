{
  description = "Kuresi Development Shell";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/00f459c69df1086032022158cfe9f8f8530b1788";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        rpkgs = builtins.attrValues {
          inherit (pkgs.rPackages)
            devtools
            ggplot2
            ggpubr
            hdf5r
            knitr
            lintr
            Matrix
            patchwork
            RColorBrewer
            rix
            rjson
            rlang
            rmarkdown
            roxygen2
            RUnit
            Seurat
            styler
            testthat;
        };

        Vesalius = (pkgs.rPackages.buildRPackage {
          name = "Vesalius";
          src = pkgs.fetchgit {
            url = "https://github.com/WonLab-CS/Vesalius";
            rev = "bb4aa3d";
            sha256 = "sha256-Cx9CYwkIjdAlVQhUXR479sY1TsetZ0ExEj1z9TU95t8=";
          };
          propagatedBuildInputs = builtins.attrValues {
            inherit (pkgs.rPackages)
              deldir
              sp
              tvR
              Matrix
              RColorBrewer
              future
              future_apply
              imagerExtra
              Signac
              ggpubr
              lmtest
              infix
              Seurat
              SeuratObject
              imager
              DESeq2
              RANN
              dplyr
              edgeR
              ggplot2
              igraph
              pwr
              purrr
              patchwork
              ggnewscale
              kohonen
              Rcpp
              RcppEigen
              TreeDist;
          };
        });

        system_packages = builtins.attrValues {
          inherit (pkgs)
            glibcLocales
            nix
            R
            pandoc
            which;
        };
      in
      {
        devShells.default = pkgs.mkShell {
          LOCALE_ARCHIVE = if pkgs.stdenv.hostPlatform.system == "x86_64-linux" then "${pkgs.glibcLocales}/lib/locale/locale-archive" else "";
          LANG = "en_US.UTF-8";
          LC_ALL = "en_US.UTF-8";
          LC_TIME = "en_US.UTF-8";
          LC_MONETARY = "en_US.UTF-8";
          LC_PAPER = "en_US.UTF-8";
          LC_MEASUREMENT = "en_US.UTF-8";

          buildInputs = [ Vesalius ] ++ rpkgs ++ system_packages;
        };
      }
    );
}