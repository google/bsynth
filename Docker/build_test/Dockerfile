FROM rocker/r-base:latest
LABEL maintainer="Ignacio Martinez <martinezig@google.com>"

RUN apt-get update -qq \
 && apt-get install -y \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libv8-dev \
        texlive-latex-extra \
	libharfbuzz-dev \
	libfribidi-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev \
	libpq5 \
	pandoc \
 && install.r \
  magrittr \
  R6 \
  Rcpp \
  cubelyr \
  dplyr \
  ggplot2 \
  glue \
  magrittr \
  methods \
  purrr \
  rstan \
  rstantools \
  scales \
  tibble \
  tidyr \
  vizdraws \
  devtools \
  BH \
  RcppEigen \
  methods \
  rmarkdown \
  vizdraws \
  knitr \
  rmarkdown \
  roxygen2 \
  gsynth
  
