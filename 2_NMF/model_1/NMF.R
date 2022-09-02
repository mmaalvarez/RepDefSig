library(tidyverse)
library(data.table)
library(dtplayer)
#devtools::install_github("zdebruine/RcppML")
library(RcppML)


results_regressions = fread("../1_parser_and_regressions/res/results.tsv")

# very fast NMF --> could not install it in the R conda environment in agendas, ask
RcppML::nmf()
