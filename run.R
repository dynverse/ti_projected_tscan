#!/usr/local/bin/Rscript

task <- dyncli::main()

library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(TSCAN)
library(igraph)

#   ____________________________________________________________________________
#   Load data                                                               ####

counts <- task$counts
parameters <- task$parameters

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# preprocess counts
cds_prep <- TSCAN::preprocess(
  t(as.matrix(counts)),
  takelog = TRUE,
  logbase = 2,
  pseudocount = 1,
  clusternum = NULL,
  minexpr_value = parameters$minexpr_value,
  minexpr_percent = parameters$minexpr_percent,
  cvcutoff = parameters$cvcutoff
)

# cluster the data
cds_clus <- TSCAN::exprmclust(
  cds_prep,
  clusternum = parameters$clusternum,
  modelNames = parameters$modelNames,
  reduce = TRUE
)

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

# process output
milestone_network <- cds_clus$MSTtree %>%
  igraph::as_data_frame() %>%
  rename(length = weight) %>%
  mutate(directed = FALSE)
dimred <- cds_clus$pcareduceres
dimred_milestones <- cds_clus$clucenter
rownames(dimred_milestones) <- as.character(seq_len(nrow(dimred_milestones)))
colnames(dimred_milestones) <- colnames(dimred)

#   ____________________________________________________________________________
#   Save output                                                             ####

output <- dynwrap::wrap_data(cell_ids = rownames(dimred)) %>%
  dynwrap::add_dimred_projection(
    milestone_network = milestone_network,
    dimred = dimred,
    dimred_milestones = dimred_milestones,
    grouping = cds_clus$clusterid
  ) %>%
  dynwrap::add_timings(timings = checkpoints)

output %>% dyncli::write_output(task$output)
