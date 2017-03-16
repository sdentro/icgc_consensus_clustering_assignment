args = commandArgs(T)
samplename = args[1]

source("~/repo/moritz_mut_assignment/MutationTime.R")
source("~/repo/moritz_mut_assignment/util.R")

load(paste0(samplename, "_assignment.RData"))


##############################################
#' Get summary table status




