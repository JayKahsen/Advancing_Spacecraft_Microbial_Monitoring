################################################################################
# clean work space, set a seed, avoid scientific nptation
################################################################################
rm(list = ls())
set.seed(20240725) # Date
options(scipen = 999)
################################################################################
# Description
################################################################################
# Global Script used as a starting point for other scripts
################################################################################

################################################################################
# set working directory from script location
# load helper file/functions
# load libraries
################################################################################
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
  }
  getwd()
}

wd <- get_script_dir()
project_folder <- basename(wd)
setwd(wd)
source("helperJ.R")
