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
wd <-dirname(normalizePath(rstudioapi::getSourceEditorContext()$path))
project_folder <- basename(wd)
setwd(wd)
source("helperJ.R")