################################################################################
# clean work space, set a seed
################################################################################
rm(list = ls())
set.seed(20251117) # Date
################################################################################
# Description
################################################################################
# Use to set up github
################################################################################
################################################################################
# set working directory from script location
# load helper file/functions
################################################################################
wd <-dirname(normalizePath(rstudioapi::getSourceEditorContext()$path))
project_folder <- basename(wd)
setwd(wd)
getwd()

source("helperJ.R")
################################################################################
# This is in new project folder
################################################################################
notes='
remove all large files, do minimal setup first time
max upload < 500 mB
max file size < 100 mb
'

################################################################################
# After project creation, create github repository
################################################################################
full_manuscript_name='Advancing Spacecraft Microbial Monitoring through Absolute Quantification and Shotgun Metagenomic Diversity Profiling'

notes=paste('
set repository name to', project_folder,'
do not check: readme,gitignore,license

create repository
')

description=paste('R-scripts and data for',full_manuscript_name)
description
################################################################################
# create text to run in terminal'
################################################################################
github_user <- "JayKahsen"  

proj_dir  <- getwd()
repo_name <- basename(proj_dir)
remote    <- sprintf("https://github.com/%s/%s.git", github_user, repo_name)

cmds <- c(
  sprintf('cd "%s"', proj_dir),
  
  # Initialize the repository (only needed once)
  "git init",
  
  # Set remote (only needed once)
  sprintf("git remote add origin %s", remote),
  
  # Stage and commit all files
  "git add .",
  'git commit -m "Initial commit"',
  
  # Ensure branch name is main
  "git branch -M main",
  
  # FIRST push — registers the upstream branch
  "git push --set-upstream origin main"
)

# Write to a file you can open & copy from
writeLines(cmds, "git_setup_commands.txt")

# Also print to the R console for quick copy-paste
cat(paste(cmds, collapse = "\n"))


notes='
copy paste with mouse
hit enter to run last line
'

################################################################################
# updating'
# git push --set-upstream origin main
################################################################################

update_description='added readme.txt'

update_cmds <- c(
  'git add .',
  paste0('git commit -m "',update_description,'"'),
  'git push'
  
)


cat(paste(update_cmds, collapse = "\n"))




