################################################################################
message("10 # KSC SHARED SETUP")
################################################################################

set.seed(20240725)
options(scipen = 999)


#=========================== Script path resolution ============================#
message("20 # Script path resolution")

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


#=================================== END SETUP ===================================#


################################################################################
message("40 # LOAD HELPERS AND PACKAGES")
################################################################################

wd <- get_script_dir()
setwd(wd)
source(file.path(dirname(wd), "helperJ.R"))

required_libraries <- c("tidyverse", "readxl", "ggtext")
load_libraries(required_libraries)


#============================== Input tables ==============================#
message("55 # Input tables")

qLoad("data_tables/meta.csv")
meta_names <- names(meta)


#============================= Shared objects =============================#
message("65 # Shared objects")

taxa_levels <- c("Phylum", "Family", "Genus", "Species")
taxa_plural <- c("Phyla", "Families", "Genera", "Species")

output_plot <- "output_plot/"
create_directory(output_plot)

dpi <- 600
point_size <- 3.5
treatment_order <- c("No PMA", "PMA")
location_order <- c("PHSF Highbay floor", "Field Control", "PHSF Airlock floor")


#=========================== Plot defaults ===========================#
message("80 # Plot defaults")

palette_label <- c(
  "Clipper" = "Clipper",
  "PHSF Highbay floor" = "PHSF Highbay floor",
  "Field Control" = "Field Control",
  "PHSF Airlock floor" = "PHSF Airlock floor",
  "No PMA" = "No PMA",
  "PMA" = "PMA",
  palette_helper_label
)

palette_color <- c(
  "Clipper" = palette_common[1],
  "PHSF Highbay floor" = palette_common[2],
  "Field Control" = palette_common[3],
  "PHSF Airlock floor" = palette_common[4],
  "PMA" = palette_blue[4],
  "No PMA" = palette_red[4],
  "panel_outline" = "black",
  "general_background" = "white",
  "panel_background" = "white",
  "plot_outline" = "white",
  "default" = "white",
  palette_helper_color
)

palette_shape <- c(
  "PHSF Highbay floor" = 22,
  "Field Control" = 24,
  "PHSF Airlock floor" = 21,
  "default" = 21
)

theme_global <- function(base_size = 11) {
  theme_gray(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(
        fill = palette_color["general_background"],
        color = palette_color["general_background"],
        linewidth = 1
      ),
      panel.background = element_rect(
        fill = palette_color["panel_background"],
        color = palette_color["panel_outline"],
        linewidth = 1
      ),
      legend.background = element_rect(fill = palette_color["general_background"]),
      legend.key = element_rect(color = "black", linewidth = 0.3),
      legend.text = element_text(size = rel(1.0)),
      legend.title = ggtext::element_markdown(size = rel(1.2), face = "bold"),
      axis.title = ggtext::element_markdown(size = rel(1.2), face = "bold"),
      axis.text = ggtext::element_markdown(size = rel(1.0), face = "bold", color = "black"),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = rel(1.1), face = "bold"),
      plot.title = ggtext::element_markdown(size = rel(1.4), face = "bold", hjust = 0),
      plot.caption = element_text(size = rel(1.0), face = "bold")
    )
}


################################################################################
message("120 # END KSC SHARED SETUP")
################################################################################
