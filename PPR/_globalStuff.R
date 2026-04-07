################################################################################
message("10 # PPR SHARED SETUP")
################################################################################

set.seed(20240725)
options(scipen = 999)

#============================= Script description =============================#

script_title <- "global"

script_description <- "
Loads shared packages, helper functions, input tables, plotting defaults,
ordering vectors, palettes, and theme objects used by the PPR scripts.
"


#=========================== Script path resolution ============================#

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
project_folder <- basename(wd)
setwd(wd)
source(file.path(dirname(wd), "helperJ.R"))

required_libraries <- c(
  "tidyverse",
  "patchwork",
  "ggtext",
  "vegan",
  "ggh4x",
  "ggpmisc"
)
load_libraries(required_libraries)


################################################################################
message("70 # INPUT TABLES")
################################################################################

current_date <- format(Sys.Date(), "%Y%m%d")

qLoad("data_tables/original/Contaminant_Id_Taxa.xlsx")
qLoad("data_tables/meta.csv")
qLoad("data_tables/matrix_names.csv")
qLoad("data_tables/ASV_taxa.csv")

if (exists("meta")) {
  meta_names <- names(meta)
  qPrint(meta_names)
}


################################################################################
message("100 # SHARED OBJECTS")
################################################################################

taxa_levels <- c("Phylum", "Family", "Genus", "Species")
taxa_plural <- c("Phyla", "Families", "Genera", "Species")

if (exists("matrix_names")) {
  data_set_order <- data_sets <- unique(matrix_names$data_sets)
}

output_data <- "output_data/"
output_plot <- "output_plot/"
dpi <- 300

treatment_order <- c("No PMA", "PMA")
extraction_method_order <- c("mag_beads", "spin_column")
location_order <- c("PHSF Highbay floor", "Field Control", "PHSF Airlock floor")
experiment_order <- c("PPR")
sample_type_order <- c("reagent_control", "2.00E+1", "2.00E+2", "2.00E+3", "2.00E+4", "2.00E+5", "2.00E+6", "2.00E+7", "water_control", "positive_control", "no_template_control", "100000", "NTC", "NA")
bar_name_order <- c("Cotton Device Control", "Macrofoam Device Control", "Cotton Environmental Control", "Macrofoam Environmental Control", "Cotton Dry + Metal", "Macrofoam Dry + Metal", "Cotton Wet + Metal", "Macrofoam Wet + Metal", "Water Only", "NTC")
bar_color_order <- c("Cotton", "Macrofoam", "Water", "NTC")
bar_axis_2_order <- c("Cotton", "Macrofoam", "Water", "NTC")
sampling_device_order <- c("cotton", "macrofoam", "Polyester Wipe", "SALSA")
bar_type_order <- c("Mastermix", "Environmental", "Device", "Metal")
bar_type_2_order <- c("Mastermix", "Environmental", "Device", "Metal Dry", "Metal Wet")
method_order <- c("qpcr", "dpcr")
experiment_type_order <- c("metal_deposition", "swab_head_retention")
figure_order <- c("2", "3A", "3B", "4A", "4B", "4C", "5")
order_vector <- unique(c(
  treatment_order,
  experiment_order,
  location_order,
  bar_color_order,
  extraction_method_order,
  sampling_device_order,
  bar_type_order,
  method_order,
  experiment_type_order
))
custom_name <- ""

atcc_genera <- c(
  "Bacillus",
  "Bifidobacterium",
  "Clostridium",
  "Deinococcus",
  "Enterococcus",
  "Escherichia",
  "Lactobacillus",
  "Cereibacter",
  "Staphylococcus",
  "Streptococcus"
)


################################################################################
message("140 # PLOT DEFAULTS")
################################################################################

log_limits <- c(NA, NA)
zoom <- c(NA, NA)
point_size <- 2.5
bar_color <- "gray30"
magnify <- 1
text_size <- 10
margin_size <- 10

palette_label <- c(
  "<i>Escherichia</i> " = "<i>Escherichia</i> ",
  "<i>Cereibacter</i> " = "<i>Cereibacter</i> ",
  "<i>Bacillus</i> " = "<i>Bacillus</i> ",
  "<i>Lactobacillus</i> " = "<i>Lactobacillus</i> ",
  "<i>Clostridium</i> " = "<i>Clostridium</i> ",
  "<i>Streptococcus</i> " = "<i>Streptococcus</i> ",
  "<i>Enterococcus</i> " = "<i>Enterococcus</i> ",
  "<i>Bifidobacterium</i> " = "<i>Bifidobacterium</i> ",
  "<i>Staphylococcus</i> " = "<i>Staphylococcus</i> ",
  "<i>Deinococcus</i> " = "<i>Deinococcus</i> ",
  "Other" = "Other",
  "reagent_control" = "Reagent<br>Control",
  "mag_beads" = "Automatic Magnetic Beads",
  "spin_column" = "Manual Spin Column",
  "5A" = "5A",
  "5B" = "5B",
  "ng100" = "<100",
  "g100" = ">100",
  "g1K" = ">1K",
  "g10K" = ">10K",
  "g100K" = ">100K",
  "metal_deposition" = "Recovery from Metal Surface",
  "swab_head_retention" = "Swab Head Released",
  "unique_features" = "Unique Features",
  "chimeras" = "Chimeras",
  "less_than_0" = "0",
  "log0" = "1",
  "log4" = "e4",
  "log6" = "e6",
  "log8" = "e8",
  "log10" = "e10",
  "0P" = "0%",
  "20P" = "10%",
  "40P" = "20%",
  "60P" = "40%",
  "80P" = "80%",
  "g0" = "0",
  "g10" = "10",
  "g20" = "20",
  "g40" = "40",
  "g80" = "80",
  "4A" = "A\\) DNA Released from Swab Heads",
  "4B" = "B\\) DNA Recovered by Swabs",
  "4C" = "C\\) DNA Rec. by Polyester Wipes and SALSA",
  "No PMA" = "No PMA",
  "PMA" = "PMA",
  "PPR" = "PPR",
  "Clipper" = "Clipper",
  "PHSF Highbay floor" = "PHSF Highbay floor",
  "Field Control" = "Field Control",
  "PHSF Airlock floor" = "PHSF Airlock floor",
  "ATCC" = "ATCC",
  "Cotton Wet" = "Cotton Wet",
  "Macrofoam Wet" = "Macrofoam Wet",
  "Cotton Dry" = "Cotton Dry",
  "Macrofoam Dry" = "Macrofoam Dry",
  "Cotton" = "Cotton",
  "Macrofoam" = "Macrofoam",
  "Water" = "Water",
  "NTC" = "NTC",
  "positive_control" = "Positive Control",
  "negative_control" = "Negative Control",
  "water_control" = "Water Control",
  "non_template_control" = "NTC",
  "no_template_control" = "NTC",
  "cotton" = "Cotton",
  "macrofoam" = "Macrofoam",
  "Polyester Wipe" = "Polyester Wipe",
  "SALSA" = "SALSA",
  "Mastermix" = "Mastermix",
  "Environmental" = "Environmental",
  "Device" = "Swab",
  "Metal" = "Metal",
  "Metal Dry" = "Metal Dry",
  "Metal Wet" = "Metal Wet",
  "Cotton Device Control" = "Cotton Device Control",
  "Macrofoam Device Control" = "Macrofoam Device Control",
  "Cotton Environmental Control" = "Cotton Environmental Control",
  "Macrofoam Environmental Control" = "Macrofoam Environmental Control",
  "Cotton Dry + Metal" = "Cotton Dry + Metal",
  "Macrofoam Dry + Metal" = "Macrofoam Dry + Metal",
  "Cotton Wet + Metal" = "Cotton Wet + Metal",
  "Macrofoam Wet + Metal" = "Macrofoam Wet + Metal",
  "Water Only" = "Water Only",
  "Positive Control" = "Positive Control",
  "Negative Control" = "Negative Control",
  "Reagent Control" = "Reagent<br>Control",
  "2.00E+2" = "2.00E+2",
  "2.00E+3" = "2.00E+3",
  "2.00E+4" = "2.00E+4",
  "2.00E+5" = "2.00E+5",
  "2.00E+6" = "2.00E+6",
  "ATP_RLU" = "ATP RLU",
  "Fungal_CFU" = "Fungal CFU",
  "Fungal_dPCR" = "Fungal dPCR",
  "Bacterial_CFU" = "Bacterial CFU",
  "Bacterial_dPCR" = "Bacterial dPCR",
  "gray1" = "gray1",
  "gray2" = "gray2",
  "gray3" = "gray3",
  "gray4" = "gray4",
  "gray5" = "gray5",
  "black" = "black",
  "default" = "default",
  palette_helper_label
)

palette_gray <- c("#BBBBBB", "#999999", "#777777", "#555555", "#333333")
palette_green <- c("#D2EDB2", "#89C189", "#60A740", "#107810", "#004B00")
palette_gold <- c("#FFDA65", "#FFCB25", "#DAA520", "#A37B18", "#7D6220")
palette_purple <- c("#D2AEFA", "#A267E8", "#7A1FD2", "#5B17A1", "#3C1070")
palette_red <- c("#FFCCCC", "#E47C7C", "#C94141", "#993333", "#6A2B2B")
palette_blue <- c("#C8D6FF", "#92AEFF", "#5C86FF", "#3A62BE", "#20407F")

assign_palette(order_vector)

palette_color <- c(
  "<i>Escherichia</i> " = palette_common[1],
  "<i>Cereibacter</i> " = palette_common[2],
  "<i>Bacillus</i> " = palette_common[3],
  "<i>Lactobacillus</i> " = palette_common[4],
  "<i>Clostridium</i> " = palette_common[5],
  "<i>Streptococcus</i> " = palette_common[6],
  "<i>Enterococcus</i> " = palette_common[7],
  "<i>Bifidobacterium</i> " = palette_common[8],
  "<i>Staphylococcus</i> " = palette_common[9],
  "<i>Deinococcus</i> " = palette_common[10],
  "Other" = palette_gray[1],
  "metal_deposition" = palette_gold[4],
  "swab_head_retention" = palette_blue[4],
  "spin_column" = palette_blue[3],
  "mag_beads" = palette_red[4],
  "macrofoam" = palette_green[2],
  "Macrofoam" = palette_green[2],
  "cotton" = palette_purple[1],
  "Cotton" = palette_purple[1],
  "Water" = palette_blue[4],
  "SALSA" = palette_red[4],
  "Polyester Wipe" = palette_green[4],
  "PMA" = palette_blue[4],
  "No PMA" = palette_red[4],
  "positive_control" = palette_blue[2],
  "negative_control" = palette_red[2],
  "reagent_control" = palette_purple[5],
  "NTC" = "grey40",
  "alt_group1" = palette_gray[1],
  "alt_group2" = palette_gray[3],
  "Clipper" = palette_common[1],
  "PHSF Highbay floor" = palette_common[2],
  "Field Control" = palette_common[3],
  "PHSF Airlock floor" = palette_common[4],
  "ATCC" = palette_common[5],
  "2.00E+2" = palette_common[22],
  "2.00E+3" = palette_common[23],
  "2.00E+4" = palette_common[24],
  "2.00E+5" = palette_common[25],
  "2.00E+6" = palette_common[26],
  "gray1" = palette_gray[1],
  "gray2" = palette_gray[2],
  "gray3" = palette_gray[3],
  "gray4" = palette_gray[4],
  "gray5" = palette_gray[5],
  "Cotton Wet" = palette_purple[5],
  "Macrofoam Wet" = palette_blue[3],
  "Cotton Dry" = palette_purple[5],
  "Macrofoam Dry" = palette_blue[3],
  "PPR" = palette_common[3],
  "Mastermix" = palette_common[17],
  "Environmental" = palette_common[18],
  "Device" = palette_common[19],
  "Metal" = palette_common[20],
  "Metal Dry" = palette_common[20],
  "Metal Wet" = palette_common[20],
  "dpcr" = palette_common[21],
  "qpcr" = palette_common[22],
  "Cotton Device Control" = palette_purple[1],
  "Macrofoam Device Control" = palette_green[2],
  "Cotton Environmental Control" = palette_purple[1],
  "Macrofoam Environmental Control" = palette_green[2],
  "Cotton Dry + Metal" = palette_purple[1],
  "Macrofoam Dry + Metal" = palette_green[2],
  "Cotton Wet + Metal" = palette_purple[1],
  "Macrofoam Wet + Metal" = palette_green[2],
  "Water Only" = palette_blue[4],
  "Positive Control" = palette_blue[2],
  "Negative Control" = palette_red[2],
  "Reagent Control" = palette_purple[5],
  "ATP_RLU" = "#DC0000FF",
  "Fungal_CFU" = "#3C5488FF",
  "Fungal_dPCR" = "#4DBBD5FF",
  "Bacterial_CFU" = "#00A087FF",
  "Bacterial_dPCR" = "#91D1C2FF",
  "unique_features" = palette_common[18],
  "chimeras" = palette_common[19],
  "ns" = "gray50",
  "sig" = "red",
  "panel_outline" = "black",
  "general_background" = "white",
  "panel_background" = "white",
  "plot_outline" = "white",
  "black" = "black",
  "default" = "white",
  palette_helper_color
)

palette_linetype <- c(
  "ATP_RLU" = "solid",
  "Fungal_CFU" = "solid",
  "Fungal_dPCR" = "dotted",
  "Bacterial_CFU" = "solid",
  "Bacterial_dPCR" = "solid",
  "default" = "solid"
)

assign_default(order_vector, 21)

palette_shape <- c(
  "metal_deposition" = 23,
  "swab_head_retention" = 21,
  "Device" = 23,
  "PHSF Highbay floor" = 22,
  "Field Control" = 24,
  "PHSF Airlock floor" = 21,
  "ATCC" = 23,
  "2.00E+2" = 21,
  "2.00E+3" = 21,
  "2.00E+4" = 21,
  "2.00E+5" = 21,
  "2.00E+6" = 21,
  "NTC" = 21,
  "magnetic_bead" = 21,
  "spin_column" = 23,
  "ATP_RLU" = 17,
  "Fungal_CFU" = 15,
  "Fungal_dPCR" = 19,
  "Bacterial_CFU" = 15,
  "Bacterial_dPCR" = 22,
  "PMA" = 23,
  "No PMA" = 21,
  "PPR" = 21,
  "Clipper" = 21,
  "Cotton Wet" = 21,
  "Macrofoam Wet" = 21,
  "Cotton Dry" = 21,
  "Macrofoam Dry" = 21,
  "Cotton" = 21,
  "Macrofoam" = 21,
  "Water" = 21,
  "mag_beads" = 21,
  "cotton" = 21,
  "macrofoam" = 21,
  "Polyester Wipe" = 21,
  "SALSA" = 21,
  "Mastermix" = 21,
  "Environmental" = 21,
  "Metal" = 21,
  "Metal Dry" = 21,
  "Metal Wet" = 21,
  "dpcr" = 21,
  "qpcr" = 21,
  "default" = 21
)

assign_default(order_vector, 3)

palette_size <- c(
  "metal_deposition" = 3,
  "swab_head_retention" = 3,
  "0P" = 1,
  "20P" = 3,
  "40P" = 6,
  "60P" = 9,
  "80P" = 12,
  "less_than_0" = 1,
  "log0" = 1,
  "log4" = 3,
  "log6" = 6,
  "log8" = 9,
  "log10" = 12,
  "ng100" = 1,
  "g100" = 2,
  "g1K" = 4,
  "g10K" = 6,
  "g100K" = 8,
  "g0" = 1,
  "g10" = 3,
  "g20" = 6,
  "g40" = 9,
  "g80" = 12,
  "gray1" = 3,
  "gray2" = 3,
  "gray3" = 3,
  "gray4" = 3,
  "gray5" = 3,
  "PMA" = 3,
  "No PMA" = 3,
  "NTC" = 3,
  "magnetic_bead" = 3,
  "spin_column" = 3,
  "SALSA" = 3,
  "Polyester Wipe" = 3,
  "Cotton" = 3,
  "Macrofoam" = 3,
  "Positive Control" = 3,
  "Water" = 3,
  "Cotton Device Control" = 3,
  "Macrofoam Device Control" = 3,
  "Cotton Environmental Control" = 3,
  "Macrofoam Environmental Control" = 3,
  "Cotton Dry + Metal" = 3,
  "Macrofoam Dry + Metal" = 3,
  "Cotton Wet + Metal" = 3,
  "Macrofoam Wet + Metal" = 3,
  "Water Only" = 3,
  "ATP_RLU" = 3,
  "Fungal_CFU" = 3,
  "Fungal_dPCR" = 3,
  "Bacterial_CFU" = 3,
  "Bacterial_dPCR" = 3,
  "PPR" = 3,
  "Clipper" = 3,
  "PHSF Highbay floor" = 3,
  "Field Control" = 3,
  "PHSF Airlock floor" = 3,
  "mag_beads" = 3,
  "cotton" = 3,
  "macrofoam" = 3,
  "Mastermix" = 3,
  "Environmental" = 3,
  "Device" = 3,
  "Metal" = 3,
  "Metal Dry" = 3,
  "Metal Wet" = 3,
  "dpcr" = 3,
  "qpcr" = 3,
  "default" = 3
)

theme_global <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.tag = element_text(size = rel(2), face = "bold"),
      strip.background = element_rect(fill = NA, color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
      plot.background = element_rect(fill = palette_color["general_background"], color = palette_color["general_background"], linewidth = 1),
      panel.background = element_rect(fill = palette_color["panel_background"], color = palette_color["panel_outline"], linewidth = 1),
      axis.ticks = element_line(linewidth = 0.4),
      legend.background = element_rect(fill = palette_color["general_background"]),
      axis.title = element_markdown(size = rel(1.4), face = "bold"),
      axis.text = element_markdown(size = rel(1.2), face = "bold", color = "black"),
      legend.title = element_markdown(size = rel(1.2), face = "bold"),
      legend.text = element_text(size = rel(1.2)),
      legend.key = element_rect(color = "black", linewidth = 0.3),
      strip.text = element_text(size = rel(1.2), face = "bold"),
      plot.title = element_markdown(size = rel(1.6), face = "bold", hjust = 0),
      plot.subtitle = element_text(size = rel(1.2)),
      plot.caption = element_text(size = rel(1.2), face = "bold")
    )
}

theme_common <- theme_global(base_size = 11)

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
