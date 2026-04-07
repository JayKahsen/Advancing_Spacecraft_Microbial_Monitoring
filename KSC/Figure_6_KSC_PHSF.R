################################################################################
message("10 # SETUP")
################################################################################


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

source(file.path(get_script_dir(), "_globalStuff.R"))


#============================= Script description =============================#

script_description <- "
Creates Figure_6_KSC_PHSF.
Builds the KSC-PHSF bacterial treatment comparison plot and writes
PNG, PDF, and RDS outputs to output_plot.
"


#=============================== Plot settings ================================#

script_title <- "Figure_6_KSC_PHSF"
plot_title <- "Figure_6A_KSC_PHSF"
plot_width <- 10
plot_height <- 8


#=================================== END SETUP ===================================#


################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

theme_common <- theme_global(base_size = 11) +
  theme(
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown()
  )

gPlot <- function(p) {
  p <- p +
    theme_common +
    labs(
      fill = "Treatment Type",
      shape = "Location"
    ) +
    scale_y_continuous(labels = log10_labels_bold, breaks = seq(-12, 12, by = 2)) +
    scale_x_discrete(labels = palette_label) +
    scale_color_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
    scale_fill_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
    scale_shape_manual(values = palette_shape, labels = palette_label)

  print(p)
  p
}


################################################################################
message("90 # MAIN")
################################################################################


#........................... transform plot input data ..........................#

df_plot1 <- meta %>%
  mutate(value_non_log10 = value_bacterial) %>%
  mutate(
    value = ifelse(
      value_non_log10 < 0,
      NA,
      log10(ifelse(value_non_log10 == 0, 1, value_non_log10))
    )
  ) %>%
  filter(!is.na(value)) %>%
  filter(treatment %in% treatment_order) %>%
  filter(location %in% location_order) %>%
  ungroup()


#...................... summarize statistics and bracket positions .............#

results <- summary_brackets(
  df = df_plot1,
  axis_column = "treatment",
  value_column = "value",
  facet_columns = "location",
  testing_col = "treatment",
  test_type = NULL,
  adjust_method = "BH",
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE
)

df_summary <- results$summary
df_segments <- results$segments
df_tests <- results$tests %>%
  left_join(
    df_segments %>%
      filter(part == "top") %>%
      select(location, y)
  ) %>%
  mutate(label_value = y) %>%
  mutate(significance = mann_raw_sig) %>%
  mutate(p.value = mann_raw_val)


#.............................. build and save plot .............................#

p <- ggplot(df_plot1, aes(x = treatment, y = value, fill = treatment)) +
  geom_jitter(width = 0.2, size = point_size, alpha = 0, fill = "white") +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = df_tests,
    aes(x = label_position, y = label_value, label = significance),
    inherit.aes = FALSE,
    size = 3,
    vjust = -0.5
  ) +
  geom_text(
    data = df_tests,
    aes(x = label_position, y = label_value, label = scales::scientific(p.value)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 1.5
  ) +
  geom_errorbar(
    data = df_summary,
    aes(x = treatment, ymin = mean_value, ymax = mean_value, color = treatment),
    width = 0.6,
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  geom_jitter(aes(shape = location), width = 0.2, size = point_size) +
  labs(
    fill = "Treatment Type",
    shape = "Location",
    title = "A) KSC-PHSF facility surface- (m<sup>2</sup>) - Bacteria",
    y = "log10(16S rRNA gene copies) / m<sup>2</sup>",
    x = ""
  ) +
  guides(
    color = "none",
    alpha = "none",
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  facet_grid(~location)

p <- gPlot(p)
qSave(plot_title, ext = c("png", "pdf", "rds"))


################################################################################
message("170 # END MAIN")
################################################################################
