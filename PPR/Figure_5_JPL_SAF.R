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
Creates the Figure 5 JPL-SAF comparison plots for the PPR dataset.
Builds the Figure 5A and 5B panels, writes the selection table, and saves the
PNG, PDF, and RDS outputs used for review.
"

#=============================== Plot settings ================================#

options(scipen = 999)

script_title <- "Figure_5"
figure_number <- "5"
plot_width <- 6
plot_height <- 8

#=================================== END SETUP ===================================#

################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

theme_common <- theme_global(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.y = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown()
  )

gPlot <- function(p) {
  p <- p +
    theme_common +
    scale_fill_manual(values = palette_color, labels = palette_label) +
    scale_color_manual(values = palette_color, labels = palette_label) +
    scale_linetype_manual(values = palette_linetype, labels = palette_label) +
    scale_shape_manual(values = palette_shape, labels = palette_label) +
    scale_y_continuous(labels = log10_labels_bold, expand = expansion(mult = c(0, 0.1))) +
    scale_x_discrete(labels = palette_label)

  print(p)
  p
}

build_summary_results <- function(df, facet_columns, value_column = "value_non_log10") {
  summary_brackets(
    df = df,
    axis_column = "sampling_device",
    value_column = value_column,
    facet_columns = facet_columns,
    testing_col = "sampling_device",
    test_type = NULL,
    adjust_method = "BH",
    width = 0.8,
    log = FALSE,
    space_adj = 1,
    tip_adj = 1,
    run_tests = TRUE
  )
}

################################################################################
message("90 # MAIN")
################################################################################

#.............................. prepare input data .............................#

df <- meta %>%
  filter(str_detect(Figure, figure_number)) %>%
  select(
    sample_name,
    sample_id,
    copies_original_qpcr_log10,
    copies_original_dpcr_log10,
    copies_original_qpcr,
    copies_original_dpcr,
    sampling_device,
    figure,
    Figure
  )

figure_5_selection <- df %>%
  mutate(
    figure_5 = case_when(
      sampling_device %in% c("cotton", "macrofoam") ~ "5A",
      sampling_device %in% c("SALSA", "Polyester Wipe") ~ "5B",
      TRUE ~ "5"
    )
  ) %>%
  select(sample_id, figure_5)

write.csv(figure_5_selection, "data_tables/figure_5_selection.csv", row.names = FALSE)

df_non_log10 <- df %>%
  mutate(
    copies_original_dpcr = as.numeric(copies_original_dpcr),
    copies_original_qpcr = as.numeric(copies_original_qpcr)
  ) %>%
  pivot_longer(
    cols = c(copies_original_qpcr, copies_original_dpcr),
    names_to = "method",
    values_to = "value_non_log10"
  ) %>%
  mutate(method = str_replace(method, "copies_original_", "")) %>%
  select(sample_name, method, value_non_log10) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  distinct()

df_plot <- df %>%
  pivot_longer(
    cols = c(copies_original_qpcr_log10, copies_original_dpcr_log10),
    names_to = "method",
    values_to = "value"
  ) %>%
  mutate(
    method = str_replace(method, "copies_original_", ""),
    method = str_replace(method, "_log10", "")
  ) %>%
  distinct() %>%
  left_join(df_non_log10, by = c("sample_name", "method")) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  distinct()

################################################################################
message("180 # BUILD PANEL A")
################################################################################

df_plot_all <- df_plot %>%
  filter(figure %in% c("5A", "5B"))

panel_tests_all <- build_summary_results(
  df = df_plot_all,
  facet_columns = c("figure", "method")
)$tests

df_plot_A <- df_plot %>%
  filter(figure == "5A")

panel_tests_A <- build_summary_results(
  df = df_plot_A,
  facet_columns = "method"
)$tests

panel_summary_A <- build_summary_results(
  df = df_plot_A,
  facet_columns = "method",
  value_column = "value"
)

df_summary_A <- panel_summary_A$summary
df_segments_A <- panel_summary_A$segments

df_tests_A <- panel_tests_all %>%
  filter(figure == "5A") %>%
  select(-figure) %>%
  left_join(
    df_segments_A %>%
      filter(part == "top") %>%
      select(method, y),
    by = "method"
  ) %>%
  mutate(
    label_value = y,
    adj.significance = mann_log_sig,
    adj.p.value = mann_log_val
  )

plot_title <- "Figure_5A_JPL_SAF"

p_A <- ggplot(df_plot_A, aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_jitter(width = 0.2, size = point_size, alpha = 0.8, shape = 21) +
  geom_segment(
    data = df_segments_A,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  facet_grid(
    method ~ .,
    scale = "free_x",
    labeller = labeller(.cols = palette_label, .rows = palette_label)
  ) +
  geom_text(
    data = df_tests_A,
    aes(x = label_position, y = label_value, label = adj.significance),
    inherit.aes = FALSE,
    size = 4,
    vjust = -0.5
  ) +
  geom_text(
    data = df_tests_A,
    aes(x = label_position, y = label_value, label = scales::scientific(adj.p.value)),
    inherit.aes = FALSE,
    size = 4,
    vjust = 1.5
  ) +
  geom_errorbar(
    data = df_summary_A,
    aes(x = sampling_device, ymin = mean_minus_sd, ymax = mean_plus_sd),
    width = 0.1,
    linewidth = 1,
    inherit.aes = FALSE,
    color = bar_color
  ) +
  geom_errorbar(
    data = df_summary_A,
    aes(x = sampling_device, ymin = mean_value, ymax = mean_value),
    width = 0.4,
    linewidth = 1.5,
    inherit.aes = FALSE,
    color = bar_color
  ) +
  coord_cartesian(ylim = zoom) +
  labs(
    fill = "",
    color = "",
    title = "A) JPL-Spacecraft Assembly Facility<br><span style=\"color:white\">A) </span>(25 cm<sup>2</sup>)",
    y = "log10(16S rRNA gene copies) / 25 cm<sup>2</sup>",
    x = ""
  )

plotA <- gPlot(p_A)

qSave(plot_title, ext = c("png", "pdf", ".rds"))

################################################################################
message("420 # FINISHED PANEL A")
################################################################################

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
