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
Creates the Figure 5 dPCR-vs-qPCR comparison plots for the PPR dataset.
Builds the Figure 5A and 5B dqPCR panels, saves the combined figure, and writes
panel-specific PDF and RDS outputs to output_plot.
"

#=============================== Plot settings ================================#

options(scipen = 999)

script_title <- "alt_Figure_5_dqPCR"
figure_number <- "5"

#=================================== END SETUP ===================================#

################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

theme_common <- theme_global(base_size = 11) +
  theme(
    legend.position = "none"
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

sig_stars <- function(p) {
  case_when(
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

build_summary_results <- function(df) {
  list(
    non_log = summary_brackets(
      df = df,
      axis_column = "method",
      value_column = "value_non_log10",
      facet_columns = "sampling_device",
      testing_col = "method",
      test_type = NULL,
      adjust_method = "BH",
      width = 0.8,
      log = FALSE,
      space_adj = 1,
      tip_adj = 1,
      run_tests = TRUE
    ),
    log = summary_brackets(
      df = df,
      axis_column = "method",
      value_column = "value",
      facet_columns = "sampling_device",
      testing_col = "method",
      test_type = NULL,
      adjust_method = "BH",
      width = 0.8,
      log = FALSE,
      space_adj = 1,
      tip_adj = 1,
      run_tests = TRUE
    )
  )
}

build_wilcox_results <- function(df_panel) {
  df_panel %>%
    group_by(sampling_device) %>%
    do({
      df_facet <- .

      df_wide <- df_facet %>%
        select(sample_name, method, value) %>%
        pivot_wider(names_from = method, values_from = value)

      if (ncol(df_wide) != 3) {
        return(tibble(
          sampling_device = unique(df_facet$sampling_device),
          p.value = NA_real_
        ))
      }

      method_cols <- colnames(df_wide)[2:3]
      v1 <- df_wide[[method_cols[1]]]
      v2 <- df_wide[[method_cols[2]]]

      wtest <- wilcox.test(v1, v2, paired = TRUE, exact = FALSE)

      tibble(
        sampling_device = unique(df_facet$sampling_device),
        p.value = wtest$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      wilcox_val = p.value,
      wilcox_sig = sig_stars(wilcox_val)
    )
}

build_panel_plot <- function(df_panel, df_tests_all, figure_id, panel_title, y_title) {
  summary_results <- build_summary_results(df_panel)

  df_summary <- summary_results$log$summary
  df_segments <- summary_results$log$segments

  df_tests <- df_tests_all %>%
    filter(figure == figure_id) %>%
    select(-figure) %>%
    left_join(
      df_segments %>%
        filter(part == "top") %>%
        select(sampling_device, y),
      by = "sampling_device"
    ) %>%
    mutate(label_value = y)

  df_wilcox <- build_wilcox_results(df_panel)

  df_tests <- df_tests %>%
    left_join(df_wilcox %>% select(sampling_device, wilcox_val, wilcox_sig), by = "sampling_device")

  pos_dodge <- position_dodge(width = 0.6)

  p <- ggplot(df_panel, aes(x = method, y = value, fill = method)) +
    geom_line(
      aes(group = sample_name, color = NULL),
      position = pos_dodge,
      linewidth = 0.6,
      alpha = 0.2
    ) +
    geom_point(
      aes(group = sample_name),
      position = pos_dodge,
      size = point_size,
      alpha = 0.8,
      shape = 21
    ) +
    geom_segment(
      data = df_segments,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE
    ) +
    facet_grid(
      sampling_device ~ .,
      scale = "free_x",
      labeller = labeller(.cols = palette_label, .rows = palette_label)
    ) +
    geom_text(
      data = df_tests,
      aes(x = label_position, y = label_value, label = wilcox_sig),
      inherit.aes = FALSE,
      size = 4,
      vjust = -0.5
    ) +
    geom_text(
      data = df_tests,
      aes(x = label_position, y = label_value, label = scales::scientific(wilcox_val)),
      inherit.aes = FALSE,
      size = 4,
      vjust = 1.5
    ) +
    geom_errorbar(
      data = df_summary,
      aes(x = method, ymin = mean_minus_sd, ymax = mean_plus_sd),
      width = 0.1,
      linewidth = 1,
      inherit.aes = FALSE,
      color = bar_color
    ) +
    geom_errorbar(
      data = df_summary,
      aes(x = method, ymin = mean_value, ymax = mean_value),
      width = 0.4,
      linewidth = 1.5,
      inherit.aes = FALSE,
      color = bar_color
    ) +
    coord_cartesian(ylim = zoom) +
    scale_y_log10(
      labels = log_labels_bold,
      breaks = log_breaks,
      limits = log_limits,
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      fill = "",
      color = "",
      title = panel_title,
      y = y_title,
      x = ""
    )

  gPlot(p)
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
message("180 # BUILD PANEL TESTS")
################################################################################

df_plot_all <- df_plot %>%
  filter(figure %in% c("5A", "5B"))

df_tests_all <- summary_brackets(
  df = df_plot_all,
  axis_column = "method",
  value_column = "value_non_log10",
  facet_columns = c("figure", "sampling_device"),
  testing_col = "method",
  test_type = NULL,
  adjust_method = "BH",
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE
)$tests

################################################################################
message("260 # BUILD PANELS")
################################################################################

df_plot_A <- df_plot %>%
  filter(figure == "5A")

plotA <- build_panel_plot(
  df_panel = df_plot_A,
  df_tests_all = df_tests_all,
  figure_id = "5A",
  panel_title = "A) JPL-Spacecraft Assembly Facility<br><span style=\"color:white\">A) </span>(25 cm<sup>2</sup>)",
  y_title = "log10(16S rRNA gene copies) / 25 cm<sup>2</sup>"
)

df_plot_B <- df_plot %>%
  filter(figure == "5B")

plotC <- build_panel_plot(
  df_panel = df_plot_B,
  df_tests_all = df_tests_all,
  figure_id = "5B",
  panel_title = "B) JPL-Spacecraft Assembly Facility<br><span style=\"color:white\">B) </span>(1000 cm<sup>2</sup>)",
  y_title = "log10(16S rRNA gene copies) / 1000 cm<sup>2</sup>"
)

################################################################################
message("460 # SAVE OUTPUTS")
################################################################################

final_plot <- plotA | plotC

ggsave(paste0(output_plot, script_title, "_dqPCR_", current_date, ".png"), final_plot, width = 12, height = 8)
ggsave(paste0(output_plot, script_title, "A_dqPCR_", current_date, ".png"), plotA, width = 6, height = 8)
ggsave(paste0(output_plot, script_title, "B_dqPCR_", current_date, ".png"), plotC, width = 6, height = 8)
#saveRDS(plotA, paste0(output_plot, script_title, "A_dqPCR_", current_date, ".rds"))
#saveRDS(plotC, paste0(output_plot, script_title, "B_dqPCR_", current_date, ".rds"))

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
