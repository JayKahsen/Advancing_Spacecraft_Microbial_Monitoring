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
Creates the Figure 4 release-and-recovery comparison plots for the PPR dataset.
Builds the point-summary panel, sequence-read insets, combined Figure 4 output,
and range-summary review files, then writes the PNG, PDF, and CSV outputs.
"

#=============================== Plot settings ================================#

options(scipen = 999)

script_title <- "Figure 4"
figure_number <- "4"
experiment_type_order <- c("metal_deposition", "swab_head_retention")
plot_height <- 12
plot_width <- 20

new_labels <- c(
  "cotton" = "Cotton Swab",
  "macrofoam" = "Macrofoam Swab"
)
palette_label <- c(new_labels, palette_label)

#=================================== END SETUP ===================================#

################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

filtered_labels <- palette_label[sampling_device_order]

theme_common <- theme_global(base_size = 11) +
  theme(
    strip.text = ggtext::element_markdown(size = rel(1.2), face = "bold"),
    axis.text.y = ggtext::element_markdown(face = "bold"),
    axis.title.y = ggtext::element_markdown(face = "bold"),
    legend.position = "bottom"
  )

gPlot <- function(p) {
  p <- p +
    theme_common +
    scale_fill_manual(
      breaks = sampling_device_order,
      values = palette_color,
      labels = filtered_labels,
      na.translate = FALSE
    ) +
    scale_color_manual(values = palette_color, labels = palette_label) +
    scale_linetype_manual(values = palette_linetype, labels = palette_label) +
    scale_shape_manual(values = palette_shape, labels = palette_label)

  print(p)
  p
}

plot_reads_by_figure <- function(df, figure_id) {
  df_plot_data <- df %>%
    filter(str_detect(Figure, figure_id)) %>%
    mutate(sampling_device = factor(sampling_device, levels = sampling_device_order))

  p <- ggplot(df_plot_data, aes(x = sampling_device, y = sample_reads, fill = sampling_device)) +
    facet_grid2(~figure, scales = "free_x", axes = "all", remove_labels = "none") +
    scale_y_log10(labels = log_labels_bold) +
    stat_summary(fun = mean, geom = "col", color = "black", width = 0.5) +
    geom_jitter(width = 0.1, alpha = 0.7, size = 2, shape = 21) +
    labs(
      y = "Mean Sequence Reads",
      x = "",
      color = "",
      fill = ""
    )

  p_out <- gPlot(p) +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_blank()
    ) +
    coord_cartesian(ylim = c(10^4.5, 10^7.65)) +
    ggtitle("")

  p_out
}

################################################################################
message("90 # MAIN")
################################################################################

#.............................. prepare input data .............................#

df <- meta %>%
  filter(experiment_type %in% experiment_type_order) %>%
  filter(!is.na(percent_recovery_dpcr)) %>%
  filter(sampling_device %in% sampling_device_order) %>%
  select(
    sample_id,
    percent_recovery_qpcr,
    percent_recovery_dpcr,
    sampling_device,
    experiment_type,
    figure,
    Figure,
    bar_fill,
    sample_reads
  )

figure_4_selection <- df %>%
  mutate(
    figure_4 = case_when(
      experiment_type == "swab_head_retention" ~ "4A",
      experiment_type == "metal_deposition" & sampling_device %in% c("macrofoam", "cotton") ~ "4B",
      sampling_device %in% c("SALSA", "Polyester Wipe") ~ "4C",
      TRUE ~ NA_character_
    )
  ) %>%
  select(sample_id, figure_4)

write.csv(figure_4_selection, "data_tables/figure_4_selection.csv", row.names = FALSE)

df_plot <- df %>%
  filter(str_detect(Figure, figure_number)) %>%
  pivot_longer(
    cols = c(percent_recovery_dpcr, percent_recovery_qpcr),
    names_to = "method",
    values_to = "value"
  ) %>%
  mutate(method = str_replace(method, "percent_recovery_", "")) %>%
  mutate(method = factor(method, levels = method_order))

################################################################################
message("180 # BUILD SUMMARY OBJECTS")
################################################################################

results <- summary_brackets(
  df = df_plot,
  axis_column = "sampling_device",
  value_column = "value",
  facet_columns = c("figure", "method"),
  testing_col = "sampling_device",
  test_type = NULL,
  adjust_method = "BH",
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE
)

df_summary <- results$summary %>%
  filter(figure %in% c("4A", "4B", "4C")) %>%
  left_join(meta %>% distinct(bar_fill, sampling_device)) %>%
  mutate(method = factor(method, levels = method_order))

df_segments <- results$segments %>%
  filter(figure %in% c("4A", "4B", "4C")) %>%
  mutate(method = factor(method, levels = method_order))

df_tests <- results$tests %>%
  filter(figure %in% c("4A", "4B", "4C")) %>%
  mutate(
    p.value = welch_raw_val,
    significance = welch_raw_sig,
    method = factor(method, levels = method_order)
  )

################################################################################
message("250 # BUILD POINT SUMMARY")
################################################################################

p <- ggplot(df_plot, aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_jitter(width = 0.2, size = point_size + 1, alpha = 0, shape = 21) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  facet_grid(
    method ~ figure,
    scale = "free_x",
    labeller = labeller(.cols = palette_label, .rows = palette_label)
  ) +
  geom_text(
    data = df_tests,
    aes(x = label_position, y = label_value, label = significance),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = -0.5
  ) +
  geom_text(
    data = df_tests,
    aes(x = label_position, y = label_value, label = scales::scientific(p.value)),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = 1.5
  ) +
  geom_errorbar(
    data = df_summary,
    aes(x = sampling_device, ymin = mean_value, ymax = mean_value, color = sampling_device),
    width = 0.6,
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  geom_jitter(width = 0.2, size = point_size + 1, alpha = 0.8, shape = 21) +
  geom_errorbar(
    data = df_summary,
    aes(x = sampling_device, ymin = mean_minus_sd, ymax = mean_plus_sd),
    width = 0.1,
    linewidth = 1,
    inherit.aes = FALSE,
    color = bar_color
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  scale_x_discrete(labels = palette_label) +
  labs(
    fill = "Sampling Device",
    color = "",
    title = "",
    y = "Percent Recovery (%)",
    x = ""
  ) +
  guides(color = "none")

plot_A_13 <- gPlot(p) +
  theme(legend.position = "right")

################################################################################
message("340 # BUILD READS PANELS")
################################################################################

df_reads <- df %>%
  filter(str_detect(Figure, figure_number)) %>%
  mutate(sampling_device = factor(sampling_device, levels = sampling_device_order))

reads_plot_4A <- plot_reads_by_figure(df = df_reads, figure_id = "4A")
reads_plot_4B <- plot_reads_by_figure(df = df_reads, figure_id = "4B")
reads_plot_4C <- plot_reads_by_figure(df = df_reads, figure_id = "4C")

################################################################################
message("420 # ASSEMBLE FINAL OUTPUT")
################################################################################

st_plot <- readRDS(paste0(output_plot, "Figure_", figure_number, "_st_plot.rds")) +
  scale_x_log10_shift(
    shift = 2,
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = log10_labels_percent,
    limits = c(0.01, 100)
  )

theme_custom <- theme(
  axis.text.y = ggtext::element_markdown(size = rel(1.2), face = "bold"),
  legend.text = element_text(size = rel(1.5)),
  legend.title = ggtext::element_markdown(size = rel(1.5)),
  strip.text = ggtext::element_markdown(size = rel(1.4), face = "bold"),
  axis.title.x = ggtext::element_markdown(size = rel(1.2), face = "bold", color = "black"),
  axis.text = ggtext::element_markdown(size = rel(1.2), face = "bold")
)

plot_A_13 <- plot_A_13 +
  theme_custom +
  theme(legend.box.margin = margin(l = -20)) +
  ggtitle("")

st_plot <- st_plot +
  theme_custom +
  theme(
    legend.position = "none",
    plot.margin = margin(l = 54.5, r = 43.5)
  )

row1 <- wrap_elements(plot_A_13)
row2 <- wrap_elements(st_plot)

combined1 <- row1 / row2 +
  plot_layout(heights = c(10, 12))

top1 <- 0.55
bottom1 <- 0.088
width1 <- 0.08
width2 <- 0.274
left2_adj <- 0.00

left1 <- 0.208
left2 <- left1 + width2 - left2_adj
left3 <- left2 + width2 + left2_adj

combined <- combined1 +
  inset_element(reads_plot_4A, left = left1, bottom = bottom1, right = left1 + width1, top = top1) +
  inset_element(reads_plot_4B, left = left2, bottom = bottom1, right = left2 + width1, top = top1) +
  inset_element(reads_plot_4C, left = left3, bottom = bottom1, right = left3 + width1, top = top1)

ggsave(
  paste0(output_plot, "Figure_", figure_number, "_Release_Recovery_", current_date, ".png"),
  combined,
  width = plot_width,
  height = plot_height,
  dpi = dpi
)

ggsave(
  paste0(output_plot, "Figure_", figure_number, "_Release_Recovery_BH_corrected_", current_date, ".pdf"),
  combined,
  width = plot_width,
  height = plot_height,
  dpi = dpi
)

################################################################################
message("520 # SAVE RANGE SUMMARY")
################################################################################

range_summary <- df_plot %>%
  group_by(method, figure, sampling_device) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(range_summary, paste0(output_plot, "fig_4_ranges.csv"), row.names = FALSE)

p_range <- p +
  geom_text(
    data = range_summary,
    aes(
      x = sampling_device,
      y = 50,
      label = paste0(
        "min=", round(min_value, 2),
        "\nmax=", round(max_value, 2),
        "\nmean=", round(mean_value, 2)
      )
    ),
    size = 3,
    color = "black",
    lineheight = 0.9,
    vjust = 0
  )

ggsave(paste0(output_plot, "Figure_", figure_number, "_Release_Recovery_ranges.png"), p_range, width = 12, height = 6)

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
