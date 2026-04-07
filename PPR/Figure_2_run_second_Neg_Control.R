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
Creates the Figure 2 negative-control comparison plots for the PPR dataset.
Builds the range summary plot, the sequence-read bar plot, and the combined
negative-control figure, then writes review outputs to output_plot.
"

#=============================== Plot settings ================================#

figure_number <- "2"
options(scipen = 999)

script_title <- "Figure 2"

#=================================== END SETUP ===================================#

################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

theme_common <- theme_global(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom"
  )

gPlot <- function(p) {
  p <- p +
    scale_color_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
    scale_size_manual(values = palette_size, labels = palette_label) +
    labs(fill = "", size = "Sequence Reads") +
    theme_common

  print(p)
  p
}

################################################################################
message("90 # MAIN")
################################################################################

#.............................. prepare input data .............................#

df <- meta %>%
  filter(str_detect(Figure, figure_number)) %>%
  select(sample_id, Figure, bar_name, copies_rxn_dpcr, copies_rxn_qpcr, bar_color, bar_type_2, bar_axis_2, sample_reads)

figure_2_selection <- df %>%
  mutate(figure_2 = "2") %>%
  select(sample_id, figure_2)

write.csv(figure_2_selection, "data_tables/figure_2_selection.csv", row.names = FALSE)

sample_read_size_order <- c("ng100", "g100", "g1K", "g10K", "g100K")

df_plot <- df %>%
  mutate(sample_read_size = case_when(
    is.na(sample_reads) ~ "ng100",
    sample_reads > 100000 ~ "g100K",
    sample_reads > 10000 ~ "g10K",
    sample_reads > 1000 ~ "g1K",
    sample_reads > 100 ~ "g100",
    TRUE ~ "ng100"
  )) %>%
  mutate(
    copies_rxn_dpcr = as.numeric(copies_rxn_dpcr),
    copies_rxn_qpcr = as.numeric(copies_rxn_qpcr)
  ) %>%
  pivot_longer(
    cols = c(copies_rxn_dpcr, copies_rxn_qpcr),
    names_to = "method",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(method = str_replace(method, "copies_rxn_", "")) %>%
  mutate(bar_type_2 = factor(bar_type_2, levels = bar_type_2_order)) %>%
  mutate(bar_name = factor(bar_name, levels = bar_name_order)) %>%
  mutate(bar_axis_2 = factor(bar_axis_2, levels = bar_axis_2_order)) %>%
  mutate(bar_color = factor(bar_color, levels = bar_color_order)) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  mutate(sample_read_size = factor(sample_read_size, levels = sample_read_size_order))

df_summary <- df_plot %>%
  group_by(bar_name, bar_type_2, bar_color, bar_axis_2, method) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(bar_type_2 = factor(bar_type_2, levels = bar_type_2_order)) %>%
  mutate(bar_name = factor(bar_name, levels = bar_name_order)) %>%
  mutate(bar_axis_2 = factor(bar_axis_2, levels = bar_axis_2_order)) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  mutate(bar_color = factor(bar_color, levels = bar_color_order))

################################################################################
message("180 # BUILD RANGE PLOT")
################################################################################

filtered_labels <- palette_label[bar_color]

p <- ggplot(df_plot, aes(x = bar_axis_2)) +
  geom_point(aes(x = bar_axis_2, y = value, color = bar_color), inherit.aes = FALSE, shape = 15, alpha = 0) +
  geom_errorbar(data = df_summary, aes(ymin = mean_value, ymax = mean_value, color = bar_color), width = 0.6, linewidth = 1) +
  geom_point(
    aes(x = bar_axis_2, y = value, fill = bar_color, size = sample_read_size),
    position = position_jitter(width = 0.1),
    inherit.aes = FALSE,
    color = "gray18",
    shape = 21,
    alpha = 0.8
  ) +
  geom_errorbar(data = df_summary, aes(ymin = mean_value - sd, ymax = mean_value + sd), width = 0.2, color = "gray18", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_x_discrete(labels = palette_label) +
  scale_fill_manual(
    breaks = bar_color,
    values = palette_color,
    labels = filtered_labels,
    na.translate = FALSE
  ) +
  scale_color_manual(
    values = palette_color,
    labels = palette_label,
    na.translate = FALSE
  ) +
  scale_size_manual(
    values = palette_size,
    labels = palette_label,
    na.translate = FALSE
  ) +
  facet_grid(
    method ~ bar_type_2,
    scale = "free_x",
    labeller = labeller(.cols = palette_label, .rows = palette_label)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(shape = 15, linetype = NA, alpha = 1, size = 8)
    ),
    fill = "none"
  ) +
  labs(
    title = "Recovery of 16S rRNA Gene Copies from Negative Controls",
    color = "Control Type",
    fill = "",
    size = "# of Sequence Reads",
    y = "Log10 (16S dPCR copies) / reaction",
    x = "Negative Controls"
  )

p1 <- p

range_summary <- df_plot %>%
  group_by(bar_axis_2, method, bar_type_2) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

p_range <- p1 +
  geom_text(
    data = range_summary,
    aes(
      x = bar_axis_2,
      y = 100,
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

write.csv(range_summary, paste0(output_plot, "Figure_", figure_number, "_ranges_", current_date, ".csv"), row.names = FALSE)
ggsave(paste0(output_plot, "Figure_", figure_number, "_ranges_", current_date, ".png"), p_range, width = 12, height = 12, dpi = dpi)

################################################################################
message("300 # BUILD SAMPLE BAR PLOT")
################################################################################

df_sample_counts <- df_plot %>%
  select(sample_id, bar_color, bar_type_2, sample_reads) %>%
  distinct()

sample_bar1 <- ggplot(df_sample_counts, aes(x = bar_color, y = sample_reads)) +
  facet_grid2(~bar_type_2, scale = "free") +
  stat_summary(
    aes(fill = bar_color),
    fun = mean,
    geom = "bar"
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.2
  ) +
  scale_y_log10(labels = log_labels_bold) +
  labs(x = "", y = "", title = "Sequence Reads") +
  scale_color_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
  scale_fill_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
  theme_global(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_markdown(),
    legend.position = "bottom"
  )

ggsave(paste0(output_plot, "Figure-", figure_number, "_sample_bar.png"), sample_bar1, width = 12)

################################################################################
message("420 # ASSEMBLE FINAL OUTPUT")
################################################################################

theme_custom <- theme(
  legend.key.size = unit(0.9, "cm"),
  plot.title = element_markdown(size = rel(1.6), face = "bold", hjust = 0),
  axis.text.y = element_markdown(size = rel(1.2), face = "bold", margin = margin(r = 10)),
  strip.text = element_text(size = rel(1.35), face = "bold"),
  legend.text = element_text(size = rel(1.35)),
  legend.title = element_markdown(size = rel(1.35), face = "bold"),
  axis.text = element_markdown(size = rel(1.2), face = "bold")
)

st_plot1 <- readRDS(paste0(output_plot, "Figure-", figure_number, "_st_plot.rds"))

st_plot <- st_plot1 +
  theme_custom +
  theme(plot.margin = margin(t = 10, l = 60)) +
  guides(color = "none", fill = "none") +
  ggtitle("B) Percent Relative Abundance") +
  scale_x_log10(breaks = c(1, 10, 100), labels = log10_labels_percent, limits = c(1, 100))

p <- gPlot(p1) +
  theme_custom +
  ggtitle("A) Recovery of 16S rRNA Gene Copies and Sequence Reads from Negative Controls") +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(b = 0),
    legend.margin = margin(t = 20, r = 90, b = 0, l = 25)
  )

sample_bar <- sample_bar1 +
  theme_custom +
  labs(y = "Sequence Reads") +
  theme(
    plot.title = element_blank(),
    strip.text = element_blank(),
    plot.margin = margin(t = 0),
    legend.position = "none"
  )

row1 <- (p / sample_bar) +
  plot_layout(heights = c(2, 1), guides = "collect")

row2 <- wrap_elements(st_plot)
row1 <- wrap_elements(row1)

combined <- (row1 / row2) +
  plot_layout(heights = c(5, 5))

ggsave(paste0(output_plot, "Figure_", figure_number, "_Neg_Control_", current_date, ".pdf"), combined, width = 16, height = 16, dpi = dpi)

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
