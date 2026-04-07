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
Creates the Figure 3 limit-of-detection plots for the PPR dataset.
Builds the Figure 3 LOD panel, the regression review plot, and the final
three-panel output, then writes the PNG and PDF files to output_plot.
"

#=============================== Plot settings ================================#

options(scipen = 999)

script_title <- "Figure_3"
figure_number <- "3"
r <- 3
y_adj <- 10
plot_height <- 10
plot_width <- 20

#=================================== END SETUP ===================================#

################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

theme_common <- theme_global(base_size = 11) +
  theme(
    axis.text.y = ggtext::element_markdown(face = "bold"),
    axis.title.y = ggtext::element_markdown(face = "bold"),
    legend.position = "none"
  )

gPlot <- function(p) {
  p <- p +
    theme_common +
    scale_linetype_manual(values = palette_linetype, labels = palette_label) +
    scale_shape_manual(values = palette_shape, labels = palette_label) +
    labs(caption = "", fill = "", color = "") +
    scale_x_discrete(labels = palette_label) +
    guides(fill = "none", color = "none") +
    theme_common

  print(p)
  p
}

################################################################################
message("90 # MAIN")
################################################################################

#.............................. prepare input data .............................#

df_matrix1 <- read.csv(matrix_names[r, "file_path"], check.names = FALSE, row.names = 1) %>%
  t() %>%
  as.data.frame()

#................................. build metadata ..............................#

df <- meta %>%
  select(
    unique_replicate_id,
    sample_name,
    sample_id,
    sample_type,
    copies_rxn_dpcr,
    copies_rxn_dpcr_log10,
    copies_rxn_qpcr,
    copies_rxn_qpcr_log10,
    extraction_method,
    experiment_type,
    figure,
    Figure
  ) %>%
  filter(str_detect(Figure, figure_number))

df_bar_plot1 <- df %>%
  mutate(copies_rxn_dpcr = as.numeric(copies_rxn_dpcr)) %>%
  mutate(copies_rxn_qpcr = as.numeric(copies_rxn_qpcr)) %>%
  pivot_longer(
    cols = c(copies_rxn_dpcr, copies_rxn_qpcr),
    names_to = "method",
    values_to = "value_non_log10"
  ) %>%
  mutate(method = str_replace(method, "copies_rxn_", "")) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = method_order)) %>%
  select(sample_name, method, value_non_log10)

df_bar_plot <- df %>%
  pivot_longer(
    cols = c(copies_rxn_dpcr_log10, copies_rxn_qpcr_log10),
    names_to = "method",
    values_to = "value"
  ) %>%
  group_by(method) %>%
  mutate(
    LOD = max(
      value[
        grepl("control", sample_type, ignore.case = TRUE) &
          experiment_type == "limit_of_detection"
      ],
      value[
        grepl("NTC", sample_type, ignore.case = TRUE) &
          figure == "3B"
      ],
      na.rm = TRUE
    )
  ) %>%
  mutate(method = str_replace(method, "copies_rxn_", "")) %>%
  mutate(method = str_replace(method, "_log10", "")) %>%
  left_join(df_bar_plot1, by = c("sample_name", "method")) %>%
  filter(!is.na(value)) %>%
  select(experiment_type, method, sample_type, value, LOD, everything()) %>%
  filter(!sample_type %in% c("NTC")) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = method_order))

################################################################################
#.............................. filter figure 3A ..............................#
################################################################################

figure_3_selection <- df %>%
  mutate(figure_3 = "3") %>%
  mutate(figure_3 = case_when(
    experiment_type == "limit_of_detection" ~ "3B",
    experiment_type == "extract_efficiency" ~ "3A",
    TRUE ~ "3B"
  )) %>%
  select(sample_id, figure_3)

write.csv(figure_3_selection, "data_tables/figure_3_selection.csv", row.names = FALSE)

################################################################################
#............................ filter matrix figure 3B ............................#
################################################################################

df_matrix <- df_matrix1 %>%
  xPlode_sample_name() %>%
  filter(str_detect(figure, "3B")) %>%
  filter(!str_detect(sample_name, "NTC")) %>%
  filter(str_detect(Figure, figure_number)) %>%
  column_to_rownames(var = "sample_name") %>%
  select(-any_of(meta_names)) %>%
  mutate_all(as.numeric) %>%
  select(which(colSums(.) > 0))

feature_names <- names(df_matrix)

df_long <- df_matrix %>%
  xPlode_sample_name() %>%
  pivot_longer(cols = any_of(feature_names), names_to = "feature", values_to = "counts") %>%
  mutate(old_feature = feature) %>%
  select(sample_name, feature, old_feature, counts, everything()) %>%
  mutate(feature = case_when(
    str_detect(feature, atcc_genera[1]) ~ atcc_genera[1],
    str_detect(feature, atcc_genera[2]) ~ atcc_genera[2],
    str_detect(feature, atcc_genera[3]) ~ atcc_genera[3],
    str_detect(feature, atcc_genera[4]) ~ atcc_genera[4],
    str_detect(feature, atcc_genera[5]) ~ atcc_genera[5],
    str_detect(feature, atcc_genera[6]) ~ atcc_genera[6],
    str_detect(feature, atcc_genera[7]) ~ atcc_genera[7],
    str_detect(feature, atcc_genera[8]) ~ atcc_genera[8],
    str_detect(feature, atcc_genera[9]) ~ atcc_genera[9],
    str_detect(feature, atcc_genera[10]) ~ atcc_genera[10],
    TRUE ~ "other"
  ))

df_percent_other <- df_long %>%
  group_by(sample_name) %>%
  mutate(sample_count = sum(counts)) %>%
  ungroup() %>%
  filter(feature == "other") %>%
  group_by(sample_name) %>%
  mutate(other_count = sum(counts)) %>%
  ungroup() %>%
  mutate(percent_other = 100 * other_count / sample_count) %>%
  distinct(sample_name, percent_other)

df_long_atcc <- df_long %>%
  filter(feature != "other") %>%
  group_by(sample_name) %>%
  mutate(sample_counts = sum(counts)) %>%
  group_by(sample_name, feature) %>%
  mutate(feature_count = sum(counts)) %>%
  ungroup() %>%
  mutate(percent_feature = 100 * feature_count / sample_counts) %>%
  ungroup() %>%
  left_join(df_percent_other) %>%
  mutate(relative_detection = percent_feature / 10) %>%
  mutate(log2_value = log2(relative_detection)) %>%
  mutate(feature_old = feature) %>%
  mutate(feature = feature_old) %>%
  mutate(feature = ifelse(feature_old == "other", "Other", paste0("*", feature_old, "*"))) %>%
  select(rev(everything())) %>%
  select(sample_name, feature, percent_feature, percent_other, relative_detection, log2_value, everything())

################################################################################
message("250 # SUMMARY BRACKETS")
################################################################################

plot_title <- "dilutions bar 3B_1"
dodge_offset <- 0.2
position_adj <- 4

df_plot2 <- df_bar_plot %>%
  filter(figure %in% c("3B", "2_3B")) %>%
  mutate(extraction_method = factor(extraction_method, levels = extraction_method_order)) %>%
  arrange(extraction_method) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  mutate(bar_fill = ifelse(as.integer(extraction_method) %% 2 == 1, "alt_group1", "alt_group2")) %>%
  distinct()

results <- summary_brackets(
  df = df_plot2,
  axis_column = "sample_type",
  dodge_column = "extraction_method",
  value_column = "value_non_log10",
  facet_columns = "method",
  testing_col = "extraction_method",
  test_type = NULL,
  adjust_method = "BH",
  dodge_width = 0.8,
  width = 0.8,
  log = TRUE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE
)

df_tests1 <- results$tests %>%
  filter(sample_type != "NTC") %>%
  mutate(method = factor(method, levels = method_order)) %>%
  mutate(adj.significance = welch_log_BH_sig) %>%
  mutate(adj.p.value = welch_log_BH_val) %>%
  ungroup()

results <- summary_brackets(
  df = df_plot2,
  axis_column = "sample_type",
  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = "method",
  testing_col = "extraction_method",
  test_type = NULL,
  adjust_method = "BH",
  dodge_width = 0.8,
  width = 0.8,
  log = FALSE,
  space_adj = 1.1,
  tip_adj = 0.8,
  run_tests = TRUE
)

df_summary <- results$summary %>%
  mutate(extraction_method = factor(extraction_method, levels = extraction_method_order)) %>%
  arrange(extraction_method) %>%
  mutate(extraction_method = factor(extraction_method, levels = unique(extraction_method))) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  mutate(bar_fill = ifelse(as.integer(extraction_method) %% 2 == 1, "alt_group1", "alt_group2"))

df_segments <- results$segments %>%
  mutate(method = factor(method, levels = method_order)) %>%
  filter(sample_type != "NTC") %>%
  ungroup()

df_tests <- df_tests1 %>%
  left_join(
    df_segments %>%
      filter(part == "top") %>%
      mutate(method = factor(method, levels = method_order)) %>%
      select(sample_type, method, y)
  ) %>%
  mutate(adj.significance = ifelse(is.na(adj.significance), "na", adj.significance))

################################################################################
message("340 # BUILD LOD PANEL")
################################################################################

breaks <- as.character(unique(df_plot2$extraction_method))
break_labels <- palette_label[breaks]

p <- ggplot(df_plot2, aes(x = sample_type, y = value, fill = extraction_method, group = extraction_method)) +
  theme(
    legend.key.height = unit(1.5, "cm"),
    legend.text.align = 0.5
  ) +
  geom_hline(aes(yintercept = LOD), linetype = "dashed", color = "red") +
  geom_text(aes(x = "2.00E+5", y = LOD), label = "LOD", hjust = 0.5, vjust = -0.5, color = "red") +
  facet_grid(
    method ~ .,
    labeller = labeller(.cols = palette_label, .rows = palette_label)
  ) +
  geom_errorbar(
    data = df_summary,
    aes(
      x = sample_type,
      ymin = mean_value,
      ymax = mean_value,
      group = extraction_method,
      color = extraction_method
    ),
    position = position_dodge(width = 0.8),
    width = 0.6,
    size = 1.5,
    inherit.aes = FALSE
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      jitter.height = 0,
      dodge.width = 0.8
    ),
    size = point_size,
    alpha = 0.8,
    shape = 21
  ) +
  geom_errorbar(
    data = df_summary,
    aes(
      x = sample_type,
      ymin = mean_minus_sd,
      ymax = mean_plus_sd,
      group = extraction_method
    ),
    position = position_dodge(width = 0.8),
    width = 0.1,
    size = 1,
    inherit.aes = FALSE,
    color = bar_color
  ) +
  coord_cartesian(ylim = c(0.1, NA)) +
  scale_y_continuous(
    labels = log10_labels_bold,
    expand = expansion(mult = c(0, 0.1))
  ) +
  coord_cartesian(ylim = zoom) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = y,
      label = adj.significance
    ),
    inherit.aes = FALSE,
    size = 4,
    vjust = -0.5
  ) +
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = y,
      label = scales::scientific(adj.p.value)
    ),
    inherit.aes = FALSE,
    size = 4,
    vjust = 2.2
  ) +
  labs(
    fill = "",
    color = "",
    title = "Determining the Limit of Detection for Two<br>Distinct DNA Extraction Methods",
    y = "log10(16S rRNA gene copies) / reaction",
    x = "Concentration of MMC"
  ) +
  scale_fill_manual(
    breaks = breaks,
    values = palette_color,
    labels = break_labels,
    na.translate = FALSE
  ) +
  scale_color_manual(
    breaks = breaks,
    values = palette_color,
    labels = break_labels,
    na.translate = FALSE
  ) +
  theme_common +
  scale_x_discrete(labels = palette_label)

plot_3B_1 <- gPlot(p) + theme(legend.position = "none")

################################################################################
message("450 # REGRESSION REVIEW")
################################################################################

plot_regression <- plot_3B_1 +
  geom_smooth(
    data = df_plot2 %>% filter(sample_type != "reagent_control"),
    aes(color = extraction_method),
    method = "lm",
    se = FALSE,
    formula = y ~ x
  ) +
  stat_poly_eq(
    data = df_plot2 %>% filter(sample_type != "reagent_control"),
    aes(
      color = extraction_method,
      label = after_stat(paste(eq.label, rr.label, p.value.label, sep = "~~~"))
    ),
    formula = y ~ x,
    parse = TRUE
  )

ggsave(paste0(output_plot, "Figure_3_regression_", current_date, ".png"), plot_regression, width = 12, height = 12, dpi = dpi)

################################################################################
message("500 # ASSEMBLE FINAL OUTPUTS")
################################################################################

plot_title <- script_title

theme_custom <- theme(
  legend.text = element_markdown(size = rel(1.5)),
  strip.text = element_text(size = rel(1.35), face = "bold"),
  axis.title.x = element_markdown(size = rel(1.2), face = "bold", color = "black"),
  axis.text = element_markdown(size = rel(1.2), face = "bold")
)

rel_log_plot1 <- readRDS(paste0(output_plot, "Figure-", figure_number, "_rel_log_plot.rds"))
counts_log_plot1 <- readRDS(paste0(output_plot, "Figure-", figure_number, "_counts_log_plot.rds"))

df_tag1 <- data.frame(
  method = "dpcr",
  extraction_method = "mag_beads",
  label = "A",
  x = -Inf,
  y = Inf
) %>%
  mutate(method = factor(method, levels = method_order))

tag_size <- 6
x_adj <- -0.4
y_adj <- 1 - x_adj

df_tag <- df_tag1

plot3B1 <- plot_3B_1 +
  labs(
    title = "",
    y = "A) log10(16S rRNA gene copies) / reaction",
    x = "Input level of MMC (cells)"
  ) +
  geom_text(
    data = df_tag,
    aes(x = x, y = y, label = label),
    hjust = x_adj,
    vjust = y_adj,
    fontface = "bold",
    size = tag_size,
    inherit.aes = FALSE
  ) +
  theme_custom +
  theme(
    legend.position = c(0.5, 0.55),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing = unit(0, "pt"),
    legend.spacing.x = unit(0, "pt"),
    legend.spacing.y = unit(0, "pt"),
    plot.margin = margin(0, 0, 0, 0),
    axis.text.x = element_markdown(angle = 45, hjust = 1)
  )

df_tag <- df_tag1 %>%
  mutate(label = "B")

counts_log_plot <- counts_log_plot1 +
  labs(
    y = "B) Sequencing depth (reads)",
    x = "Input level of MMC (cells)"
  ) +
  theme_custom +
  scale_y_log10(labels = log_labels_bold) +
  geom_text(
    data = df_tag,
    aes(x = x, y = y, label = label),
    hjust = x_adj,
    vjust = y_adj,
    fontface = "bold",
    size = tag_size,
    inherit.aes = FALSE
  ) +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    legend.position = "none"
  )

df_tag <- df_tag1 %>%
  mutate(label = "C")

rel_log_plot <- rel_log_plot1 +
  labs(
    y = "C) Relative abundance (%)",
    x = "Input level of MMC (cells)"
  ) +
  geom_text(
    data = df_tag,
    aes(x = x, y = y, label = label),
    hjust = x_adj,
    vjust = y_adj,
    fontface = "bold",
    size = tag_size,
    inherit.aes = FALSE
  ) +
  theme_custom +
  theme(
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.margin = margin(l = -20)
  )

final_plot <- plot3B1 | counts_log_plot | rel_log_plot

ggsave(paste0(output_plot, plot_title, "_LOD_", current_date, ".png"), final_plot, width = plot_width, height = plot_height, dpi = dpi)
ggsave(paste0(output_plot, plot_title, "_LOD_BH_corrected_", current_date, ".pdf"), final_plot, width = plot_width, height = plot_height, dpi = dpi)

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
