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
Creates the Figure 3 community-composition review plots for the PPR dataset.
Builds the relative-abundance and sequencing-depth ATCC mock-community plots,
then writes the review PNG and RDS outputs to output_plot.
"

#=============================== Plot settings ================================#

script_title <- "abundance"
figure_number <- "3"
zoom <- 3
starting_r <- 1
ending_r <- 1
testing <- "yes"
testing_r <- 3
selected_data_class <- "doubleton"
selected_data_class_name <- paste(selected_data_class, collapse = "_")
use_custom_labels <- "no"
extraction_method_order <- rev(c("mag_beads", "spin_column"))
figure_order <- c("3", "3A", "3B")

#============================= Parameter settings =============================#

parameter_sets <- list(
  set1 = list(
    filter1_group = "figure",
    y_axis_group = "extraction_method",
    x_facet_group = "sample_type",
    y_facet_group = "experiment",
    number_size = 3
  )
)

#=================================== END SETUP ===================================#

################################################################################
message("90 # MAIN")
################################################################################

Loop_Group_order <- Loop_Group <- "run_once"
if (exists("Run_Group")) {
  Loop_Group <- Run_Group
  Loop_Group_order <- Run_Group_order
}

if (!exists("ending_r")) {
  ending_r <- nrow(matrix_names)
}

if (!exists("starting_r")) {
  starting_r <- 1
}

if (testing == "yes") {
  starting_r <- testing_r
  ending_r <- testing_r
}

for (lp in Loop_Group_order) {
  for (r in starting_r:ending_r) {
    taxa_levs <- matrix_names[r, "taxa_levs"]
    data_set <- matrix_names[r, "data_sets"]
    taxa_plural <- matrix_names[r, "taxa_plural"]

    params <- parameter_sets[[paste0("set", match(data_set, data_set_order))]]
    list2env(params, envir = environment())

    groups <- c("filter1_group", "y_axis_group", "x_facet_group", "y_facet_group")
    for (var in groups) {
      group_value <- get(var)
      assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
    }

    grouping_columns <- c(y_axis_group, x_facet_group, y_facet_group)
    run_suffix <- ""
    if (exists("Run_Group")) {
      run_suffix <- paste0("_", lp)
    }

    p_title <- paste0(taxa_levs, "_", script_title, custom_name, "_", data_set, run_suffix)
    qPrint(p_title)

    ################################################################################
    message("124 # LOAD AND FILTER DATA")
    ################################################################################

    df_matrix1 <- read.csv(matrix_names[r, "file_path"], check.names = FALSE, row.names = 1) %>%
      t() %>%
      as.data.frame() %>%
      xPlode_sample_name() %>%
      filter(.data[[filter1_group]] %in% filter1_group_order) %>%
      filter(.data[[y_facet_group]] %in% y_facet_group_order) %>%
      filter(.data[[x_facet_group]] %in% x_facet_group_order) %>%
      filter(.data[[y_axis_group]] %in% y_axis_group_order) %>%
      filter(!sample_type %in% c("reagent_control", "NTC", "positive_control")) %>%
      filter(str_detect(Figure, figure_number))

    if (exists("Run_Group")) {
      df_matrix1 <- df_matrix1 %>%
        filter(.data[[Loop_Group]] %in% lp)
    }

    df_matrix2 <- df_matrix1 %>%
      imPlode_sample_name() %>%
      mutate_all(as.numeric) %>%
      select(which(colSums(.) > 0))

    df_matrix <- combine_other(df = df_matrix2, patterns = atcc_genera, starts_with = TRUE)
    feature_names <- names(df_matrix)

    ################################################################################
    message("170 # DATA PROCESSING")
    ################################################################################

    feature_labels_df <- feature_labels(y_facet_group, averaging_group = x_facet_group)
    y_axis_group_order <- rev(y_axis_group_order)

    plot_ind_df <- df_matrix %>%
      xPlode_sample_name() %>%
      pivot_longer(cols = any_of(feature_names), names_to = "feature", values_to = "counts") %>%
      left_join(feature_labels_df %>% select(-c(taxa)), by = c("feature", y_facet_group)) %>%
      mutate(
        label_new = ifelse(
          feature_label %in% atcc_genera,
          paste0("<i>", feature_label, "</i> "),
          feature_label
        )
      ) %>%
      group_by(feature, across(any_of(y_facet_group))) %>%
      mutate(grouped_feature_counts = sum(counts)) %>%
      filter(grouped_feature_counts > 0) %>%
      group_by(sample_name) %>%
      mutate(sample_counts = sum(counts)) %>%
      ungroup() %>%
      mutate(rel_abun = counts / sample_counts) %>%
      group_by(across(any_of(y_facet_group))) %>%
      mutate(group_counts = mean(sample_counts)) %>%
      ungroup() %>%
      mutate(counts_log10 = log10(counts + 1e-9)) %>%
      mutate(rel_abun_log10 = log10(rel_abun + 1e-9)) %>%
      mutate(group_counts_log10 = log10(group_counts + 1e-9)) %>%
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
      mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_facet_group_order)) %>%
      mutate(!!sym(y_axis_group) := factor(!!sym(y_axis_group), levels = y_axis_group_order)) %>%
      mutate(!!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order)) %>%
      ungroup() %>%
      arrange(plot_order) %>%
      filter(label_new != "Other") %>%
      distinct()

    label_new_order <- plot_ind_df %>%
      select(sample_type, extraction_method, rel_abun, label_new) %>%
      filter(sample_type == "2.00E+6") %>%
      filter(extraction_method == "spin_column") %>%
      group_by(label_new) %>%
      summarize(mean_abun = mean(rel_abun), .groups = "drop") %>%
      arrange(desc(mean_abun)) %>%
      pull(label_new) %>%
      unique()

    plot_ind_df <- plot_ind_df %>%
      mutate(label_new = factor(label_new, levels = label_new_order))

    plot_df <- plot_ind_df %>%
      group_by(feature, across(any_of(grouping_columns))) %>%
      summarize(
        group_mean_counts = mean(counts),
        sample_counts = mean(sample_counts),
        mean_rel_abun = mean(rel_abun),
        percent_rel_abun = 100 * mean(rel_abun),
        sd_rel_abun = sd(rel_abun, na.rm = TRUE),
        group_counts = mean(group_counts),
        mean_rel_log10 = mean(rel_abun_log10),
        sd_rel_log10 = sd(rel_abun_log10),
        mean_counts_log10 = mean(counts_log10),
        sd_counts_log10 = sd(counts_log10),
        .groups = "drop"
      ) %>%
      mutate(xmax_log10 = log10(group_mean_counts + sd_rel_abun + 1e-9)) %>%
      mutate(mean_rel_abun_log10 = log10(mean_rel_abun + 1e-9)) %>%
      mutate(group_counts_log10 = log10(group_counts + 1e-9)) %>%
      mutate(sd_rel_abun_log10 = log10(sd_rel_abun + 1e-9)) %>%
      left_join(feature_labels_df %>% select(-c(taxa)), by = c("feature", y_facet_group)) %>%
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
      mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_facet_group_order)) %>%
      mutate(!!sym(y_axis_group) := factor(!!sym(y_axis_group), levels = y_axis_group_order)) %>%
      mutate(!!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order)) %>%
      mutate(
        label_new = ifelse(
          feature_label %in% atcc_genera,
          paste0("<i>", feature_label, "</i> "),
          paste0(feature_label, " ")
        )
      ) %>%
      arrange(plot_order) %>%
      ungroup()

    ################################################################################
    message("260 # THEME AND FUNCTIONS")
    ################################################################################

    margin_size <- 10
    annotate_text_size <- 6 * magnify
    axis_text_size <- 16 * magnify
    axis_title_text_size <- 18 * magnify
    strip_text_size <- 20 * magnify
    legend_text_size <- 20 * magnify
    plot_width <- 10
    plot_height <- 10
    margin_buffer <- 30

    theme_common <- theme_global(base_size = 11) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = ggtext::element_markdown(face = "bold"),
        axis.title.y = ggtext::element_markdown(face = "bold"),
        legend.text = element_markdown(size = rel(1.2))
      )

    gPlot <- function(p) {
      p <- p +
        scale_color_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
        scale_fill_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
        labs(fill = "", color = "", title = "", caption = "") +
        guides(color = guide_legend(override.aes = list(linewidth = 2.5))) +
        theme_common

      print(p)
      p
    }

    ################################################################################
    message("320 # STRIP BACKGROUNDS")
    ################################################################################

    unique_y_labels <- unique(plot_df[[y_facet_group]])
    outline_color_y <- "black"
    s <- paste0(
      'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
      paste0(unique_y_labels, collapse = '"]),element_rect(color=outline_color_y,fill = palette_color["'),
      '"]))'
    )
    print(s)
    eval(parse(text = s))

    if (use_custom_labels == "yes") {
      custom_y_labels <- c("slope", "1500TF", "x3FC", "slope", "1500TF", "x3FC", "24cycles", rep("Soil", 3), rep("Zymo", 3), "NTC")
      unique_y_labels <- custom_y_labels
      s <- paste0(
        'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
        paste0(unique_y_labels, collapse = '"]),element_rect(color=outline_color_y,fill = palette_color["'),
        '"]))'
      )
      print(s)
      eval(parse(text = s))
    }

    y_strip <- as.data.frame(unique_y_labels) %>%
      mutate(text_color = ifelse(unique_y_labels == "Soil", "white", "black")) %>%
      mutate(text_color = "white") %>%
      mutate(face = "bold") %>%
      mutate(text_size = strip_text_size) %>%
      ungroup()

    unique_x_labels <- unique(plot_df[[x_facet_group]])
    outline_color_x <- "black"
    s <- paste0(
      'backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',
      paste0(unique_x_labels, collapse = '"]),element_rect(color=outline_color_x,fill = palette_color["'),
      '"]))'
    )
    print(s)
    eval(parse(text = s))

    if (use_custom_labels == "yes") {
      custom_x_labels <- c("slope", "1500TF", "x3FC", "slope", "1500TF", "x3FC", "24cycles", rep("Soil", 3), rep("Zymo", 3), "NTC")
      unique_x_labels <- custom_x_labels
      s <- paste0(
        'backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',
        paste0(unique_x_labels, collapse = '"]),element_rect(color=outline_color_x,fill = palette_color["'),
        '"]))'
      )
      print(s)
      eval(parse(text = s))
    }

    x_strip <- as.data.frame(unique_x_labels) %>%
      mutate(text_color = ifelse(unique_x_labels == "Soil", "white", "black")) %>%
      mutate(text_color = "white") %>%
      mutate(face = "bold") %>%
      mutate(text_size = strip_text_size) %>%
      ungroup()

    ################################################################################
    message("400 # BUILD PLOTS")
    ################################################################################

    sc <- abs(floor(min(plot_df$mean_rel_abun_log10)))
    x_min <- log10(2)
    x_alt <- 1
    magnify <- 3
    dodge_pos <- position_dodge(width = 0.8)

    assign_palette(levels(plot_ind_df$label_new))
    new_labels <- c("spin_column" = "Manual", "mag_beads" = "Automated")
    palette_label <- c(new_labels, palette_label)

    rel_log_plot <- ggplot(plot_ind_df, aes(x = !!sym(x_facet_group), y = rel_abun / 0.01, color = label_new, group = label_new)) +
      facet_grid(
        rows = vars(!!sym(y_axis_group)),
        labeller = labeller(!!sym(y_axis_group) := as_labeller(palette_label))
      ) +
      stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.07, linewidth = 0.7) +
      scale_y_log10() +
      labs(
        x = "Concentration of MMC",
        y = "Percent Relative Abundance"
      )

    rel_log_plot <- gPlot(rel_log_plot)

    counts_log_plot <- ggplot(plot_ind_df, aes(x = !!sym(x_facet_group), y = counts, color = label_new, group = label_new)) +
      facet_grid(
        rows = vars(!!sym(y_axis_group)),
        labeller = labeller(!!sym(y_axis_group) := as_labeller(palette_label))
      ) +
      stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
      stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.07, linewidth = 0.7) +
      scale_y_log10() +
      labs(
        x = "Concentration of MMC",
        y = "Number of Sequence Reads"
      )

    counts_log_plot <- gPlot(counts_log_plot)

    saveRDS(rel_log_plot, paste0(output_plot, "Figure-", figure_number, "_rel_log_plot.rds"))
    ggsave(
      paste0(output_plot, "Figure-", figure_number, "_rel_log_plot.png"),
      rel_log_plot,
      width = plot_width * 2,
      height = plot_height
    )

    saveRDS(counts_log_plot, paste0(output_plot, "Figure-", figure_number, "_counts_log_plot.rds"))
    ggsave(
      paste0(output_plot, "Figure-", figure_number, "_counts_log_plot.png"),
      counts_log_plot,
      width = plot_width * 2,
      height = plot_height
    )
  }
}

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
