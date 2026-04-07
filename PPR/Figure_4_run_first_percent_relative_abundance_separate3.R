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
Creates the Figure 4 percent-relative-abundance plot for the PPR dataset.
Builds the stacked composition panel used in the combined Figure 4 release and
recovery output, then writes the review RDS and PNG files to output_plot.
"

#=============================== Plot settings ================================#

options(scipen = 999)

script_title <- "abundance for Figure 4"
figure_number <- "4"
starting_r <- 1
ending_r <- 1
testing <- "yes"
testing_r <- 3

#============================= Parameter settings =============================#

parameter_sets <- list(
  set1 = list(
    filter1_group = "figure",
    y_axis_group = "sampling_device",
    x_facet_group = "figure",
    y_facet_group = "experiment",
    number_size = 3
  )
)

#=================================== END SETUP ===================================#

################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

theme_common <- theme_global(base_size = 11) +
  theme(
    axis.text.y.right = ggtext::element_markdown(face = "bold"),
    axis.title.y.right = ggtext::element_markdown(face = "bold"),
    axis.text.y = ggtext::element_markdown(face = "bold"),
    axis.title.y = ggtext::element_markdown(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  )

gPlot <- function(p) {
  p <- p +
    scale_color_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
    scale_fill_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
    labs(fill = "", title = "", caption = "", y = "MMC Constituents") +
    guides(color = "none") +
    theme_common

  print(p)
  p
}

################################################################################
message("90 # MAIN")
################################################################################

#................................. loop settings ..............................#

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

    set_idx <- match(data_set, data_set_order)
    if (is.na(set_idx)) {
      stop("data_set not found in data_set_order: ", data_set)
    }

    params <- parameter_sets[[paste0("set", set_idx)]]
    list2env(params, envir = environment())

    groups <- c("filter1_group", "y_axis_group", "x_facet_group", "y_facet_group")
    for (var in groups) {
      group_value <- get(var)
      ord_name <- paste0(group_value, "_order")
      if (!exists(ord_name)) {
        stop("Missing required order vector: ", ord_name)
      }
      assign(paste0(var, "_order"), get(ord_name))
    }

    grouping_columns <- c(y_axis_group, x_facet_group, y_facet_group)
    run_suffix <- ""
    if (exists("Run_Group")) {
      run_suffix <- paste0("_", lp)
    }

    p_title <- paste0(taxa_levs, "_", script_title, custom_name, "_", data_set, run_suffix)
    qPrint(p_title)

    ################################################################################
    message("180 # LOAD AND FILTER")
    ################################################################################

    df_matrix1 <- read.csv(matrix_names[r, "file_path"], check.names = FALSE, row.names = 1) %>%
      t() %>%
      as.data.frame() %>%
      xPlode_sample_name() %>%
      filter(.data[[filter1_group]] %in% filter1_group_order) %>%
      filter(.data[[y_facet_group]] %in% y_facet_group_order) %>%
      filter(.data[[x_facet_group]] %in% x_facet_group_order) %>%
      filter(.data[[y_axis_group]] %in% y_axis_group_order) %>%
      filter(str_detect(Figure, figure_number))

    if (exists("Run_Group")) {
      df_matrix1 <- df_matrix1 %>%
        filter(.data[[Loop_Group]] %in% lp)
    }

    df_matrix2 <- df_matrix1 %>%
      imPlode_sample_name() %>%
      mutate(across(everything(), as.numeric)) %>%
      select(which(colSums(.) > 0))

    df_matrix <- combine_low_abundance(
      df = df_matrix2,
      threshold = 0.01,
      name = "(Consolidated < 1%)"
    )
    df_matrix <- combine_other(df = df_matrix, patterns = atcc_genera, starts_with = TRUE)
    feature_names <- names(df_matrix)

    ################################################################################
    message("250 # BUILD LONG DATA")
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
      mutate(
        rel_abun = counts / sample_counts,
        percent_rel_abun = rel_abun / 0.01,
        log10_percent_rel_abun = log10(percent_rel_abun),
        rel_abun_log10 = log10(rel_abun + 1e-9)
      ) %>%
      group_by(across(any_of(y_facet_group))) %>%
      mutate(group_counts = mean(sample_counts)) %>%
      ungroup() %>%
      mutate(
        group_counts_log10 = log10(group_counts + 1e-9),
        !!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_facet_group_order),
        !!sym(y_axis_group) := factor(!!sym(y_axis_group), levels = y_axis_group_order),
        !!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order),
        new_plot_group = ifelse(figure == "4A", "4A", "4BC"),
        label_new = factor(label_new, levels = rev(c("Other", setdiff(unique(label_new), "Other"))))
      ) %>%
      arrange(plot_order) %>%
      filter(label_new != "Other")

    plot_df <- plot_ind_df %>%
      group_by(feature, label_new, across(any_of(grouping_columns))) %>%
      summarize(
        group_mean_counts = mean(counts),
        sample_counts = mean(sample_counts),
        mean_log10_percent_rel_abun = mean(log10_percent_rel_abun),
        mean_percent_rel_abun = mean(percent_rel_abun),
        group_counts = mean(group_counts),
        .groups = "drop"
      ) %>%
      left_join(feature_labels_df %>% select(-c(taxa)), by = c("feature", y_facet_group)) %>%
      mutate(
        !!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_facet_group_order),
        !!sym(y_axis_group) := factor(!!sym(y_axis_group), levels = y_axis_group_order),
        !!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order)
      ) %>%
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
      select(rev(everything()))

    ################################################################################
    message("360 # BUILD st_plot")
    ################################################################################

    axis_text_size <- 16 * 3
    dodge_pos <- position_dodge(width = 0.8)
    facet_labels <- c(
      "4A" = "D\\) MMC Released from Swab Heads",
      "4B" = "E\\) MMC Recovered by Swabs",
      "4C" = "F\\) MMC Rec. by Polyester Wipes and SALSA"
    )
    palette_label_plot <- c(facet_labels, palette_label)
    y_labels <- setNames(as.character(plot_df$label_new), plot_df$plot_order)

    st_plot_base <- ggplot(plot_df, aes(x = mean_percent_rel_abun, y = plot_order)) +
      geom_col(
        aes(fill = !!sym(y_axis_group)),
        width = 0.6,
        position = dodge_pos,
        alpha = 0.7
      ) +
      geom_jitter(
        data = plot_ind_df,
        aes(y = plot_order, x = rel_abun / 0.01, fill = !!sym(y_axis_group)),
        position = position_jitterdodge(
          dodge.width = 0.8,
          jitter.width = 0.4,
          jitter.height = 0
        ),
        color = "black",
        shape = 21,
        size = 2
      ) +
      facet_grid2(
        formula(paste("~", x_facet_group)),
        labeller = labeller(
          .rows = as_labeller(palette_label_plot),
          .cols = as_labeller(palette_label_plot)
        ),
        space = "free"
      )

    st_plot <- st_plot_base +
      geom_text(
        aes(
          label = round(mean_percent_rel_abun, 1),
          x = ifelse(mean_percent_rel_abun < 50, 90, mean_percent_rel_abun),
          color = !!sym(y_axis_group),
          group = !!sym(y_axis_group)
        ),
        size = number_size,
        hjust = 1.1,
        fontface = "bold",
        show.legend = FALSE,
        position = position_dodge(width = 0.8)
      ) +
      scale_y_discrete(labels = y_labels, position = "right") +
      theme(
        axis.text.y = ggtext::element_markdown(face = "bold", size = axis_text_size),
        strip.text = ggtext::element_markdown(size = rel(1.2), face = "bold")
      )

    ################################################################################
    message("450 # APPLY gPlot AND SAVE")
    ################################################################################

    st_plot <- gPlot(st_plot) +
      labs(title = "", x = "Mean Percent Relative Abundance") +
      guides(
        color = guide_legend(reverse = TRUE),
        fill = guide_legend(reverse = TRUE)
      )

    saveRDS(st_plot, paste0(output_plot, "Figure_", figure_number, "_st_plot.rds"))
    ggsave(paste0(output_plot, "Figure_", figure_number, "_st_plot.png"), st_plot,width=15)
  }
}

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
