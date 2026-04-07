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
Creates the Figure 2 percent-relative-abundance plot for the PPR dataset.
Builds the stacked abundance panel from the selected working matrix and writes
the RDS output used by the companion Figure 2 script.
"

#=============================== Plot settings ================================#

options(scipen = 999)

script_title <- "abundance for figure 2"
figure_number <- "2"
starting_r <- 1
ending_r <- 1
testing <- "yes"
testing_r <- 3

#============================= Parameter settings =============================#

parameter_sets <- list(
  set1 = list(
    filter1_group = "figure",
    y_axis_group = "bar_color",
    x_facet_group = "bar_type_2",
    y_facet_group = "experiment",
    number_size = 3
  )
)

#=================================== END SETUP ===================================#

################################################################################
message("60 # THEME AND FUNCTIONS")
################################################################################

gPlot <- function(p) {
  p <- p +
    scale_color_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
    scale_fill_manual(values = palette_color, labels = palette_label, na.translate = FALSE) +
    labs(fill = "", title = "", caption = "", y = "") +
    guides(color = "none") +
    theme_common

  print(p)
  p
}

################################################################################
message("90 # MAIN")
################################################################################

#................................. sanity checks ................................#

need_objs <- c("matrix_names", "data_set_order")
for (nm in need_objs) {
  if (!exists(nm)) {
    stop("Missing required object: ", nm)
  }
}

need_funs <- c("xPlode_sample_name", "imPlode_sample_name", "combine_low_abundance", "feature_labels", "theme_global")
for (fn in need_funs) {
  if (!exists(fn)) {
    stop("Missing required function: ", fn)
  }
}

#............................... theme settings ...............................#

magnify <- 3
axis_text_size <- 16 * magnify

theme_common <- theme_global(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  )

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

  custom_name <- if (exists("custom_name")) custom_name else ""
  p_title <- paste0(taxa_levs, "_", script_title, custom_name, "_", data_set)
  qPrint(p_title)

  ################################################################################
  message("200 # LOAD AND FILTER")
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

  df_matrix2 <- df_matrix1 %>%
    imPlode_sample_name() %>%
    mutate(across(everything(), as.numeric)) %>%
    select(which(colSums(.) > 0))

  df_matrix <- combine_low_abundance(df = df_matrix2, threshold = 0.01, name = "(Consolidated < 1%)")
  feature_names <- names(df_matrix)

  ################################################################################
  message("240 # BUILD LONG DATA")
  ################################################################################

  feature_labels_df <- feature_labels(y_facet_group, averaging_group = x_facet_group)
  y_axis_group_order <- rev(y_axis_group_order)

  plot_ind_df <- df_matrix %>%
    xPlode_sample_name() %>%
    pivot_longer(cols = any_of(feature_names), names_to = "feature", values_to = "counts") %>%
    left_join(feature_labels_df %>% select(-c(taxa)), by = c("feature", y_facet_group)) %>%
    group_by(feature, across(any_of(y_facet_group))) %>%
    mutate(grouped_feature_counts = sum(counts)) %>%
    filter(grouped_feature_counts > 0) %>%
    group_by(sample_name) %>%
    mutate(sample_counts = sum(counts)) %>%
    ungroup() %>%
    mutate(rel_abun = counts / sample_counts) %>%
    ungroup() %>%
    arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
    mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_axis_group_order)) %>%
    mutate(!!sym(y_axis_group) := factor(!!sym(y_axis_group), levels = y_axis_group_order)) %>%
    mutate(!!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order)) %>%
    distinct()

  plot_df <- plot_ind_df %>%
    group_by(feature, across(any_of(c(y_axis_group, x_facet_group, y_facet_group)))) %>%
    summarize(
      mean_rel_abun = mean(rel_abun),
      .groups = "drop"
    ) %>%
    left_join(feature_labels_df %>% select(-c(taxa)), by = c("feature", y_facet_group)) %>%
    arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
    mutate(
      label_new = ifelse(
        feature_label %in% c(
          "(Consolidated < 1%)",
          "Unclassified",
          "Unclassified Bacilli",
          "Unclassified Bacteria",
          "Unclassified Paracoccaceae",
          "Unclassified Propionibacteriaceae",
          "Unmapped"
        ),
        feature_label,
        paste0("*", feature_label, "*")
      )
    ) %>%
    mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_axis_group_order)) %>%
    mutate(!!sym(y_axis_group) := factor(!!sym(y_axis_group), levels = y_axis_group_order)) %>%
    mutate(!!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order))

  ################################################################################
  message("340 # BUILD st_plot")
  ################################################################################

  dodge_pos <- position_dodge(width = 0.8)
  y_labels <- setNames(plot_df$label_new, plot_df$plot_order)

  st_plot <- ggplot(plot_df, aes(x = mean_rel_abun / 0.01, y = plot_order)) +
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
        .rows = as_labeller(palette_label),
        .cols = as_labeller(palette_label)
      ),
      space = "free"
    ) +
    geom_text(
      aes(
        label = round(mean_rel_abun / 0.01, 1),
        x = ifelse(mean_rel_abun / 0.01 < 90, 90, mean_rel_abun / 0.01),
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
    theme(axis.text.y = ggtext::element_markdown(face = "bold", size = axis_text_size))

  ################################################################################
  message("410 # APPLY gPlot + SAVE")
  ################################################################################

  st_plot <- gPlot(st_plot) +
    labs(title = "Percent Relative Abundance", x = "")

  saveRDS(st_plot, paste0(output_plot, "Figure-", figure_number, "_st_plot.rds"))
  print(st_plot)
}

################################################################################
message(paste("620 # FINISHED", script_title))
################################################################################
