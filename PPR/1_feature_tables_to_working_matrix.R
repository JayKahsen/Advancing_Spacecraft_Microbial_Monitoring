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
Creates data_tables/matrix_names.csv as an index of the PPR working matrices.

Reads the merged feature tables, aligns them to the metadata,
and writes standardized working matrices for downstream analysis.

Summarizes per-sample reads in data_tables/sample_reads.csv.
"


#=============================== Plot settings ================================#

script_title <- "feature_tables_to_working_matrix"


#=================================== END SETUP ===================================#


################################################################################
message("60 # USER CONTROLS")
################################################################################

making_matrix_names <- "yes"
making_working_matrix <- "yes"
working_matrix_path <- "data_tables/working_matrix"
create_directory(working_matrix_path)

active_taxa_levels <- intersect(c("Phylum", "Family", "Genus", "Species"), taxa_levels)

if (!exists("data_sets")) {
  data_sets <- "s0"
  s0 <- "s0"
}


################################################################################
message("80 # MATRIX NAMES")
################################################################################

if (making_matrix_names == "yes") {
  names(taxa_plural) <- taxa_levels

  matrix_names <- expand.grid(taxa_levs = active_taxa_levels, data_sets = data_sets) %>%
    mutate(
      taxa_plural = case_when(
        taxa_levs %in% names(taxa_plural) ~ taxa_plural[taxa_levs],
        TRUE ~ taxa_levs
      )
    ) %>%
    mutate(working_path = file.path(working_matrix_path, paste0("working_matrix_", taxa_levs, ".csv"))) %>%
    mutate(file_name = working_path) %>%
    mutate(file_path = working_path) %>%
    arrange(taxa_levs)

  write.csv(matrix_names, "data_tables/matrix_names.csv", row.names = FALSE)
  message("created data_tables/matrix_names.csv")
}


################################################################################
message("110 # WORKING MATRICES")
################################################################################

if (making_working_matrix == "yes") {
  sample_reads_df2 <- NULL

  for (taxa_levs in active_taxa_levels) {
    qPrint(taxa_levs)
    source_file <- paste0("data_tables/original/Merged_", taxa_levs, "_raw_abund.csv")

    if (!file.exists(source_file)) {
      message(paste("missing source table:", source_file))
      next
    }

    df_matrix3 <- read.csv(source_file, row.names = 1, check.names = FALSE) %>%
      column_to_rownames(var = "ID") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = "sample_id_PPR") %>%
      left_join(meta) %>%
      select(sample_id_PPR, sample_name, everything())

    df_matrix2 <- df_matrix3 %>%
      filter(!is.na(sample_name)) %>%
      imPlode_sample_name()

    df_matrix2[is.na(df_matrix2)] <- 0

    sample_reads <- rowSums(df_matrix2) %>%
      as.data.frame() %>%
      setNames("Reads") %>%
      mutate(taxa_level = paste0(taxa_levs, "_Reads")) %>%
      xPlode_sample_name()

    sample_reads_df2 <- rbind(sample_reads_df2, sample_reads)

    df_matrix2 %>%
      t() %>%
      as.data.frame() %>%
      write.csv(file.path(working_matrix_path, paste0("working_matrix_", taxa_levs, ".csv")))
  } # end loop over taxa levels

  if (!is.null(sample_reads_df2)) {
    sample_reads_df <- sample_reads_df2 %>%
      pivot_wider(names_from = taxa_level, values_from = Reads) %>%
      imPlode_sample_name()

    write.csv(sample_reads_df, "data_tables/sample_reads.csv")
    message("created data_tables/sample_reads.csv")
  }

  message("created working matrices in data_tables/working_matrix")
}


################################################################################
message("170 # END")
################################################################################
