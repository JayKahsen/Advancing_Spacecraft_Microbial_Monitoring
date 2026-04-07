################################################################################
message("10 # SETUP")
################################################################################


#=========================== Script path resolution ============================#
message("20 # Script path resolution")

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


#============================= Description =============================#
script_description='
Creates data_tables/matrix_names.csv as an index of the KSC working matrices.

Reads the merged feature tables, aligns them to the metadata,
and writes standardized working matrices for downstream analysis.

Summarizes per-sample reads in data_tables/sample_reads.csv.
'

script_title <- "feature_tables_to_working_matrix"
message("The source preprocessing steps are usually already completed.")
message("Rerun this script only when merged feature tables or metadata change.")


#================================ END SETUP ================================#


################################################################################
message("60 # USER CONTROLS")
################################################################################

making_matrix_names <- "yes"
making_working_matrix <- "yes"
working_matrix_path <- "data_tables/working_matrix"

create_directory(working_matrix_path)


################################################################################
message("80 # MATRIX NAMES")
################################################################################

if (making_matrix_names == "yes") {
  names(taxa_plural) <- taxa_levels

  matrix_names <- tibble(taxa_levs = taxa_levels) %>%
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

  for (taxa_levs in taxa_levels) {
    qPrint(taxa_levs)
    source_file <- file.path("data_tables", "original", paste0("Merged_", taxa_levs, "_raw_abund.csv"))

    if (!file.exists(source_file)) {
      message(paste("missing source table:", source_file))
      next
    }

    df_matrix2 <- read.csv(source_file, check.names = FALSE, row.names = 1) %>%
      column_to_rownames(var = "ID")

    df_matrix2[is.na(df_matrix2)] <- 0

    df_matrix1 <- df_matrix2 %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = "sample_id") %>%
      left_join(meta) %>%
      select(any_of(meta_names), everything()) %>%
      filter(!is.na(sample_name)) %>%
      imPlode_sample_name() %>%
      ungroup()

    sample_reads <- rowSums(df_matrix1) %>%
      as.data.frame() %>%
      setNames("Reads") %>%
      mutate(taxa_level = paste0(taxa_levs, "_Reads")) %>%
      xPlode_sample_name()

    sample_reads_df2 <- rbind(sample_reads_df2, sample_reads)

    df_matrix1 %>%
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
