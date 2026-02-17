# Feature-table standardization workflow (JPL)
#
# This script keeps only the required operations:
#   1) build/update data_tables/matrix_names.csv
#   2) transform merged original feature tables into standardized raw matrices
#
# It intentionally does NOT generate filtered/rarefied outputs.

script_dir <- tryCatch(
  dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),
  error = function(e) getwd()
)
source(file.path(script_dir, '_globalStuff.R'))

script_title <- 'feature_tables_to_raw_matrix'

# run switches
making_matrix_names <- 'yes'
making_raw_matrix <- 'yes'

# output paths used in this reduced workflow
raw_matrix_path <- 'data_tables/raw_matrix'
create_directory(raw_matrix_path)

# use only canonical merged taxonomic inputs
active_taxa_levels <- c('Phylum', 'Family', 'Genus', 'Species')
active_taxa_levels <- intersect(active_taxa_levels, taxa_levels)

if (!exists('data_sets')) {
  data_sets <- 's0'
  s0 <- 's0'
}

if (making_matrix_names == 'yes') {
  names(taxa_plural) <- taxa_levels

  matrix_names <- expand.grid(taxa_levs = active_taxa_levels, data_sets = data_sets) %>%
    mutate(taxa_plural = case_when(
      taxa_levs %in% names(taxa_plural) ~ taxa_plural[taxa_levs],
      TRUE ~ taxa_levs
    )) %>%
    mutate(raw_path = file.path(raw_matrix_path, paste0('raw_matrix_', taxa_levs, '.csv'))) %>%
    mutate(file_name = raw_path) %>%
    mutate(file_path = raw_path) %>%
    arrange(taxa_levs)

  write.csv(matrix_names, 'data_tables/matrix_names.csv', row.names = FALSE)
  message('created data_tables/matrix_names.csv')
}

if (making_raw_matrix == 'yes') {
  sample_reads_df2 <- NULL

  for (taxa_levs in active_taxa_levels) {
    source_file <- file.path('data_tables', 'original', paste0('Merged_', taxa_levs, '_raw_abund.csv'))

    if (!file.exists(source_file)) {
      warning('Missing source table: ', source_file)
      next
    }

    df_matrix2 <- read.csv(source_file, check.names = FALSE, row.names = 1) %>%
      column_to_rownames(var = 'ID')

    df_matrix2[is.na(df_matrix2)] <- 0

    df_matrix1 <- df_matrix2 %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = 'sample_id') %>%
      left_join(meta) %>%
      select(any_of(meta_names), everything()) %>%
      filter(!is.na(sample_name)) %>%
      imPlode_sample_name() %>%
      ungroup()

    sample_reads <- rowSums(df_matrix1) %>%
      as.data.frame() %>%
      setNames('Reads') %>%
      mutate(taxa_level = paste0(taxa_levs, '_Reads')) %>%
      xPlode_sample_name()

    sample_reads_df2 <- rbind(sample_reads_df2, sample_reads)

    df_matrix1 %>%
      t() %>%
      as.data.frame() %>%
      write.csv(file.path(raw_matrix_path, paste0('raw_matrix_', taxa_levs, '.csv')))
  }

  if (!is.null(sample_reads_df2)) {
    sample_reads_df <- sample_reads_df2 %>%
      pivot_wider(names_from = taxa_level, values_from = Reads) %>%
      imPlode_sample_name()

    write.csv(sample_reads_df, 'data_tables/sample_reads.csv')
    message('created data_tables/sample_reads.csv')
  }

  message('created standardized raw matrices in ', raw_matrix_path)
}
