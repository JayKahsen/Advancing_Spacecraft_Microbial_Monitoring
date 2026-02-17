source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
################################################################################
# Description
################################################################################
script_title='ASV_taxa'
# make ASV_taxa for other scripts
# need ASV and taxa data from
# QIIME2
# classified_rep_seqs_filtered_v2.qzv download as metadata resave as taxa_data.tsv
################################################################################
# get sequence data
################################################################################
  #library(insect)

  # sequence_data=readFASTA('data_tables/sequences.fasta')%>%
  #   dna2char()%>%
  #   as.data.frame()%>%
  #   `names<-`("Sequence") %>%
  #   as.data.frame()%>%
  #   rownames_to_column(var = "ASV")
  ################################################################################

  ASV_taxa=read.delim('data_tables/original/taxa_data.tsv')%>%
    slice(-1)%>%
    setNames(c('ASV','taxa','Confidence'))%>%
    select(ASV,taxa)%>%
    mutate(Domain = str_extract(taxa, "d__([^;\\s]*)"),
           Phylum = str_extract(taxa, "p__([^;\\s]*)"),
           Class = str_extract(taxa, "c__([^;\\s]*)"),
           Order = str_extract(taxa, "o__([^;\\s]*)"),
           Family = str_extract(taxa, "f__([^;\\s]*)"),
           Genus = str_extract(taxa, "g__([^;\\s]*)"),
           Species = str_extract(taxa, "s__([^;\\s]*)"))%>%
    mutate(default = case_when(
      !is.na(Species) ~ Species,
      !is.na(Genus) ~ Genus,
      !is.na(Family) ~ Family,
      !is.na(Order) ~ Order,
      !is.na(Class) ~ Class,
      !is.na(Phylum) ~ Phylum,
      !is.na(Domain) ~ Domain,
      TRUE ~ ASV))%>%
    mutate(taxa = case_when(is.na(taxa) ~ ASV,TRUE ~ taxa),
           Domain = case_when(is.na(Domain) ~ taxa,TRUE ~ Domain),
           Phylum = case_when(is.na(Phylum) ~ Domain,TRUE ~ Phylum),
           Class = case_when(is.na(Class) ~ Phylum,TRUE ~ Class),
           Order = case_when(is.na(Order) ~ Class,TRUE ~ Order),
           Family = case_when(is.na(Family) ~ Order,TRUE ~ Family),
           Genus = case_when(is.na(Genus) ~ Family,TRUE ~ Genus),
           Species = case_when(is.na(Species) ~ Genus,TRUE ~ Species))
  
  write.csv(ASV_taxa,'data_tables/ASV_taxa.csv',row.names = FALSE)
  


