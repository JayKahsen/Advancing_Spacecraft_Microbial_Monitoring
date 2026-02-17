source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
################################################################################
# Description
script_title='meta_data'

################################################################################
# Read in file names from Phylum table

# create meta.csv (meta data) for other scripts


################################################################################
# meta_data<<-c('sample_name','Exp','taxa_lev','Type','Rpl','variant','group','primer_label','Primers','Primer')
# 
# primer_names<-c('StandardO','TruncatedFO','Standard','TruncatedR','TruncatedF','TruncatedA','V1.2-12.CC','V2.2-13.CA','V3.3-14.TC','V4.4-15.TA','V5.5-16.xA','V6.6-15.xC','V7.7-17.xA','V8.8-17.xC')
# 
# primer_labels<-c('St','Tr','St','Tr','Tr','Tr','V1','V2','V3','V4','V5','V6','V7','V8')
# 
# primer_groups <- c(rep('Pool', 6), rep('Ind', 8))
# 
# recode_labels<-setNames(primer_labels,primer_names)
# recode_groups<-setNames(primer_groups,primer_names)
################################################################################
if (file.exists('data_tables/fastq_names.csv')) {
  fastq_names=read.csv('data_tables/fastq_names.csv')
  }else{
    print('getting filenames from fastq')
fastq_names <- data.frame(file_name = list.files(path = paste0(wd, '/demultiplexed_seqs'), pattern = "\\.fastq$", full.names = FALSE))
write.csv(fastq_names,'data_tables/fastq_names.csv',row.names=FALSE)
}

fastq_df=fastq_names%>%
  mutate(sample_ID = str_remove(file_name, "\\.fastq$"))%>%
  mutate(Temp=sample_ID)%>%
  separate(Temp, into=c('merge_numper','experiment','mastermix','sample_type','primer_concentration','plate_type','well'),sep='_')%>%
  mutate(short_name=paste(experiment,mastermix,sample_type,primer_concentration,plate_type,sep='_'))%>%
  mutate(sample_name=paste(short_name,well,sep='_'))

write.csv(fastq_df,'data_tables/meta.csv')
stop()


%>%
  #mutate(working = str_remove(sample_ID, "_.*"))%>%
  mutate(Exp = sapply(strsplit(sample_ID, "_"), function(x) x[1]))%>%
  mutate(well = sapply(strsplit(sample_ID, "_"), function(x) x[2]))%>%
  mutate(column = as.numeric(gsub("[^0-9]", "", well)))%>%
  mutate(row= gsub("[0-9]", "", well))%>%
  mutate(Temperature = '52')%>%
  mutate(primer = 'Standard')%>%
  mutate(sample_type = case_when(
    column %in% c(1,3,5,7,9,11) &  row %in% c('A','B','C','D')~ 'Zymo',
    column %in% c(2,4,6,8,10,12) &  row %in% c('A','B','C','D')~ 'Soil',
    TRUE~'NTC'))%>%
  mutate(mastermix = case_when(
    column %in% c(1,2,5,6,9,10)~ 'AMPED',
    column %in% c(3,4,7,8,11,12)~ 'QB',
    TRUE~'empty'))%>%
  mutate(replicate = case_when(
    row %in% c('A')~ '1',
    row %in% c('B')~ '2',
    row %in% c('C')~ '3',
    row %in% c('D')~ '4',
    row %in% c('E')~ 'N1',
    row %in% c('F')~ 'N2',
    row %in% c('G')~ 'N3',
    row %in% c('H')~ 'N4'
    
))%>%
  mutate(normalization = case_when(
    row %in% c('A','B','C','D') & column %in% c(1:4)~ 'slope',
    row %in% c('A','B','C','D') & column %in% c(5:8)~ 'x3FC',
    row %in% c('A','B','C','D') & column %in% c(9:12)~ '1500TF',
    TRUE~'24cycles'
  ))%>%
  mutate(short_name=paste(mastermix,sample_type,normalization,sep='_'))%>%
  mutate(sample_name=paste(mastermix,sample_type,normalization,well,sep='_'))%>%
    ungroup()
 
 

df=fastq_df


# 
# df <- read.delim(paste0('data_tables/original/Table_Phylum.tsv'),check.names = FALSE,row.names = 1, skip = 1)%>%
#   t%>%as.data.frame()%>%
#   rownames_to_column(var='sample_name')%>%
#   select(sample_name)%>%
#   left_join(fastq_df)%>%
#   # mutate(Exp = sapply(strsplit(sample_name, "_"), function(x) x[1]))%>%
#   # mutate(Type = sapply(strsplit(sample_name, "_"), function(x) x[2]))%>%
#   # mutate(Primer = sapply(strsplit(sample_name, "_"), function(x) x[3]))%>%
#   # mutate(Rpl = sapply(strsplit(sample_name, "_"), function(x) x[4]))%>%
# 
#   mutate(primer_label = recode(Primer, !!!recode_labels))%>%
#   mutate(Group = recode(Primer, !!!recode_groups))%>%
#   mutate(Forward_Primer=paste0(primer_label,'-tr515F'))%>%
#   mutate(Forward_Primer=case_when(
#     primer_label=='St'~ '515F',
#     primer_label=='Tr'~ 'tr515F',
#     TRUE~paste0(primer_label,'-tr515F')
#   ))%>%
#   mutate(Reverse_Primer = case_when(
#     primer_label == "Tr" & Exp %in% c("E1", "E2") ~ "tr806R",
#     TRUE ~ "806R"
#   ))%>%
#   mutate(Type=case_when(
#     Type=='Fecal'~'Feces',
#     Type=='Water'~'Wastewater',
#     TRUE~Type
#   ))%>%
#   mutate(sample_name=paste(Exp,Type,Forward_Primer,Reverse_Primer,Rpl,sep='_'))%>%
#   mutate(cName=paste(Exp,Type,primer_label,sep='_'))%>%
#   arrange(Primer)
#   
write.csv(df,paste0('data_tables/meta.csv'),row.names=FALSE)

################################################################################
