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
library(readxl)

# sample_id_XT=read.delim('data_tables/original/PPR_XT/Project_1790_PPR_XT_Output.family.allfilter.abund.tsv',row.names = 1,check.names = FALSE) %>% 
#   t() %>% as.data.frame() %>% 
#   rownames_to_column(var='sample_id_XT') %>% 
#   mutate(sample_id=sample_id_XT) %>% 
#   select(sample_id_XT,sample_id) %>% 
#   mutate(
#     sample_id = str_remove(sample_id, "^.*?_PPR_"),
#     sample_id = str_remove(sample_id, "_S\\d{3}.*$")
#   )
# 
# sample_id_XT_2=read.delim('data_tables/original/PPR_XT_2/Project_1790_PPR_XT_2_Output.family.allfilter.abund.tsv',row.names = 1,check.names = FALSE) %>% 
#   t() %>% as.data.frame() %>% 
#   rownames_to_column(var='sample_id_XT_2') %>% 
#   mutate(sample_id=sample_id_XT_2) %>% 
#   select(sample_id_XT_2,sample_id) %>% 
#   mutate(
#     sample_id = str_remove(sample_id, "^.*?_PPR_"),
#     sample_id = str_remove(sample_id, "_S\\d{3}.*$")
#   )

sample_id_PPR=read.csv('data_tables/original/Merged_Phylum_raw_abund.csv',row.names = 1,check.names = FALSE) %>% 
  column_to_rownames(var='ID') %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column(var='sample_id_PPR') %>%
  mutate(sample_name=sample_id_PPR) %>%
  select(sample_name,sample_id_PPR) %>%
  mutate(
    sample_name = str_remove(sample_id_PPR, "^.*?_PPR_"),
    sample_name = str_remove(sample_name, "_S\\d+.*$"),
    sample_name=str_remove(sample_name, "^1790XT_[0-9]+_")
  ) %>% 
  mutate(sample_name=ifelse(sample_name=='12_pcc_400_Qi_R_4','12_PCc_400_Qi_R_4',sample_name)) %>% 
  mutate(sample_name = case_when(
    sample_id_PPR == "1790XT_59_3_A3_R_1_S124" ~ "3_A2_R_1",
  #   sample_name == "PC4" ~ "PC4",
  #   sample_name == "NTC10" ~ "NTC_10",
  #   sample_name == "NTC11" ~ "NTC_11",
  #   sample_name == "NTC12" ~ "NTC_12",
  #   sample_name == "PC5" ~ "PC5",
  #   sample_name == "NTC13" ~ "NTC_13",
  #   sample_name == "NTC14" ~ "NTC_14",
  #   sample_name == "NTC15" ~ "NTC_15",
  #   sample_name == "PC6" ~ "PC6",
  #   sample_name == "NTC16" ~ "NTC_16",
  #   sample_name == "NTC17" ~ "NTC_17",
  #   sample_name == "NTC18" ~ "NTC_18",
  #   sample_name == "PC1" ~ "PC1",
  #   sample_name == "NTC01" ~ "NTC_01",
  #   sample_name == "NTC02" ~ "NTC_02",
  #   sample_name == "NTC03" ~ "NTC_03",
  #   sample_name == "PC2" ~ "PC2",
  #   sample_name == "NTC04" ~ "NTC_04",
  #   sample_name == "NTC05" ~ "NTC_05",
  #   sample_name == "NTC06" ~ "NTC_06",
  #   sample_name == "MaxBlank" ~ "MaxBlank",
  #   sample_name == "MaxWaterBlank" ~ "MaxWaterBlank",
  #   sample_name == "PC3" ~ "PC3",
  #   sample_name == "NTC07" ~ "NTC_07",
  #   sample_name == "NTC08" ~ "NTC_08",
  #   sample_name == "NTC09" ~ "NTC_09",
  #   sample_name == "zymo_pc4" ~ "zymo_pc4",
  #   sample_name == "zymo_pc5" ~ "zymo_pc5",
  #   sample_name == "NTC6" ~ "NTC6",
  #   sample_name == "NTC7" ~ "NTC7",
  #   sample_name == "NTC8" ~ "NTC8",
  #   sample_name == "NTC9" ~ "NTC9",
  #   sample_name == "2_PCM_R_2_repeat" ~ "2_PCM_R_2_repeat",
  #   sample_name == "3_PCM_R_1_repeat" ~ "3_PCM_R_1_repeat",
  #   sample_name == "NTC2" ~ "NTC2",
  #   sample_name == "NTC3" ~ "NTC3",
  #   sample_name == "NTC4" ~ "NTC4",
  #   sample_name == "NTC5" ~ "NTC5",
  #   sample_name == "zymo_pc2" ~ "zymo_pc2",
  #   sample_name == "zymo_pc3" ~ "zymo_pc3",
  #   sample_name == "1790XT_NTC1" ~ "1790XT_NTC1",
  #   sample_name == "1790XT_zymo_pc1" ~ "1790XT_zymo_pc1",
     TRUE ~ sample_name   # keeps everything else unchanged
  ))

sample_id_PPR$sample_name[duplicated(sample_id_PPR$sample_name)]
meta$sample_name[duplicated(meta$sample_name)]

sample_reads=read.csv('data_tables/sample_reads.csv') %>% 
  select(X,Phylum_Reads) %>% 
  rename(sample_name=X,sample_reads=Phylum_Reads)
  




#qLoad('data_tables/sample_reads.csv')

df0<- read_excel("data_tables/RUSH Master Sample Metadata 15JUL25.xlsx",sheet = "Sample Manifest") 
names(df0)

df2<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 2")
names(df2)

df_barnames <- df2 %>%
  select(`Bar Name...29`,unique_replicate_id) %>%
    rename(bar_name = `Bar Name...29`) %>% 
  mutate(bar_name = ifelse(
    str_detect(coalesce(unique_replicate_id, ''), 'NTC'),
    'NTC',
    bar_name
  )) %>% 
  #mutate(sample_id = str_replace_all(sample_id, "-", "_")) %>% 
   
 # filter(!is.na(bar_name)) %>% 
  ungroup()


names(df2)

df3<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 3") 
df4<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 4") 

df4A <- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx", sheet = "Figure 4A") %>%
  select(`RUSH ID#`, `Percent Recovery...17`, `Percent Recovery...23`) %>%
  rename(
    sample_id = `RUSH ID#`,
    percent_recovery_qpcr = `Percent Recovery...17`,
    percent_recovery_dpcr = `Percent Recovery...23`
  ) %>% 
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))


df4B<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 4B") %>% 
  select(`RUSH ID#`, `Percent Recovery...23`, `Percent Recovery...27`) %>%
  rename(
    sample_id = `RUSH ID#`,
    percent_recovery_qpcr = `Percent Recovery...23`,
    percent_recovery_dpcr = `Percent Recovery...27`
  ) %>% 
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))

df4C<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 4C") %>% 
  select(`RUSH ID#`, `Percent Recovery...23`, `Percent Recovery...27`) %>%
  rename(
    sample_id = `RUSH ID#`,
    percent_recovery_qpcr = `Percent Recovery...23`,
    percent_recovery_dpcr = `Percent Recovery...27`
  ) %>% 
  mutate(sample_id = str_replace_all(sample_id, "-", "_"))


df5<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 5") 
df6<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 6") 
df6A<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 6A") 
df6B<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 6B") 

df_percent_recovery=bind_rows(df4A,df4B,df4C) %>% 
  filter(sample_id!='#N/A') %>% 
  filter(percent_recovery_dpcr!='NA') %>% 
  filter(!is.na(percent_recovery_qpcr)) 

qLoad('data_tables/figure_2_selection.csv')
qLoad('data_tables/figure_3_selection.csv')
qLoad('data_tables/figure_4_selection.csv')
qLoad('data_tables/figure_5_selection.csv')

figure_selection <- meta%>% 
  select(sample_id) %>% 
left_join(figure_2_selection) %>% 
  left_join(figure_3_selection) %>% 
  left_join(figure_4_selection) %>% 
  left_join(figure_5_selection) %>% 
  mutate(
    figure = pmap_chr(
      select(., -sample_id),
      ~ paste(na.omit(c(...)), collapse = "_")
    )
  ) %>% 
  select(sample_id,figure) %>% 
  distinct()





df_meta=df0 %>% 
  mutate(sample_id=`RUSH ID #`,
         sample_id = str_replace_all(sample_id, "-", "_"),
         sample_id=ifelse(is.na(sample_id),unique_replicate_id,sample_id),
         copies_rxn_qpcr=`16s_rrna_copy_rxn`, #fig 3
         copies_rxn_dpcr=`Total Copies in Undiluted Sample (1 µL) [Dil*dPCR Copies]`,#fig 3
         copies_original_qpcr=`16s_rrna_copy_sample`,
         copies_original_dpcr=`Total Copies in Orginal Sample (50uL) [dPCR copies*dil*50]`
  ) %>%
mutate(across(
  c(copies_rxn_qpcr, copies_rxn_dpcr,
    copies_original_qpcr, copies_rxn_dpcr),
  ~ suppressWarnings(ifelse(is.na(.) | !is.numeric(.) | . <= 0, NA, log10(as.numeric(.)))),
  .names = "{.col}_log10"
)) %>% 
  mutate(sample_type=case_when(
    cell_concentration=='20'~'2.00E+1',
    cell_concentration=='200'~'2.00E+2',
    cell_concentration=='2000'~'2.00E+3',
    cell_concentration=='20000'~'2.00E+4',
    cell_concentration=='200000'~'2.00E+5',
    cell_concentration=='2000000'~'2.00E+6',
    cell_concentration=='20000000'~'2.00E+7',
    control_type %in% c('reagent_control','water_control','no_template_control')~control_type,
    TRUE~cell_concentration
  ))%>%
  mutate(across(
    c(copies_rxn_qpcr, copies_rxn_dpcr,
      copies_original_qpcr, copies_original_dpcr),
    ~ {
      num_val <- suppressWarnings(as.numeric(.))
      num_val[num_val < 0] <- NA         # negatives become NA
      num_val[num_val == 0] <- 1         # zeros become 1 so log10 = 0
      log10(num_val)
    },
    .names = "{.col}_log10"
  )) %>% 

  mutate(sample_type=ifelse(experiment_type=='extract_efficiency',control_type,sample_type)) %>% 
left_join(df_percent_recovery)%>%
left_join(df_barnames) %>% 
 left_join(figure_selection) %>%
  mutate(bar_name = ifelse(
    str_detect(coalesce(Figure, ''), '2') & str_detect(coalesce(sample_id, ''), 'NTC'),
    'NTC',
    bar_name
  )) %>% 
  mutate(sample_type=ifelse((str_detect(Figure, '3') & str_detect(sample_id, 'NTC')),'NTC',sample_type)) %>% 
  mutate(bar_axis = case_when(
    str_detect(bar_name, "^Cotton Dry") ~ "Cotton Dry",
    str_detect(bar_name, "^Macrofoam Dry") ~ "Macrofoam Dry",
    str_detect(bar_name, "^Cotton Wet") ~ "Cotton Wet",
    str_detect(bar_name, "^Macrofoam Wet") ~ "Macrofoam Wet",
    str_detect(bar_name, "^Cotton") ~ "Cotton",
    str_detect(bar_name, "^Macrofoam") ~ "Macrofoam",
    str_detect(bar_name, "^Water") ~ "Water",
    str_detect(bar_name, "^NTC") ~ "NTC",
    TRUE ~ NA_character_
  )) %>% 
  mutate(bar_axis_2 = case_when(
    str_detect(bar_name, "^Cotton Dry") ~ "Cotton",
    str_detect(bar_name, "^Macrofoam Dry") ~ "Macrofoam",
    str_detect(bar_name, "^Cotton Wet") ~ "Cotton",
    str_detect(bar_name, "^Macrofoam Wet") ~ "Macrofoam",
    str_detect(bar_name, "^Cotton") ~ "Cotton",
    str_detect(bar_name, "^Macrofoam") ~ "Macrofoam",
    str_detect(bar_name, "^Water") ~ "Water",
    str_detect(bar_name, "^NTC") ~ "NTC",
    TRUE ~ NA_character_
  )) %>% 
  mutate(bar_color = case_when(
    str_detect(bar_name, "^Cotton") ~ "Cotton",
    str_detect(bar_name, "^Macrofoam") ~ "Macrofoam",
    str_detect(bar_name, "^Water") ~ "Water",
    str_detect(bar_name, "^NTC") ~ "NTC",
    TRUE ~ NA_character_
  )) %>% 
  mutate(bar_type = case_when(
    str_detect(bar_name, "Environmental") ~ "Environmental",
    str_detect(bar_name, "Device") ~ "Device",
    str_detect(bar_name, "NTC") ~ "Mastermix",
    str_detect(bar_name, "Metal") ~ "Metal",
    str_detect(bar_name, "Water") ~ "Mastermix",
    TRUE ~ NA_character_
  )) %>% 
  mutate(bar_type_2 = case_when(
    bar_type=="Metal" & str_detect(bar_axis, "Wet") ~ "Metal Wet",
    bar_type=="Metal" & str_detect(bar_axis, "Dry") ~ "Metal Dry",
    TRUE ~ bar_type
  )) %>% 
  
  mutate(bar_fill = case_when(
    str_detect(bar_name, "^Cotton") ~ "gray1",
    str_detect(bar_name, "^Macrofoam") ~ "gray3",
    str_detect(bar_name, "^Water") ~ "gray1",
    str_detect(bar_name, "^NTC") ~ "gray3",
    sampling_device %in% c('cotton','Polyester Wipe')~'gray1',
    sampling_device %in% c('macrofoam','SALSA')~'gray3',
    TRUE ~ NA_character_
  )) %>% 
  mutate(sample_name=sample_id) %>%
# left_join(sample_id_XT) %>%
 #left_join(sample_id_XT_2) %>%
  left_join(sample_id_PPR) %>%
  left_join(sample_reads) %>%
  select(rev(everything())) %>% 
  
# select(
#   sample_id,sample_id_XT,sample_id_XT_2,sample_reads, sample_type,
#   experiment_type,figure,Figure,
#   percent_recovery_qpcr, percent_recovery_dpcr,
#   copies_rxn_qpcr, copies_rxn_dpcr,
#   copies_original_qpcr,  copies_original_dpcr,
#   copies_rxn_qpcr_log10, copies_rxn_dpcr_log10,
#   copies_original_qpcr_log10,  copies_original_dpcr_log10,
#   extraction_method,sampling_device,
#   bar_name, bar_type, bar_color, bar_fill,bar_axis, control_type,sample_name
#   ,everything()
#   )%>% 

  select(sample_id_PPR,sample_name,unique_replicate_id,`RUSH ID #`,batch_id,sample_id,replicate_id,
       #  `Sequencing Name ID `,
         `RUSH Project Batch Name`,
         figure,Figure,`George Figure number`,
         unique_replicate_id,experiment_type,sampling_device,bar_name,location,environment,extraction_method,control_type,environmental_sample,control_sample,dna_kit,microbial_type,everything()) %>% 
  arrange((experiment_type)) %>% 
  mutate(experiment='PPR') %>% 
  distinct() %>% 
  ungroup()

names(df_meta)
df_check=df_meta %>% 
  filter(str_detect(figure, '2'))
names(df_meta)
names(df4A)
names(df4B)
write.csv(df_meta,'data_tables/meta.csv',row.names=FALSE)

not_matched=setdiff(sample_id_PPR$sample_id_PPR,df_meta$sample_id_PPR)

sample_check=df_meta %>% 
  filter(is.na(sample_id_PPR)) %>% 
  write.csv('unmatched_meta.csv',row.names=FALSE)

sample_check$sample_name

remaining=sample_id_PPR %>% 
  filter(sample_id_PPR %in% not_matched) %>% 
  write.csv('sample_id_PPR_not_matched.csv',row.names=FALSE)


meta$sample_reads
unique(df_meta$bar_type)
unique(df_meta$bar_type_2)

stop()
################################################################################
message('get sample count for figures')
################################################################################
names(df_meta)
unique(df_meta$figure)
df_meta2=df_meta %>% 
  select(sample_id_PPR,sample_name,sample_id,figure,`George Figure number`,sample_type) %>% 
  filter(!is.na(figure)) %>% 
  filter(str_detect(figure, paste(c("2","3","4","5"), collapse = "|"))) %>% 
  group_by(figure) %>%
  summarise(n = n())
df_meta3=df_meta %>% 
  select(sample_id_PPR,sample_name,sample_id,figure,`George Figure number`,sample_type) %>% 
  filter(!is.na(figure)) %>% 
  filter(str_detect(figure, paste(c("2","3","4","5"), collapse = "|"))) %>% 
  group_by(figure,sample_type) %>%
  summarise(n = n())
df_meta3

################################################################################
################################################################################

df0<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "MASTER") 
df1<- read_excel("data_tables/RUSH_Fig 2_to_Fig 6 Sample Tables.xlsx",sheet = "Figure 6") 

names(df0)

names(df1)

df2 <- df1 %>% 
  mutate(sample_name=`GMCF Number`) %>% 
  mutate(sample_name = str_replace_all(sample_name,"-","_")) %>% 
  mutate(value_bacterial=`Bacterial Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`)%>%
  mutate(value_fungal=`Fungall Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`)%>%
  rename(treatment=Treatment,location=Location) %>% 
  #  mutate(sample_name = sub("-[^-]*$", "", sample_name))%>%
  mutate(experiment='Clipper') %>% 
  select(sample_name,`GMCF Number`,value_bacterial,`Bacterial Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`,
         value_fungal,`Fungall Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`,treatment,location,experiment,everything()) %>% 
  ungroup()

names(df2)

meta=df2
meta_names=names(meta)

write.csv(meta,'data_tables/meta.csv',row.names=FALSE)

