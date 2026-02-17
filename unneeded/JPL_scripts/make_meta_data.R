source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
################################################################################
# Description
script_title='meta_data'

################################################################################
# Read in file names from Phylum table

# create meta.csv (meta data) for other scripts


################################################################################
sample_id <- read.csv(paste0('data_tables/original/Merged_Phylum_raw_abund.csv'),check.names = FALSE,row.names = 1)%>%
  column_to_rownames(var='ID') %>% 
  t%>%as.data.frame() %>% 
  rownames_to_column(var='sample_id')%>%
  filter(str_detect(sample_id,'Illumina')) %>% 
  mutate(sample_name =  str_remove(sample_id, "_Illumina")) %>% 
  select(sample_id,sample_name)

################################################################################
library(readxl)

qLoad('data_tables/sample_reads.csv')

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
  left_join(sample_id) %>% 
#  mutate(sample_name = sub("-[^-]*$", "", sample_name))%>%
  mutate(experiment='Clipper') %>% 
  select(sample_name,`GMCF Number`,value_bacterial,`Bacterial Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`,
         value_fungal,`Fungall Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`,treatment,location,experiment,everything()) %>% 
ungroup()

names(df2)

meta=df2
meta_names=names(meta)

write.csv(meta,'data_tables/meta.csv',row.names=FALSE)

