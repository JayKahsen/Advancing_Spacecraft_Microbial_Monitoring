
# load file for common project stuff
source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
################################################################################
# Description
script_title='rarefied_tables'
################################################################################
# Read in raw tables from QIIME2 output
# create Matrix.names.csv for other scripts
# Ready tables for use in other scripts
#Drop low sample counts, Rarefy, filter to minimum counts or relative abundance.
################################################################################
#library(scales)
################################################################################
# IMPORTANT DATA
################################################################################
if(!exists('data_sets')){data_sets='s0';s0='s0'}
if(data_sets=='s0'){s0='s0'}
qPrint(data_sets)
for(data_set in data_sets){print(get(sym(data_set)))}




################################################################################
# Switches 
################################################################################
use_minimum_count_group=making_filtered_matrix=making_matrix_names=making_raw_matrix=making_rarefied_matrix=look_at_counts=testing='no'

# 54 create matrix_names.csv to load files for other scripts
# making_matrix_names='yes'

# 79 format QIIME2 tables for standard operations
 making_raw_matrix='yes'     # make sure names match with meta

qPrint(data_sets)
for(data_set in data_sets){print(get(sym(data_set)))}
# determine min sample count
  look_at_counts='yes'

# 213 separate data by primer sets;drop low count samples and rarefy by repeated subsampling for each Type
#
# filtering rarefied matrix by feature abundance across samples by Type
#making_filtered_matrix='yes'

#testing='yes'


################################################################################
# Variables for script
################################################################################
# define minimum number of acceptable counts

filter_by_counts=0
min_relative_abundance=1/100000

raw_matrix_path='data_tables/raw_matrix'
rarefied_matrix_path='data_tables/rarefied_matrix'
filtered_matrix_path='data_tables/filtered_matrix'

################################################################################
# make matrix names
# - a file to help work on different data sets
################################################################################

if (making_matrix_names == 'yes') { qPrint(making_matrix_names)
  
  create_directory(raw_matrix_path)
  create_directory(rarefied_matrix_path)
  create_directory(filtered_matrix_path)
  
  names(taxa_plural) <- taxa_levels
  
  matrix_names <- expand.grid(taxa_levs = taxa_levels, data_sets = data_sets) %>%
    mutate(taxa_plural = case_when(
      taxa_levs %in% names(taxa_plural) ~ taxa_plural[taxa_levs],
      TRUE ~ taxa_levs))%>%
    mutate(raw_path = paste0(raw_matrix_path,'/raw_matrix_',taxa_levs,'.csv'))%>%
    mutate(rarefied_path = paste0(rarefied_matrix_path,'/rarefied_matrix_', taxa_levs, '_', data_sets, '.csv'))%>%
    mutate(filtered_path = paste0(filtered_matrix_path,'/filtered_matrix_', taxa_levs, '_', data_sets, '.csv'))%>%
    mutate(file_path = filtered_path)%>%
    arrange(taxa_levs)
  
  write.csv(matrix_names,'data_tables/matrix_names.csv',row.names = FALSE)
  print('created matrix_names.csv')
  message()   
  stop('\tdeselect making matrix names to continue') 
}

################################################################################
# make raw matrix
# - create a matrix for each data set in a standard form
################################################################################
if(making_raw_matrix=='yes'){qPrint(making_raw_matrix)
  # meta<<-read.csv('data_tables/meta.csv')
  # meta_names=names(meta)
  sample_reads_df2=NULL
  
  for (taxa_levs in taxa_levels){qPrint(taxa_levs)
    
    if(taxa_levs=='ASV'){next()
      raw_matrix <- read.delim(paste0('data_tables/original/Table_Features_PostRemoval.tsv'),check.names = FALSE,row.names = 1, skip = 1)%>%
        t%>%as.data.frame()
      feature_names=names(raw_matrix)
    }else if(taxa_levs=='Contaminated'){
      raw_matrix <- read.delim(paste0('data_tables/original/Table_Features_PreRemoval.tsv'),check.names = FALSE,row.names = 1, skip = 1)%>%
        t%>%as.data.frame()
      feature_names=names(raw_matrix)
    }else{
      file=(paste0('data_tables/original/PPR_XT/Project_1790_PPR_XT_Output.',tolower(taxa_levs),'.allfilter.abund.tsv'))
      file2=(paste0('data_tables/original/PPR_XT_2/Project_1790_PPR_XT_2_Output.',tolower(taxa_levs),'.allfilter.abund.tsv'))
    }
    df_matrix_XT<- read.delim(file,check.names = FALSE,row.names = 1) %>% 
      t() %>% as.data.frame() %>% 
      rownames_to_column(var = 'sample_id_XT') %>%
      left_join(meta)%>%
      select(sample_name,everything()) %>% 
   filter(!is.na(sample_name)) %>% 
      imPlode_sample_name()
    sum(is.na(df_matrix_XT))
    
    df_matrix_XT_2<- read.delim(file2,check.names = FALSE,row.names = 1)%>% 
      t() %>% as.data.frame() %>% 
      rownames_to_column(var = 'sample_id_XT_2') %>%
      left_join(meta) %>%
      select(sample_name,everything()) %>% 
      filter(!is.na(sample_name)) %>% 
      imPlode_sample_name()
    sum(is.na(df_matrix_XT_2))
    
    df_matrix2 <- bind_rows(df_matrix_XT, df_matrix_XT_2) %>% 
      replace(is.na(.), 0)
    

    # fix sample_name 
    sum(is.na(df_matrix2))
    #df_matrix2[is.na(df_matrix2)] <- 0

    df_matrix1 <- df_matrix2 
    
    
    sample_reads=rowSums(df_matrix1)%>%
      as.data.frame()%>%
      setNames(paste("Reads"))%>%
      mutate(taxa_level=paste0(taxa_levs,'_Reads'))%>%
      xPlode_sample_name()
    
    sample_reads_df2=rbind(sample_reads_df2,sample_reads)
    
    # transpose for quicker reading
    df_matrix1%>%
      t()%>%as.data.frame()%>%
      write.csv(paste0(raw_matrix_path,'/raw_matrix_',taxa_levs,'.csv'))
  }
  print('created raw matrixes')
  
  sample_reads_df=sample_reads_df2%>%
    pivot_wider(names_from = taxa_level,values_from = Reads)%>%
    imPlode_sample_name()
  
  write.csv(sample_reads_df,paste0('data_tables/sample_reads.csv'))
  print('created raw matrices.csv')
  stop('\tdeselect making raw matrix to continue') 
} # skipping if making_raw_matrix != yes

##############################################################################
# look at counts
##############################################################################
###############################################################################
# visualizing cut off effects
##############################################################################
if(look_at_counts=='yes'){qPrint(look_at_counts)
  taxa_levs='Phylum'
  matrix_df=read.csv(paste0(raw_matrix_path,'/raw_matrix_',taxa_levs,'.csv'),row.names=1,check.names=FALSE)%>%t%>%as.data.frame()
  
  feature_names=names(matrix_df)
  matrix_df$sample_counts=rowSums(matrix_df)
  
  df=matrix_df%>%
    select(sample_counts)%>%
    xPlode_sample_name()%>%
    # filter(sample_type %in% setdiff(sample_type_order,'NTC'))%>%
    #mutate(primer_set=factor(primer_set,levels=primer_set_order))%>%
    #arrange(normalization,mastermix,primer_set)%>%
    ungroup()%>%
    filter(sample_type!='NTC')%>%
    filter(experiment %in% c('G4','G6','G7'))%>%
    mutate(short_name=factor(short_name, levels = unique(short_name)))
  
  ggplot(df,aes(x=short_name,y=sample_counts,color=temperature,shape=sample_type))+
   # geom_point(aes(color=normalization),size=5,alpha=.2)+
    geom_point()+
    #scale_y_log10() + 
   # geom_hline(yintercept = 1000, linetype = "dashed", color = "red")+
facet_wrap(~experiment,scales='free')+
    theme(axis.text.x = element_text(angle = 90))
  
 
  ggsave(paste(output_plot,'looking at raw counts',project_folder,'.png'),width=36,height=12)
  
  
  # move to making rarefied
  #minimum_counts <- c(D4 = 150000, D5 = 100000, D6 = 250000, D7 = 250000, E1 = 50000, E2 = 300000)
  message("choose minimum counts and deselect looking to continue")
  
  stop()}
# stop to look at counts and set minimum)
##############################################################################
# subsample matrixes
if (making_rarefied_matrix == 'yes') {qPrint(making_rarefied_matrix)
  ##############################################################################
  ##############################################################################
  # Function
  ##############################################################################
  # auto_min_counts max of slice1
  auto_min_counts=function(sample_count_dfx,keep_number){
    min_count=sample_count_dfx%>%
      group_by(short_name)%>%
      arrange(desc(sample_count))%>%
      slice(1:keep_number)%>%
      ungroup()%>%
      pull(sample_count)%>%
      min()
    return(min_count)    
  }
  # rarefy function
  rarefy_function=function(dfx,min_sample_count){
    df=dfx%>%
      filter(sample_count >= min_sample_count) %>%
      select(-sample_count)%>%
      imPlode_sample_name()%>%
      as.matrix()
    set.seed(20250314) 
    rarefied_df=rrarefy(df,sample=min_sample_count)%>%
      as.data.frame()
    
    return(rarefied_df)    
  }
  ##############################################################################
  # prepping
  ##############################################################################
  minimum_counts=10000
  
  matrix_names=read.csv('data_tables/matrix_names.csv',check.names = FALSE)
  r=2
  n_rows=1
  if(testing=='no'){n_rows=nrow(matrix_names)}
  
  # load matrix from names csv
  for (r in 1:n_rows){
    
    taxa_levs=matrix_names[r,'taxa_levs']
    data_set=matrix_names[r,'data_sets']
    p_title=paste0(script_title,'_',taxa_levs,'_',data_set)
    
    data_set_group <- get(sym(data_set))
    qPrint(data_set)
    qPrint(data_set_group)
    
    matrix_df=read.csv(matrix_names[r,'raw_path'],check.names=FALSE,row.names=1)%>%
      t()%>%as.data.frame()
    
    feature_names=names(matrix_df)
    
    # matrix_df1=matrix_df;matrix_df1$sample_count=rowSums(matrix_df)
    # 
    # data_set_df=matrix_df1%>%
    #   xPlode_sample_name()%>%
    #   filter(.data[[data_set_filter_group]] %in% data_set_group)%>%
    #   ungroup()
    
    ##############################################################################
    # Filtering for data sets and selected group
    ##############################################################################
    data_set_df <- matrix_df %>%
      mutate(sample_count = rowSums(.)) %>%
      xPlode_sample_name() %>%
      filter(.data[[data_set_filter_group]] %in% data_set_group) %>%
      filter(.data[[rarefy_by_group]] %in% selected_groups) %>%
      ungroup()
    
    ##############################################################################
    # looping by rarefy_by_group if exists while rarefying
    ##############################################################################
    rarefied_df=NULL
    #rarefied_df <- data.frame()
    if (exists("rarefy_by_group")) {
      # If `rarefy_by_group` exists, loop over unique values in the column
      for (t in unique(data_set_df[[rarefy_by_group]])) {
        qPrint(t)
        
        # Subset data for the current group
        subset_df <- data_set_df %>%
          filter(.data[[rarefy_by_group]] == t)
        
        # Determine minimum sample count
        min_sample_count <- minimum_counts
        if (use_auto_minimum_count == 'yes') {
          min_sample_count <- auto_min_counts(subset_df, keep_number)
        }
        
        # Rarefy the subset
        rarefied_subset <- rarefy_function(subset_df, min_sample_count)
        
        # Bind the rarefied subset to the result data frame
        rarefied_df <- bind_rows(rarefied_df, rarefied_subset)
      }
      print(rarefied_df%>%
              mutate(sample_count = rowSums(.)) %>%
              xPlode_sample_name()%>%
              select(any_of(rarefy_by_group),sample_count)%>%
              distinct())
      
      
    } else {
      # If `rarefy_by_group` does not exist, no subsetting
      message('\t No subsetting')
      
      # Determine minimum sample count for the entire dataset
      min_sample_count <- minimum_counts
      if (use_auto_minimum_count == 'yes') {
        min_sample_count <- auto_min_counts(data_set_df, keep_number)
      }
      
      # Rarefy the entire dataset
      rarefied_df <- rarefy_function(data_set_df, min_sample_count)
      
      print(table(rowSums(rarefied_df)))  
    }
    
    
    print(paste('writing',Sys.time(),'writing'));cat('\n')
    rarefied_df%>%
      t()%>%as.data.frame()%>%
      write.csv(paste0(rarefied_matrix_path,'/rarefied_matrix_',taxa_levs,'_',data_set,'.csv'))
    
  }
} # skipping if making_rarefied_matrix != yes
##############################################################################
if (making_filtered_matrix == 'yes') {qPrint(making_filtered_matrix)
  ##############################################################################
  matrix_names=read.csv('data_tables/matrix_names.csv',check.names = FALSE)
  r=1
  n_rows=1
  if(testing=='no'){n_rows=nrow(matrix_names)}
  
  # load matrix from names csv
  for (r in 1:n_rows){
    
    taxa_levs=matrix_names[r,'taxa_levs']
    data_set=matrix_names[r,'data_sets']
    p_title=paste0(script_title,'_',taxa_levs,'_',data_set)
    
    qPrint(p_title)
    
    matrix_df=read.csv(matrix_names[r,'rarefied_path'],check.names=FALSE,row.names=1)%>%
      t()%>%as.data.frame()
    
    feature_names=names(matrix_df)
    
    df=matrix_df%>%
      select(which(colSums(.) > filter_by_counts))%>%
      t()%>%as.data.frame()
    
    write.csv(df,paste0(filtered_matrix_path,'/filtered_matrix_',taxa_levs,'_',data_set,'.csv'))
    
  }
} # skipping if making_filtered_matrix != yes

