# Load global settings and helper functions
source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)), '/_globalStuff.R'))
################################################################################
# Script Title & Output Folders
################################################################################
script_title <- 'abundance'


# for(plot_set in plot_set_order){

# folder=plot_set


# output_data <- paste0('output_data/', script_title, '/')
# output_plot <- paste0('output_plot/', script_title, '/')
# output_data <- paste0('output_data/',folder,'/')
# output_plot <- paste0('output_plot/',folder,'/')


# if(plot_set=='standard'){
#   normalization_order=c('fixed cycles')
#   data_class_order=c('doubleton')
#   normalization_name=paste(normalization_order,collapse='_')
#   data_class_name=paste(data_class_order,collapse='_')
# 
# }
# if(plot_set=='normalization'){
#   normalization_order=c('fixed cycles','targeted fluorescence')
#   data_class_order=c('doubleton')
#   normalization_name=paste(normalization_order,collapse='_')
#   data_class_name=paste(data_class_order,collapse='_')
# 
# }
# if(plot_set=='data_class'){
#   normalization_order=c('targeted fluorescence')
#   data_class_order=c('singleton','doubleton')
#   normalization_name=paste(normalization_order,collapse='_')
#   data_class_name=paste(data_class_order,collapse='_')
# 
# }

################################################################################
# Dataset selection
################################################################################
zoom=3
starting_r <- 1  # Control starting dataset index
ending_r <- 1    # Control ending dataset index
testing='yes';testing_r=3 # run single dataset index
selected_data_class='doubleton' # adding filter by
selected_data_class_name=paste(selected_data_class,collapse='_')# modified p_title # 122
################################################################################
# Data processing settings
################################################################################
#make_new_data_files <- 'yes'
use_custom_labels <- 'no' # can customize strip labels
################################################################################
# Variable sets for analysis
################################################################################

#custom_name=paste0('_',folder)# start with underscore

figure_number='4'

figure_order=c('4','4A','4B','4C')
parameter_sets <- list(
  set1 = list(
    # filter1_group = 'figure',
    # y_axis_group = 'experiment_type',
    # x_facet_group = 'sampling_device',
    # y_facet_group = 'experiment',
    # number_size = 3
    filter1_group = 'figure',
    y_axis_group = 'sampling_device',
    x_facet_group = 'figure',
    y_facet_group = 'experiment',
    number_size = 3
    #   ),
    # set2 = list(
    #   filter1_group = 'data_class',
    #   y_axis_group = 'normalization',
    #   x_facet_group = 'temperature',
    #   y_facet_group = 'sample_type',
    #   number_size = 1
    # ),
    # set3 = list(
    #   filter1_group = 'normalization',
    #   y_axis_group = 'data_class',
    #   x_facet_group = 'temperature',
    #   y_facet_group = 'sample_type',
    #   number_size = 1
  )
)


############################################################################
# optional Run_Group

# Run_Group='primer_set' # will make separate runs for each value in subgroup_order
# Run_Group_order='Standard'
############################################################################
################################################################################
# imported functions used
################################################################################

################################################################################
message(' 102 # Main Loop: Process Datasets')
################################################################################
Loop_Group_order=Loop_Group='run_once'
if(exists('Run_Group')){Loop_Group=Run_Group;Loop_Group_order=Run_Group_order}
for(lp in Loop_Group_order){
  
  if(!exists('ending_r')){ending_r=nrow(matrix_names);message('\tsetting ending r to nrows matrix names')}
  if(!exists('starting_r')){starting_r=1;message('\tsetting starting r to 1')}
  if(testing=='yes'){ r=starting_r=ending_r=testing_r;message('\trunning r= ',testing_r)}
  r=1 # for manual testing
  
  for (r in starting_r:ending_r) {
    # Extract dataset details
    taxa_levs <- matrix_names[r, 'taxa_levs']
    data_set <- matrix_names[r, 'data_sets']
    taxa_plural <- matrix_names[r, 'taxa_plural']
    
    # # Adjusting for full loop
    # data_set=plot_set
    # data_set_order=plot_set_order
    
    # Assign analysis parameters dynamically
    list2env(parameter_sets[[paste0('set', match(data_set, data_set_order))]], envir = .GlobalEnv)
    groups <- c('filter1_group', 'y_axis_group', 'x_facet_group','y_facet_group')
    for (var in groups) {
      group_value <- get(var)  # e.g., "data_class"
      assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
    }
    grouping_columns=c(y_axis_group, x_facet_group,y_facet_group)
    test_across = c(x_facet_group,y_facet_group)
    # group_within=c(y_axis_group,y_facet_group)
    
    run_suffix=''
    if(exists('Run_Group')){run_suffix=paste0('_',lp)}
    
    p_title <- paste0(taxa_levs, '_', script_title,custom_name,'_',data_set,run_suffix)
    
    qPrint(p_title)
    
    ################################################################################
    message (' 124 # Load and Filter Data')
    ################################################################################
    df_matrix1=read.csv(matrix_names[r,'file_path'],check.names=FALSE,row.names=1)%>%
      t()%>%as.data.frame()%>%
      xPlode_sample_name()%>%
      filter(.data[[filter1_group]] %in% filter1_group_order)%>%
      filter(.data[[y_facet_group]] %in% y_facet_group_order)%>%
      filter(.data[[x_facet_group]] %in% x_facet_group_order)%>%
      filter(.data[[y_axis_group]] %in% y_axis_group_order) %>% 
      filter(str_detect(Figure, figure_number))
    
    
    if(exists('Run_Group')){
      df_matrix1=df_matrix1%>%filter(.data[[Loop_Group]] %in% lp)
    }
    
    qPrint(p_title)
    
    df_matrix2=df_matrix1%>%
      imPlode_sample_name()%>%
      mutate_all(as.numeric)%>%
      select(which(colSums(.) > 0))
    rowSums(df_matrix2)
    
    
    df_matrix=combine_low_abundance(df=df_matrix2,threshold = 0.01,name='(Consolidated < 1%)')
    df_matrix=combine_other(df=df_matrix2, patterns = atcc_genera,starts_with=TRUE)
    names(df_matrix)
    
    feature_names=names(df_matrix)
    ################################################################################
    # Statistical Analysis
    ################################################################################
    
    #kw_group_results <- Kruskal_Mann_Whitney_Test(df = df_matrix, testing_group = y_axis_group, category_group = test_across)
    
    
    ################################################################################
    message (' 137   # Data Processing and Normalization')
    ################################################################################
    # get feature labels to 
    feature_labels_df=feature_labels(y_facet_group,averaging_group = x_facet_group)
    
    y_axis_group_order=rev(y_axis_group_order)
    
    plot_ind_df <- df_matrix %>%
      xPlode_sample_name()%>%
      pivot_longer(cols = any_of(feature_names), names_to = 'feature', values_to = 'counts') %>%
      left_join(feature_labels_df %>% select(-c(taxa)), by = c('feature', y_facet_group)) %>%
      # left_join(feature_labels_df %>% select(feature, taxa)) %>%
      # mutate(Domain = sapply(taxa, extract_domain)) %>%
      ##########################################################################
    mutate(
      label_new = ifelse(
        feature_label %in% atcc_genera,
        
        paste0("<i>", feature_label, "</i> "),
        feature_label
      )
    ) %>% 
     # filter(label_new!='Other') %>% 
    ############################################################################
      group_by(feature, across(any_of(y_facet_group))) %>%
      mutate(grouped_feature_counts = sum(counts))%>%
      filter(grouped_feature_counts > 0)%>%
      group_by(sample_name) %>%
      mutate(sample_counts = sum(counts)) %>%
      ungroup() %>%
      mutate(rel_abun = counts / sample_counts) %>% 
      mutate(percent_rel_abun = rel_abun/.01) %>% 
      mutate(log10_percent_rel_abun = log10(percent_rel_abun)) %>% 
      group_by(across(any_of(y_facet_group))) %>%
      mutate(group_counts = mean(sample_counts)) %>%
      ungroup() %>%
      mutate(rel_abun_log10 = log10(rel_abun+1e-9)) %>% 
      mutate(group_counts_log10 = log10(group_counts+1e-9)) %>% 
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group))))%>%
      mutate(!!sym(y_facet_group):=factor(!!sym(y_facet_group), levels = y_facet_group_order))%>%
      mutate(!!sym(y_axis_group):=factor(!!sym(y_axis_group), levels = y_axis_group_order))%>%
      mutate(!!sym(x_facet_group):=factor(!!sym(x_facet_group), levels = x_facet_group_order))%>%
      mutate(new_plot_group=ifelse(figure=='4A','4A','4BC')) %>% 
      arrange(plot_order) %>% 
      mutate(label_new=factor(label_new,levels=rev(c('Other',setdiff(unique(label_new),'Other'))))) %>% 
      ################################################################################
    filter(label_new!='Other') %>% 
    ################################################################################
      ungroup()
    
    
    
    plot_df=plot_ind_df%>%
      #group_by(feature, Domain, across(any_of(grouping_columns))) %>%
      group_by(feature,label_new, across(any_of(grouping_columns))) %>%
      summarize(group_mean_counts = mean(counts),
                sample_counts = mean(sample_counts),
                mean_log10_percent_rel_abun = mean(log10_percent_rel_abun),
                mean_percent_rel_abun = mean(percent_rel_abun),
                group_counts = mean(group_counts),
                .groups = 'drop') %>%
      #group_by(across(any_of(y_facet_group))) %>%
      #mutate(group_counts = mean(sample_counts)) %>%
      ungroup() %>%
      #  mutate(xmin_log10=log10(group_mean_counts-sd_rel_abun+1e-9)) %>% 
      # mutate(xmax_log10=log10(group_mean_counts+sd_rel_abun+1e-9)) %>% 
      # mutate(mean_rel_abun_log10 = log10(mean_rel_abun+1e-9)) %>% 
      # mutate(group_counts_log10 = log10(group_counts+1e-9)) %>% 
      # mutate(sd_rel_abun_log10 = log10(sd_rel_abun+1e-9)) %>% 
      left_join(feature_labels_df %>% select(-c(taxa)), by = c('feature', y_facet_group)) %>%
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group))))%>%
      mutate(!!sym(y_facet_group):=factor(!!sym(y_facet_group), levels = y_facet_group_order))%>%
      mutate(!!sym(y_axis_group):=factor(!!sym(y_axis_group), levels = y_axis_group_order))%>%
      mutate(!!sym(x_facet_group):=factor(!!sym(x_facet_group), levels = x_facet_group_order))%>%
      ################################################################################
    # mutate(
    #   label_new = ifelse(
    #     feature_label %in% c(
    #       "(Consolidated < 1%)", "Unclassified", "Unclassified Bacilli",
    #       "Unclassified Bacteria", "Unclassified Paracoccaceae",
    #       "Unclassified Propionibacteriaceae", "Unmapped",'Unclassified Enterobacteriaceae'
    #     ),
    #  #   paste0(feature_label, " ", rel_abun_label),                # plain
    #  #   paste0("<i>", feature_label, "</i> ", rel_abun_label)      # italicized
    #     paste0(feature_label),                # plain
    #     paste0("<i>", feature_label, "</i>")      # italicized
    #   )
    # 
    # ) %>% 
      # mutate(label_new=case_when(
      #   label_new=='Unclassified Paracoccaceae'~'Unclassified<br>Paracoccaceae',
      #   label_new=='Unclassified Enterobacteriaceae'~'Unclassified<br>Enterobacteriaceae',
      #   TRUE~label_new
      # )) %>%
      select(rev(everything())) %>% 
      ################################################################################
    ungroup()
    
    ################################################################################
    message( '164 # Generate Plots')
    ################################################################################
    ################################################################################
    # Plotting functions and variables
    ################################################################################
    margin_size=10
    annotate_text_size=6*magnify
    axis_text_size=16*magnify
    axis_title_text_size=18*magnify
    strip_text_size=20*magnify
    legend_text_size=20*magnify
    
    theme_common<-theme_global(base_size = 11)+
      theme(       
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # axis.text.x = element_markdown(angle = 45, 
        #                            #size = axis_text_size,
        #                            hjust=1,vjust=1),
        legend.position='bottom'
      )
    
    plot_width=10
    plot_height=10
    margin_buffer=30
    
    
    gPlot <- function(p) {
      p=p+
        #   scale_fill_manual(values = palette_color) +
        scale_color_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
        scale_fill_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
        labs(
          fill='',
          title='',
          caption='',
          #   title = taxa_levs,
          #   x = "Mean Reads",
          y = "MMC Constituents",
          #   caption = 'mean differential abunance'  
        )+
        guides(color = "none")+
        #scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
        
        theme_common
      print(p)
    }
    ###############################################################################
    message('207 # strip backgrounds with outlines')
    ################################################################################
    unique_y_labels <- unique(plot_df[[y_facet_group]])
    
    
    outline_color_y='black'
    
    s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),'"]))')
    print(s)
    eval(parse(t=s))
    text_color=c("white", rep("black", length(unique_y_labels) - 1))
    
    
    if (use_custom_labels=='yes'){
      custom_y_labels=c("slope","1500TF","x3FC","slope","1500TF","x3FC",'24cycles',rep('Soil',3),rep('Zymo',3),'NTC')
      
      custom_y_labels
      unique_y_labels=custom_y_labels
      s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),'"]))')
      print(s)
      eval(parse(t=s))
    }
    
    y_strip=as.data.frame(unique_y_labels)%>%
      mutate(text_color=ifelse(unique_y_labels=='Soil','white','black'))%>%
      mutate(text_color='white')%>%
      mutate(face='bold')%>%
      mutate(text_size=strip_text_size)%>%
      ungroup()
    
    ################################################################################
    unique_x_labels <- unique(plot_df[[x_facet_group]])
    
    
    outline_color_x='black'
    
    s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = palette_color["'),'"]))')
    print(s)
    eval(parse(t=s))
    text_color=c("white", rep("black", length(unique_x_labels) - 1))
    
    
    if (use_custom_labels=='yes'){
      custom_x_labels=c("slope","1500TF","x3FC","slope","1500TF","x3FC",'24cycles',rep('Soil',3),rep('Zymo',3),'NTC')
      
      custom_y_labels
      unique_y_labels=custom_y_labels
      s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),'"]))')
      print(s)
      eval(parse(t=s))
    }
    
    x_strip=as.data.frame(unique_x_labels)%>%
      mutate(text_color=ifelse(unique_x_labels=='Soil','white','black'))%>%
      mutate(text_color='white')%>%
      mutate(face='bold')%>%
      mutate(text_size=strip_text_size)%>%
      ungroup()
    
    ###############################################################################
    message('248 # plot variables')

   # sc = abs(floor(min(plot_df$mean_rel_abun_log10))) 
axis_adjust=1
    
    x_min=log10(2);x_alt=1 # modifies mean count text
    magnify=3# controls text size
    ################################################################################
    dodge_pos <- position_dodge(width = 0.8)
    
    unique(plot_df$sampling_device)
palette_label1=palette_label
    
    
    # new_labels=c(
    # 'cotton'='Cotton',
    # 'macrofoam'='Macrofoam',
    # 'Polyester Wipe'='Polyester Wipe',
    # 'SALSA'='SALSA') 
    # 
    # palette_label=c(new_labels,palette_label1,new_labels)
    plot_ind_df_alt=plot_ind_df
    
    message('st_plot')
    ################################################################################
    
    new_labels=c(
      '4A'='D\\) MMC Released from Swab Heads',
    '4B'='E\\) MMC Recovered by Swabs',
    '4C'='F\\) MMC Rec. by Polyester Wipes and SALSA')
    palette_label=c(new_labels,palette_label)
    

    
    st_plot3=ggplot(plot_df, aes(x = mean_percent_rel_abun , y = plot_order)) +
      geom_col(aes(fill = !!sym(y_axis_group)),width=.6,position = dodge_pos,alpha=.7)+
      geom_jitter(
        data = plot_ind_df,
        aes(y = plot_order,x = rel_abun/.01, fill = !!sym(y_axis_group)),
        position = position_jitterdodge(
          dodge.width = 0.8,  # same as your dodge_pos
          jitter.width = .4,   # no x jitter
          jitter.height = 0 # vertical y jitter
        ),
        color = 'black',
        shape=21,size=2
      )+
      # stat_summary(data=    plot_ind_df, # would average the logs
      #   aes(x=percent_rel_abun,fill = !!sym(y_axis_group)),
      #   fun = mean,                      # or median, sum, etc.
      #   geom = "col",
      #   width = 0.3,
      #   position = dodge_pos,
      #   alpha = 1,
      #   color='black'
      # )+
    
      # geom_jitter(
      #   data = plot_ind_df,
      #   aes(y = plot_order,x = rel_abun/.01*sc, fill = !!sym(y_axis_group)),
      #   position = position_jitterdodge(
      #     dodge.width = 0.8,  # same as your dodge_pos
      #     jitter.width = .4,   # no x jitter
      #     jitter.height = 0 # vertical y jitter
      #   ),
      #   color = 'black',
      #   shape=21,size=2
      # )+
  #    scale_y_discrete(breaks = custom_breaks, labels = custom_labels)+
      facet_grid2(formula(paste('~',x_facet_group))
                  ,  labeller = labeller(
                    .rows = as_labeller(palette_label),
                    .cols = as_labeller(palette_label)
                  )
                  ,space = "free")
    
    st_plot3
    gPlot(st_plot3)
    
    magnify=1

    y_labels <- setNames(as.character(plot_df$label_new), plot_df$plot_order)
    
    st_plot2=st_plot3+
      geom_text(
        aes(
          label = round(mean_percent_rel_abun, 1),
          x = ifelse(mean_percent_rel_abun < 50,90, mean_percent_rel_abun),
          color = !!sym(y_axis_group),
          group = !!sym(y_axis_group)
        ),
        size = number_size,
        hjust = 1.1,
        fontface = 'bold',
        show.legend = FALSE,
        position = position_dodge(width = 0.8))+
      scale_y_discrete(labels = y_labels)+ 
      
      theme(axis.text.y = element_markdown('bold',size=axis_text_size))
    st_plot2
    
    st_plot=gPlot(st_plot2)+
      scale_y_discrete(labels = y_labels,position='right')+
      labs(title='',x='Mean Percent Relative Abundance')+
     # scale_x_continuous(trans = log10_decimal(factor = 10)#
                         #,breaks=c(.1,1,10,100)
        #                 ) +
      # stat_count(data=    plot_ind_df,
      #   aes(y = plot_order,
      # 
      #       label = after_stat(paste0("n = ", count))
      #   ),
      #   inherit.aes = FALSE,
      #   geom = "text",
      #   vjust = 0,
      #   fontface = "bold",
      #   size = 3.5
      # )+
      theme(     strip.text = element_markdown(size = rel(1.2), face = 'bold'))+
      guides(color = guide_legend(reverse = TRUE),
             fill  = guide_legend(reverse = TRUE)
      )
    
    st_plot
    
    saveRDS(st_plot,paste0(output_plot,'Figure_',figure_number,'_st_plot.rds'))
    ggsave(paste0(output_plot,'Figure_',figure_number,'_st_plot.png'),st_plot)
    
    stop()
    
    
    
    
    
    
    
    
    
plot_df_alt=plot_df
unique(plot_df_alt$sampling_device)
    ################################################################################
plot_ind_df=plot_ind_df_alt %>% filter(sampling_device %in% c('cotton','macrofoam'))
plot_df=plot_df_alt %>% filter(sampling_device %in% c('cotton','macrofoam'))
plot_ind_df=plot_ind_df_alt %>% filter(new_plot_group %in% c('4A'))
plot_df=plot_df_alt %>% filter(new_plot_group %in% c('4A'))


    message('st_plot1')
    ################################################################################
    st_plot1=ggplot(plot_df, aes(x = mean_rel_abun/.01, y = plot_order)) +
      geom_col(aes(fill = !!sym(y_axis_group)),width=.6,position = dodge_pos,alpha=.7)+
      geom_jitter(
        data = plot_ind_df,
        aes(y = plot_order,x = rel_abun/.01, fill = !!sym(y_axis_group)),
        position = position_jitterdodge(
          dodge.width = 0.8,  # same as your dodge_pos
          jitter.width = .4,   # no x jitter
          jitter.height = 0 # vertical y jitter
        ),
        color = 'black',
        shape=21,size=2
      )+
      scale_y_discrete(breaks = custom_breaks, labels = custom_labels)+
      facet_grid2(formula(paste('~',x_facet_group))
                  ,  labeller = labeller(
                    .rows = as_labeller(palette_label),
                    .cols = as_labeller(palette_label)
                  )
                  ,space = "free")
    
    gPlot(st_plot1)

    magnify=1
    
    y_labels <- setNames(plot_df$label_new, plot_df$plot_order)
    
    st_plot1=st_plot1+
      geom_text(
        aes(
             label = round(mean_rel_abun/.01, 1),
          x = ifelse(mean_rel_abun/.01 < 50,90, mean_rel_abun/.01),
          color = !!sym(y_axis_group),
          group = !!sym(y_axis_group)
        ),
        size = number_size,
        hjust = 1.1,
        fontface = 'bold',
        show.legend = FALSE,
        position = position_dodge(width = 0.8))+
      scale_y_discrete(labels = y_labels)+ 

      theme(
        axis.text.y = element_markdown('bold',size=axis_text_size))

    st_plot1=gPlot(st_plot1)+
      scale_y_discrete(labels = y_labels,position='right')+
      labs(title='',x='Percent Relative Abundance')
    saveRDS(st_plot1,paste0(output_plot,'Figure_',figure_number,'_st_plot1.rds'))
    ggsave(paste0(output_plot,'Figure_',figure_number,'_st_plot1.png'),st_plot1)
    ################################################################################
    plot_ind_df=plot_ind_df_alt %>% filter(sampling_device %in% c('Polyester Wipe','SALSA'))
    plot_df=plot_df_alt %>% filter(sampling_device %in% c('Polyester Wipe','SALSA'))
    plot_ind_df=plot_ind_df_alt %>% filter(new_plot_group %in% c('4BC'))
    plot_df=plot_df_alt %>% filter(new_plot_group %in% c('4BC'))
    message('st_plot2')
    ################################################################################
    st_plot2=ggplot(plot_df, aes(x = mean_rel_abun/.01, y = plot_order)) +
      geom_col(aes(fill = !!sym(y_axis_group)),width=.6,position = dodge_pos,alpha=.7)+
      geom_jitter(
        data = plot_ind_df,
        aes(y = plot_order,x = rel_abun/.01, fill = !!sym(y_axis_group)),
        position = position_jitterdodge(
          dodge.width = 0.8,  # same as your dodge_pos
          jitter.width = .4,   # no x jitter
          jitter.height = 0 # vertical y jitter
        ),
        color = 'black',
        shape=21,size=2
      )+
      scale_x_continuous(trans = log10_factor(100))+
      scale_y_discrete(breaks = custom_breaks, labels = custom_labels)+
      facet_grid2(formula(paste('~',x_facet_group))
                  ,  labeller = labeller(
                    .rows = as_labeller(palette_label),
                    .cols = as_labeller(palette_label)
                  )
                  ,space = "free")
    
    gPlot(st_plot2)
    
    magnify=1
    
    y_labels <- setNames(plot_df$label_new, plot_df$plot_order)
    
    st_plot2=st_plot2+
      geom_text(
        aes(
          label = round(mean_rel_abun/.01, 1),
          x = ifelse(mean_rel_abun/.01 < 50,90, mean_rel_abun/.01),
          color = !!sym(y_axis_group),
          group = !!sym(y_axis_group)
        ),
        size = number_size,
        hjust = 1.1,
        fontface = 'bold',
        show.legend = FALSE,
        position = position_dodge(width = 0.8))+
      scale_y_discrete(labels = y_labels)+ 
      
      theme(
        axis.text.y = element_markdown('bold',size=axis_text_size))
    
    st_plot2=gPlot(st_plot2)+
      scale_y_discrete(labels = y_labels,position='right')+
      labs(title='',x='Percent Relative Abundance')
    saveRDS(st_plot2,paste0(output_plot,'Figure_',figure_number,'_st_plot2.rds'))
    ggsave(paste0(output_plot,'Figure_',figure_number,'_st_plot2.png'),st_plot2)
    ################################################################################
    ################################################################################
    library(cowplot)
    
    # Extract the legend from plot_A_13
    legend_st_plot1 <- get_legend(st_plot1 + theme(legend.position = "bottom"))
    
    # Wrap as a plot
    legend_plot <- ggdraw(legend_st_plot1)
    
    legend_plot 
    ggsave(paste0(output_plot,'Figure_',figure_number,'_legend_plot.png'),legend_plot)
    saveRDS(legend_plot,paste0(output_plot,'Figure_',figure_number,'_legend_plot.rds'))
    
    
    ################################################################################
    
    ################################################################################
    library(cowplot)
    
    # Extract the legend from plot_A_13
    legend_st_plot1 <- get_legend(st_plot1 + theme(legend.position = "right"))
    
    # Wrap as a plot
    legend_plot <- ggdraw(legend_st_plot1)
    
    legend_plot 
    ggsave(paste0(output_plot,'Figure_',figure_number,'_legend_plot_right.png'),legend_plot)
    saveRDS(legend_plot,paste0(output_plot,'Figure_',figure_number,'_legend_plot_right.rds'))
    
    
    ################################################################################

    
  }
}
