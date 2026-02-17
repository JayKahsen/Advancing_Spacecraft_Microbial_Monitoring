# Load global settings and helper functions
source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)), '/_globalStuff.R'))
################################################################################
# Script Title & Output Folders
################################################################################
script_title <- 'alpha_diversity'

################################################################################
# all metrics used
################################################################################
metric_order=c('unique_features','chimeras','percent_Muri','percent_Lachno',
               'percent_Propi','percent_Archaea','percent_good','ideal_score',
               'richness', 'evenness', 'shannon','n')

selected_metric=metric_order


################################################################################

if(!exists('plot_set_order')){plot_set_order='none'}
for(plot_set in plot_set_order){
  
  if(plot_set!='none'){
    
    folder=plot_set 
    output_data <- paste0('output_data/', script_title, '/')
    output_plot <- paste0('output_plot/', script_title, '/')
    output_data <- paste0('output_data/',folder,'/')
    output_plot <- paste0('output_plot/',folder,'/')
  }
  
  if(plot_set=='standard'){
    normalization_order=c('fixed cycles')
    data_class_order=c('doubleton')
    normalization_name=paste(normalization_order,collapse='_')
    data_class_name=paste(data_class_order,collapse='_')
    #metric_filter=setdiff(metric_order,c())
    selected_metric=c('percent_Muri','percent_Lachno',
                       'percent_Propi','percent_Archaea',
                       'richness', 'evenness', 'shannon','n')
  }
  if(plot_set=='normalization'){
    normalization_order=c('fixed cycles','targeted fluorescence')
    data_class_order=c('doubleton')
    normalization_name=paste(normalization_order,collapse='_')
    data_class_name=paste(data_class_order,collapse='_')
   # metric_filter=setdiff(metric_order,c())
    selected_metric=c('percent_Muri','percent_Lachno',
                       'percent_Propi','percent_Archaea',
                       'richness', 'evenness', 'shannon','n')
  }
  if(plot_set=='data_class'){
    normalization_order=c('targeted fluorescence')
    data_class_order=c('singleton','doubleton')
    normalization_name=paste(normalization_order,collapse='_')
    data_class_name=paste(data_class_order,collapse='_')
   # metric_filter=setdiff(metric_order,c())
    selected_metric=c('unique_features','percent_Muri','percent_Lachno',
                       'percent_Propi','percent_Archaea',
                       'richness', 'evenness', 'shannon','n')
    selected_metric=c('percent_Muri','percent_Lachno',
                      'percent_Propi','percent_Archaea',
                      'richness', 'evenness', 'shannon','n')
  }
  
  ################################################################################
  # Dataset selection
  ################################################################################
  testing=make_new_data_files <- 'no'
  make_new_data_files <- 'yes'
  starting_r <- 3  # Control starting dataset index
  ending_r <- 3    # Control ending dataset index
  testing='yes';testing_r=4# run single dataset index
  
  ################################################################################
  # Data processing settings
  ################################################################################
  #make_new_data_files <- 'yes'
  use_custom_labels <- 'no' # can customize strip labels
  ################################################################################
  # Variable sets for analysis
  ################################################################################
  #custom_name=paste0('_',folder)# start with underscore
  

  figure_number='6'
  parameter_sets <- list(
    set1 = list(
      filter1_group='experiment',
      x_axis_group='treatment', # what are we comparing
      x_facet_group='location',
      plot_group='experiment'
    #   ),
    # set2 = list(
    #   filter1_group='data_class',
    #   x_axis_group = 'normalization',
    #   x_facet_group = 'temperature',
    #   plot_group_group = 'sample_type'
    # ),
    # set3 = list(
    #   filter1_group='normalization',
    #   x_axis_group = 'data_class',
    #   x_facet_group = 'temperature',
    #   plot_group_group = 'sample_type'
    ))
  ############################################################################
  # # optional Run_Group
  # 
  # Run_Group='primer_set' # will make separate runs for each value in subgroup_order
  # Run_Group_order='Standard'
  # #Run_Group_order=primer_set_order
  ############################################################################
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
      
      if(plot_set!='none'){
        # Adjusting for full loop
        data_set=plot_set
        data_set_order=plot_set_order
      }
      
      # Assign analysis parameters dynamically
      list2env(parameter_sets[[paste0('set', match(data_set, data_set_order))]], envir = .GlobalEnv)
      groups <- c('filter1_group', 'x_axis_group', 'x_facet_group','plot_group')
      for (var in groups) {
        group_value <- get(var)  # e.g., "data_class"
        assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
      }
      grouping_columns=c(x_axis_group,x_facet_group,plot_group)
      test_across = c(x_facet_group,plot_group)
      #group_within=c(x_axis_group,plot_group)
      group_within=grouping_columns
      
      run_suffix=''
      if(exists('Run_Group')){run_suffix=paste0('_',lp)}
      
      filter1_names=(paste(filter1_group_order,collapse='_'))
      x_axis_names=(paste(x_axis_group_order,collapse='_'))
      x_facet_names=(paste(x_facet_group_order,collapse='_'))

      unique_names_order=unique(c(filter1_group_order,x_axis_group_order,x_facet_group_order))
      unique_names=paste(unique_names_order,collapse='_')
      
      p_title <- paste0(taxa_levs, '_', script_title,custom_name,'_',filter1_names,'_',x_axis_names,'_',x_facet_names,'_',data_set,run_suffix)
      p_title <- paste0(taxa_levs, '_', script_title,custom_name,'_',unique_names,'_',data_set,run_suffix)
      
      qPrint(p_title)
      qPrint(p_title)
      ################################################################################
      if (make_new_data_files=='yes'){ # save data to load from
        ################################################################################
        message (' 124 # Load and Filter Data')
        ################################################################################
        matrix_df1=read.csv(matrix_names[r,'file_path'],check.names=FALSE,row.names=1)%>%
          t()%>%as.data.frame()%>%
          xPlode_sample_name()%>%
          #filter(data_class %in% selected_data_class)%>%
          #filter(.data[[plot_group]] %in% selected_groups)
          filter(.data[[filter1_group]] %in% filter1_group_order)%>%
          filter(.data[[plot_group]] %in% plot_group_order)%>%
          filter(.data[[x_axis_group]] %in% x_axis_group_order)%>%
          filter(.data[[x_facet_group]] %in% x_facet_group_order)%>% 
         # filter(str_detect(Figure, figure_number)) %>% 
          ungroup()
        
        
    
        if(exists('Run_Group')){matrix_df1=matrix_df1%>%filter(.data[[Loop_Group]] %in% lp)}
        
        matrix_df=matrix_df1%>%
          imPlode_sample_name()%>%
          mutate_all(as.numeric)%>%
          select(which(colSums(.) > 0))
        
        feature_names=names(matrix_df)
        
        ################################################################################
        message('# 134 get percent Archaea and other taxa of interest')
        ################################################################################
        
        # grouping features by test across
        df2=matrix_df%>%
          xPlode_sample_name()%>%
          pivot_longer(col=all_of(feature_names),names_to='feature',values_to='counts')%>%
          group_by(feature,across(any_of(test_across)))%>%
          mutate(feature_grouped_counts=sum(counts))%>%
          filter(feature_grouped_counts>0)%>%
          mutate(taxa=feature)%>%
          group_by(sample_name) %>%
          mutate(sample_counts=sum(counts))%>%
          ungroup()
        
        num_unique_features_df=df2%>%
          select(sample_name,feature,counts,any_of(group_within))%>%
          filter(counts>0)%>%
          group_by(across(all_of(group_within)), feature) %>%
          filter(n_distinct(sample_name) == 1) %>%
          ungroup()%>%
          group_by(across(all_of(group_within)), sample_name) %>%
          summarise(smp_unique_features = n(), .groups = "drop")%>%
          select(sample_name,smp_unique_features)%>%
          distinct()
        
        if(matrix_names[r,"taxa_levs"]=='ASV'){
          
          df2=df2%>%
            rename(ASV=taxa)%>%
            left_join(ASV_taxa%>%select(ASV,taxa,Family))%>%
            ungroup()
        }
        
        Archaea_df <- extract_taxa_percentage(df2, name="Archaea", target="archaea")
        
        Propi_df <- extract_taxa_percentage(df2, name="Propi", target="g__Cutibacterium")
        
        Lachno_df <- extract_taxa_percentage(df2, "Lachno", "f__Lachnospiraceae")
        
        Muri_df <- extract_taxa_percentage(df2, "Muri", "f__Muribaculaceae")
        
        ################################################################################
        print('218 # for Zymo Testing Ideal Score')
        ################################################################################
        if("sample_type" %in% colnames(df2)) {
        if ('Zymo' %in% df2$sample_type){
          Zymo_names=c(
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Listeriaceae;g__Listeria',#
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus',#
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus',#
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus',#
            'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus',
            'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella',#
            'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;__',
            'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas')
          
          Zymo_ASV=ASV_taxa%>%filter(str_replace_all(taxa, " ", "") %in% str_replace_all(Zymo_names, " ", ""))%>%
            pull(ASV)
          feature_lable_df=feature_labels()
          
          Zymo_df=matrix_df%>%select(any_of(Zymo_ASV))%>%
            xPlode_sample_name()%>%
            pivot_longer(cols = any_of(Zymo_ASV),names_to='feature',values_to = 'counts')%>%
            left_join(feature_lable_df)%>%
            group_by(feature_label,sample_name,sample_type)%>%
            summarize(counts=sum(counts))
          
          Zymo_labels=unique(Zymo_df$feature_label)
          
          Zymo_scores=Zymo_df%>%
            pivot_wider(names_from = feature_label,values_from = counts)%>%
            mutate(smp_ideal_score = if_else(sample_type=='Zymo', qIdealScore(c_across(any_of(Zymo_labels)),y=rep(1/8,8)), NA_real_))%>%
            select(sample_name,smp_ideal_score)
          Archaea_df=Archaea_df%>%
            full_join(Zymo_scores)
        }}
        ################################################################################
        message ('210 # Generate alpha diversity metrics')
        ################################################################################
        metrics_df <- matrix_df %>%
          as.data.frame() %>%
          xPlode_sample_name()%>%
          rowwise()%>%
          mutate(
            mn = mean(c_across(any_of(feature_names))),
            sd=sd(c_across(any_of(feature_names))),
            smp_n=sum(c_across(any_of(feature_names))),
            smp_shannon=qShannon(c_across(any_of(feature_names))),
            smp_evenness=qEvenness(c_across(any_of(feature_names))),
            smp_richness = qRichness(c_across(any_of(feature_names)))
          )%>%
          select(-any_of(feature_names))%>%
          left_join(Archaea_df)%>% # piggy backing Zymo scores
          left_join(Propi_df)%>%
          left_join(Lachno_df)%>%
          left_join(Muri_df)%>%
         # left_join(meta_meta_SD)%>%
        #  left_join(num_unique_features_df)%>%
       #   mutate(smp_chimeras=p_chimeric_reads)%>%
          # mutate(smp_percent_good=percentage.of.input.non.chimeric)%>%
          group_by(across(any_of(grouping_columns))) %>%
          mutate(across(starts_with("smp_"), ~mean(., na.rm = TRUE), .names = "{sub('smp_', 'mean_', .col)}"))%>%
          rename_with(~gsub("smp_", "", .), starts_with("ratio_smp_"))%>%
          ungroup()%>%
          select(where(~!all(is.na(.))))
        
        write.csv(metrics_df,paste0(output_data,p_title,'_metrics_df.csv'),row.names = FALSE)
        ################################################################################
      } # skipping above if (make_new_data_files!='yes')
      # read data files from save
      
      if (make_new_data_files!='yes'){metrics_df=read.csv(paste0(output_data,p_title,'_metrics_df.csv'),check.names = FALSE)}
      ################################################################################
      # filter metrics of interest to different groups
      
      result_df=data_df=metrics_df%>%
        pivot_longer(cols = starts_with("smp_"), names_to = "metric")%>%
        # filter(
        #   (Type %in% c('Skin') & metric != 'smp_percent_Muri'& metric != 'smp_percent_Archaea') |
        #     (Type %in% c('Soil', 'Wastewater') & metric != 'smp_percent_Propi'& metric != 'smp_percent_Muri') |
        #     (Type %in% c('Feces') & metric != 'smp_percent_Propi'& metric != 'smp_percent_Archaea') 
        # )%>%
        # filter(metric!='smp_percent_Lachno')%>%
        
        # filter(!(sample_type %in% c('Skin') & (metric %in% c('smp_percent_Muri', 'smp_percent_Archaea'))))%>%
        # filter(!(sample_type %in% c('Soil') & (metric %in% c('smp_percent_Propi', 'smp_percent_Muri'))))%>%
        # filter(!(sample_type %in% c('Feces') & (metric %in% c('smp_percent_Propi', 'smp_percent_Archaea'))))%>%
        # filter(metric!='smp_percent_Lachno')%>%
        
        group_by(across(any_of(grouping_columns)),metric)%>%
        mutate(mean_value=mean(value, na.rm = TRUE))%>%
        ungroup()
      
      #################################################################################
      print('267 # T-test if length(x_axis_group)==2')
      #################################################################################
      unique_values <- unique(data_df[[x_axis_group]])
      
      if (length(unique(data_df[[x_axis_group]])) == 2) {
        t_test_results <- data_df %>%
          filter(!is.nan(value)) %>%
          group_by(across(any_of(test_across)), metric) %>%
          summarise(
            t_test_summary = {
              unique_subgroups <- unique(.data[[x_axis_group]])
              
              # Check if exactly two subgroups exist
              if (length(unique_subgroups) != 2) {
                list(data.frame(p.value = NA, statistic = NA))
              } else {
                # Split data into subgroups and check unique values
                subgroup_counts <- map(unique_subgroups, ~ n_distinct(.data$value[.data[[x_axis_group]] == .x]))
                
                if (any(unlist(subgroup_counts) == 1)) {
                  list(data.frame(p.value = NA, statistic = NA))  # Return NA if any subgroup lacks variation
                } else {
                  list(broom::tidy(wilcox.test(value ~ .data[[x_axis_group]])))  # Run Wilcoxon test
                }
              }
            }, 
            .groups = 'drop'
          ) %>%
          unnest(t_test_summary)
        
        # Extract necessary columns
        t_test_df <- t_test_results %>%
          select(any_of(test_across), metric, p.value)
        
        # Dynamically create the join condition
        join_by <- setNames(c(test_across, "metric"), c(test_across, "metric"))
        
        # Join the results back
        result_df <- data_df %>%
          left_join(t_test_df, by = join_by) %>%
          rename(p_value = p.value) %>%
          mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
      }
      #################################################################################
      message('320 # looping  plots by plot group')
      #################################################################################
      ################################################################################
      # Plotting functions and variables
      ################################################################################
      magnify = 1
      first_height=14
      first_width=14
      
      second_height=14
      second_width=10
      
      p_value_size=3*magnify
      geom_text_size=3*magnify
      ind_text_size=3*magnify
      strip_text_size=11*magnify
      x_text_size=9*magnify
      
      theme_common<-theme_global(base_size = 11)+
        theme(
          plot.tag = element_text(size = 36, face = "bold"),
          legend.position = "bottom",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text( face = 'plain'),
        )
      
 

      
      gPlot <- function(p) {
        p=p+
          scale_color_manual(values = palette_color,labels=palette_label)+
          scale_fill_manual(values = palette_color,labels=palette_label)+
          theme_common+
            guides(
              color = guide_legend(order = 1),
              fill = guide_legend(
                order=2,
                override.aes = list(
                  shape = 22,        
                  size = 5,           
                  color = "black",
                  alpha=1) # adding to legend
          
                ))+
          labs(
            color='',
            fill='',
            #   title = taxa_levs,
            x = '',
            y = '',
            #   caption = ''
          )
        #print(p)
      }
      ################################################################################
      message('320 # begin plotting loop')
      ################################################################################
      plot_list <- list()
      
      list_category=unique(result_df[[plot_group]])
      ################################################################################ 
      for (t in list_category) {
        
        qPrint(unique(result_df[[plot_group]]))
        qPrint(t)
        
        df=result_df%>%
          filter(.data[[plot_group]]==t) %>%
         # filter(metric %in% selected_metric) %>% 
          mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
        
        ###############################################################################
        print('339 # strip backgrounds with outlines')
        ################################################################################
        unique_x_labels <- unique(df[[plot_group]])
        unique_x_labels
        
        outline_color_x='black'
        
        s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = palette_color["'),'"]))')
        print(s)
        eval(parse(t=s))
        backgrounds_x
        
        unique_y_labels =c(rep(unique_x_labels,(length(unique(df$metric))/length(unique_x_labels))))
        unique_y_labels
        outline_color_y='black'
        
        s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),'"]))')
        eval(parse(t=s))
        
        y_strip=as.data.frame(unique_y_labels)%>%
          #mutate(text_color=ifelse(unique_y_labels=='Soil','white','black'))%>%
          mutate(text_color='white')%>%
          mutate(face='bold')%>%
          mutate(text_size=strip_text_size)%>%
          ungroup()
        x_strip=as.data.frame(unique_x_labels)%>%
          #mutate(text_color=ifelse(unique_y_labels=='Soil','white','black'))%>%
          mutate(text_color='white')%>%
          mutate(face='bold')%>%
          mutate(text_size=strip_text_size)%>%
          ungroup()
        
        # }
        ################################################################################
        print('362 # Changing strip labels')
        ################################################################################
        ################################################################################
        print('368 # Plotting')
        ################################################################################ 
        
        summary_df <- df %>%
          select(sample_name,any_of(grouping_columns),metric,value) %>%
          filter(metric=='smp_n')%>%
          distinct() %>%
          group_by(across(any_of(grouping_columns))) %>%
          mutate(max_reads=max(value))%>%
          mutate(y_label=max_reads/2)%>%
          mutate(count = n())%>%
          mutate(metric='n')%>%
         # mutate(mean_cycles=mean(cycle))%>%
          select(-c(sample_name,value)) %>%
          distinct() %>%
          filter(metric %in% selected_metric) %>% 
          mutate(metric = factor(metric, levels = rev(metric_order)))%>%
          mutate(!!x_axis_group := factor(.data[[x_axis_group]], levels = x_axis_group_order))%>%
          mutate(!!x_facet_group := factor(.data[[x_facet_group]], levels = x_facet_group_order))%>%
          mutate(!!plot_group := factor(.data[[plot_group]], levels = plot_group_order))
          #filter(metric %in% metric_filter)
        
        
        #if (length(unique(data_df[[x_axis_group]])) != 2) {
        
        plot_df=df%>%
          mutate(metric = gsub("smp_", "", metric))%>%
          filter(metric %in% selected_metric) %>% 
          mutate(metric = factor(metric, levels = rev(metric_order)))%>%
          mutate(mean_value_label=signif(mean_value,3))%>%
          mutate(mean_value_label=if_else(mean_value_label>1000,round(mean_value,0),mean_value_label))%>%
          mutate(!!x_axis_group := factor(.data[[x_axis_group]], levels = x_axis_group_order))%>%
          mutate(!!x_facet_group := factor(.data[[x_facet_group]], levels = x_facet_group_order))%>%
          mutate(!!plot_group := factor(.data[[plot_group]], levels = plot_group_order))%>%
          ungroup()
          #filter(metric %in% metric_filter)
        
        x_label=(1+length(x_axis_group_order))/2
        
        plot_title=paste(t,'\n',p_title)
        

        p= ggplot(plot_df,aes(x = get(x_axis_group), y = value
                              ,color = get(x_axis_group)
                              )) +
          geom_point(aes(fill=get(plot_group)),alpha=0)+
         # geom_boxplot(aes(fill=get(x_axis_group)),alpha=0)+
          geom_boxplot(width=.6,outlier.size = .5,outlier.shape = NA,
                       #,color='gray60'
                       alpha=.2)+
         # geom_jitter(aes(fill = get(filter1_group)), width = 0.4,shape=21,size=4)+
          # geom_point(aes(fill = get(filter1_group)),position = position_dodge(width = 0.8), 
          #            size = 3, shape=21)+
          # geom_segment(data = segments_df,
          #              aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
          #              color = "black", size = 1) +
        
          geom_text(aes(label = paste0('  ',mean_value_label)),y=0,vjust=.5,hjust=0,color='black', show.legend=FALSE,size=geom_text_size,angle=90)+
          
          # using to see how many samples there are
          geom_text(data=summary_df,aes(y=y_label,label = paste('n =',count)), vjust = -4,hjust=.5,size=geom_text_size,fontface='bold',color='black') +
         # geom_text(data=summary_df,aes(x=x_label,y=y_label,label = 'avg cycles'), vjust = 0,hjust=.5,size=geom_text_size,fontface='bold',color='gray30') +
        #  geom_text(data=summary_df,aes(y=y_label,label = round(mean_cycles,1)), vjust = 2,hjust=.5,size=geom_text_size,fontface='bold',color='gray30') +
          
          facet_grid2(metric~ get(x_facet_group),
                      strip = strip_themed(
                        background_x = backgrounds_x,background_y = backgrounds_y,
                        text_y = elem_list_text(face = y_strip$face,
                                                # size = y_strip$text_size,
                                                color = y_strip$text_color),
                        text_x = elem_list_text(face = x_strip$face,
                                                #   size = x_strip$text_size,
                                                color = x_strip$text_color)),
                      labeller = labeller(.cols = palette_label, .rows = palette_label),
                      scale = "free") +
          scale_y_continuous(limits=c(0,NA))+
          scale_x_discrete(labels = function(x) {
            color <- palette_color[x]  
            label <- palette_label[x]
            #paste0("<span style='color:", color, "'>", x, "</span>")
            paste0("<span style='color:", color, "'>", label, "</span>")
          }) +
          theme(
            axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5,face='bold',size=x_text_size)  # Rotate and enable markdown for x-axis labels
          )+
          labs(
            title = plot_title,
            x = '',
            y = '',
          )
        
          
        p
        #  plot_list[[t]] <- gPlot(p)
        
        # }
        ############################################################################
        
        if (length(unique(data_df[[x_axis_group]])) == 2) {
          
          ############################################################################
          # x_labels <- setNames(df$primer_label, df$primer_label)
          
          plot_df2=df%>%
            mutate(metric = gsub("smp_", "", metric))%>%
            filter(metric %in% selected_metric) %>% 
            # mutate(metric = case_when(metric %in% names(new_metric) ~ new_metric[metric],TRUE ~ metric))%>%
            mutate(metric = factor(metric, levels = rev(metric_order)))%>%
            #mutate(p_label=paste0('p-value = ',signif(p_value, 3)))%>%
            mutate(p_label = paste0('p = ', format(p_value, scientific = TRUE, digits = 3)))%>%
            mutate(significance=ifelse(p_value<.05,'sig','ns'))%>%
            mutate(p_label = ifelse(metric=='n', "", p_label))%>%
            mutate(p_label_colored = paste0("<span style='color:",palette_color[significance],"'>", p_label, "</span>"))%>%
            mutate(mean_value_label=signif(mean_value,3))%>%
            mutate(mean_value_label=if_else(mean_value_label>1000,round(mean_value,0),mean_value_label))%>%
            # group_by(across(any_of(test_across)),metric)%>%
            # mutate(y_p_label=mean(value)/2)%>%
            group_by(across(any_of(plot_group)),metric)%>%
            mutate(y_p_label=max(value)/2)%>%
            
            mutate(!!x_axis_group := factor(.data[[x_axis_group]], levels = x_axis_group_order))%>%
            mutate(!!x_facet_group := factor(.data[[x_facet_group]], levels = x_facet_group_order))%>%
            mutate(!!plot_group := factor(.data[[plot_group]], levels = plot_group_order))%>%
            ungroup()%>%
           # filter(metric %in% metric_filter)%>%
            filter(metric!='n')
          p
  
          q=p+
           # geom_richtext(data= plot_df2,aes(label = p_label, x = 1.5, y = y_p_label ,color=significance), hjust = .5,vjust = 1,size=p_value_size, fontface='bold',inherit.aes = FALSE,show.legend = FALSE,fill='gray90',alpah=.1)+
           # geom_richtext(data= plot_df2,aes(label = p_label, x = 1.5, y = y_p_label ,color=significance), hjust = .5,vjust = 1,size=p_value_size, fontface='bold',inherit.aes = FALSE,show.legend = FALSE,fill=NA)#
              geom_richtext(
             data = plot_df2%>%filter(significance=='sig'),
           #   data = plot_df2,
              aes(
                x = 1.5,
                y = y_p_label,
                label = p_label_colored,
                #color = significance
              #  color = ifelse(plot_df2$significance == "sig", "red", "gray50")
              ),
              color='red',
              inherit.aes = FALSE,
              fill = scales::alpha("gray95", 0.1),  # transparent background fill
         #     label.color = NA,                    # no outline/border
              show.legend = FALSE,
              hjust = 0.5,
              vjust = 1,
              size = p_value_size,
              fontface = "bold"
            ) +
            geom_richtext(
              data = plot_df2%>%filter(significance=='ns' | is.na(significance)),
              aes(
                x = 1.5,
                y = y_p_label,
                label = p_label,
                #color = significance
                #  color = ifelse(plot_df2$significance == "sig", "red", "gray50")
              ),
              color='gray50',
              inherit.aes = FALSE,
              fill = scales::alpha("gray90", 0.05),  # transparent background fill
              #     label.color = NA,                    # no outline/border
              show.legend = FALSE,
              hjust = 0.5,
              vjust = 1,
              size = p_value_size,
              fontface = "bold"
            )
      
          q
          
          p=q
        }
        plot_list[[t]] <- gPlot(p)
      }
      
      
      ############################################################################
      print('# save plots')
      ############################################################################
      # type_plots <- do.call(grid.arrange, c(plot_list, nrow = 1))
    #  type_plots <- wrap_plots(plot_list, n = 1)
      
      combined_plot <- wrap_plots(plot_list, nrow = 1) +
        plot_annotation(
          tag_levels = 'A',
         caption = "ns = not significant; * p < 0.05; ** p < 0.01; *** p < 0.001",
          theme = theme_plot
          )
        
      
      # if(plot_set=='second_set'){
      #   ggsave(paste0(output_plot,p_title,'.png'),type_plots,width=second_width,height=second_height)
      # }else{
      ggsave(paste0(output_plot,'Figure_',figure_number,'_',p_title,'.png'),combined_plot,width=first_width+10,height=first_height)
      # type_plots=grid.arrange(plot_list[['Feces']],plot_list[['Wastewater']],nrow=1,widths = c(1, 1)) 
      #ggsave(paste0(output_plot,p_title,'supplemental.png'),type_plots,width=ind_width,height=ind_height)
    }
  }
}
############################################################################
print('# End')
############################################################################

