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
testing='yes';testing_r=1 # run single dataset index
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


figure_order=c('2','2A','2B')
parameter_sets <- list(
  set1 = list(
    filter1_group = 'figure',
    y_axis_group = 'bar_name',
    x_facet_group = 'bar_color',
    y_facet_group = 'experiment',
    number_size = 1
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
    df_matrix1=read.csv(matrix_names[r,'raw_path'],check.names=FALSE,row.names=1)%>%
      t()%>%as.data.frame()%>%
      xPlode_sample_name()%>%
      filter(.data[[filter1_group]] %in% filter1_group_order)%>%
      filter(.data[[y_facet_group]] %in% y_facet_group_order)%>%
      filter(.data[[x_facet_group]] %in% x_facet_group_order)%>%
      filter(.data[[y_axis_group]] %in% y_axis_group_order)
    

    if(exists('Run_Group')){
      df_matrix1=df_matrix1%>%filter(.data[[Loop_Group]] %in% lp)
    }
    
    qPrint(p_title)
    
    df_matrix=df_matrix1%>%
      imPlode_sample_name()%>%
      mutate_all(as.numeric)%>%
      select(which(colSums(.) > 0))
    
    feature_names=names(df_matrix)
    ################################################################################
    # Statistical Analysis
    ################################################################################
    kw_group_results <- Kruskal_Mann_Whitney_Test(df = df_matrix, testing_group = y_axis_group, category_group = test_across)

    ################################################################################
    message (' 137   # Data Processing and Normalization')
    ################################################################################
    # get feature labels to 
    feature_labels_df=feature_labels(y_facet_group)
    
    y_axis_group_order=rev(y_axis_group_order)

    plot_df <- df_matrix %>%
      xPlode_sample_name()%>%
      pivot_longer(cols = any_of(feature_names), names_to = 'feature', values_to = 'counts') %>%
      # left_join(feature_labels_df %>% select(feature, taxa)) %>%
      # mutate(Domain = sapply(taxa, extract_domain)) %>%
      group_by(feature, across(any_of(y_facet_group))) %>%
      mutate(grouped_feature_counts = sum(counts))%>%
      filter(grouped_feature_counts > 0)%>%
      group_by(sample_name) %>%
      mutate(sample_counts = sum(counts)) %>%
      ungroup() %>%
      mutate(rel_abun = counts / sample_counts) %>%
      #group_by(feature, Domain, across(any_of(grouping_columns))) %>%
      group_by(feature, across(any_of(grouping_columns))) %>%
      summarize(group_mean_counts = mean(counts),
                sample_counts = mean(sample_counts),
                mean_rel_abun = mean(rel_abun), .groups = 'drop') %>%
      group_by(across(any_of(y_facet_group))) %>%
      mutate(group_counts = mean(sample_counts)) %>%
      ungroup() %>%
      left_join(feature_labels_df %>% select(-c(taxa)), by = c('feature', y_facet_group)) %>%
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group))))%>%
      mutate(!!sym(y_facet_group):=factor(!!sym(y_facet_group), levels = y_facet_group_order))%>%
      mutate(!!sym(y_axis_group):=factor(!!sym(y_axis_group), levels = y_axis_group_order))%>%
      mutate(!!sym(x_facet_group):=factor(!!sym(x_facet_group), levels = x_facet_group_order))%>%
      ungroup()

    mn_plot_df=plot_df%>% # mean normalized
      group_by(feature,across(any_of(test_across)))%>%
      mutate(group_rel_abun=mean(mean_rel_abun))%>%
      ungroup()%>%
      mutate(normalized_rel_abun=mean_rel_abun/group_rel_abun)%>%
      #select(ASV,Domain,Type,mean_normalized,plot_order,ASV_label2)%>%
      left_join(kw_group_results, by = c("feature",test_across)) %>%
      ungroup()%>%
      group_by(feature,across(any_of(test_across)))%>%
      mutate(log2_value = log2(normalized_rel_abun + 1e-6))%>%
      mutate(geo_mean=exp(mean(log(mean_rel_abun+1e-6))),
             log2_ratio=log2((mean_rel_abun+1e-6)/geo_mean)) %>% 
      ungroup() %>% 
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group))))%>%
      mutate(!!sym(y_facet_group):=factor(!!sym(y_facet_group), levels = y_facet_group_order))%>%
      mutate(!!sym(y_axis_group):=factor(!!sym(y_axis_group), levels = y_axis_group_order))%>%
      mutate(!!sym(x_facet_group):=factor(!!sym(x_facet_group), levels = x_facet_group_order))%>%
      ungroup()
    
    # df_check=mn_plot_df %>% 
    #   select(feature,normalized_rel_abun,log2_value,log2_ratio,geo_mean,everything())%>%
    #   filter(location=='EVA Module')
    
    
    names(mn_plot_df)
    names(kw_group_results)
    ################################################################################
    message( '164 # Generate Plots')
    ################################################################################
    ################################################################################
    # Plotting functions and variables
    ################################################################################
    if(taxa_levs=='Phylum'){magnify=1}
    
    margin_size=10
    annotate_text_size=6*magnify
    axis_text_size=16*magnify
    axis_title_text_size=18*magnify
    strip_text_size=20*magnify
    legend_text_size=20*magnify
    
    common_theme=  theme(
      strip.text = element_text(face='bold',size = strip_text_size),
      strip.text.y = element_text(face='bold',size = strip_text_size),
      plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
      plot.background = element_rect(fill = palette_color['general_background']),
      legend.background = element_rect(fill = palette_color['general_background']),
      # legend.title = element_blank(),
      legend.text = element_text(size = legend_text_size),
      axis.title.y = element_text(face='bold',size = axis_title_text_size,margin = margin(r = 10)),
      axis.title.x = element_text(face='bold',size = axis_title_text_size,margin = margin(t = 10)),
      #plot.title = element_blank(),
      plot.caption = element_text(hjust = 0.5))
    
    
    gPlot <- function(p) {
      p=p+
        common_theme +
        scale_fill_manual(values = palette_color,guide = guide_legend(reverse = TRUE),labels=palette_label) +
        scale_color_manual(values = palette_color,guide = guide_legend(reverse = TRUE),labels=palette_label)+
        labs(
          title = p_title,
          #   x = "",
             y = 'Phyla arranged by average realtive abundance',
          #    caption = ''  
        )
      print(p)
      p  
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
    #big_constant_number=1000 # needed?
    sc=1000 # shifts bars on left plot
    sc=1
    x_min=2/sc;x_alt=10 # modifies mean count text
    magnify=3# controls text size
    ################################################################################

    
    st_plot=ggplot(plot_df, aes(x = mean_rel_abun*group_counts*sc, y = plot_order)) +
      geom_col(aes(fill = !!sym(y_axis_group)),width=.6,position = position_dodge(width = 0.8))+
      scale_y_discrete(breaks = custom_breaks, labels = custom_labels)+
      facet_grid2(formula(paste(y_facet_group,'~',x_facet_group)),
                  strip = strip_themed(
                    background_y = backgrounds_y,
                    text_y = elem_list_text(
                      face = y_strip$face,
                      size = y_strip$text_size,
                      color = y_strip$text_color),
                    background_x = backgrounds_x,
                    text_x = elem_list_text(
                      face = x_strip$face,
                      size = x_strip$text_size,
                      color = x_strip$text_color) 
                  ),scale = "free",space = "free") +
      scale_x_log10(labels = sc_log_labels,breaks = log_breaks)
    st_plot
    
    mn_plot=ggplot(mn_plot_df, aes(x = normalized_rel_abun, y = plot_order, fill = !!sym(y_axis_group))) +
    #  mn_plot=ggplot(mn_plot_df, aes(x = 2^log2_ratio, y = plot_order, fill = !!sym(y_axis_group))) +
      geom_col(aes(fill = !!sym(y_axis_group)),width=.8,position = position_dodge(width = 0.8))+
      geom_text(aes(x=Inf,label=significance,color=adjusted_significance),hjust=1.1)+
      geom_vline(aes(xintercept=1),color='magenta',linetype='dashed')+
      facet_grid2(formula(paste(y_facet_group,'~',x_facet_group)),
                  #facet_grid2(formula(paste(primary_group,'+',subgroup2, "~.")),
                  #facet_grid2(formula(paste(primary_group, "~.")),
                  strip = strip_themed(
                    background_y = backgrounds_y,
                    text_y = elem_list_text(
                      face = y_strip$face,
                      size = y_strip$text_size,
                      color = y_strip$text_color),
                    background_x = backgrounds_x,
                    text_x = elem_list_text(
                      face = x_strip$face,
                      size = x_strip$text_size,
                      color = x_strip$text_color) 
                  ),scale = "free",space = "free") +
      scale_x_continuous(trans = "log2",labels = log_labels2) +
      scale_y_discrete(breaks = custom_breaks, labels = custom_labels,position='right')+
      coord_cartesian(xlim = c(1/2^zoom, +2^zoom))+
      theme(legend.position = "none") 
    
    

    
    
    mn_plot
    
    if(taxa_levs=='Phylum'){
      
      magnify=1
      
      y_labels <- setNames(plot_df$feature_label2, plot_df$plot_order)
      
      st_plot=st_plot+
        geom_text(
          aes(
            label = signif(group_mean_counts, 3),
            x = ifelse(group_mean_counts < x_min, x_alt, group_mean_counts * sc),
            color = ifelse(group_mean_counts < x_min, 'black', 'white'),
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
      
      
      y_labels <- setNames(mn_plot_df$feature_label2, mn_plot_df$plot_order)
      
      mn_plot=mn_plot+
        scale_y_discrete(labels = y_labels,position='right')+ 
        theme(
          axis.text.y = element_markdown('bold',size=axis_text_size),
          legend.position = "none") 
      mn_plot
      
    }
    
    
    st_plot=gPlot(st_plot)
    mn_plot=gPlot(mn_plot)
    
    combined_plot <- st_plot + mn_plot + plot_layout(nrow = 1, widths = c(1.2, 1))
    
    match(data_set, data_set_order) # optional way to change plots
    
    if(taxa_levs=='Phylum'){
      ggsave(paste0(output_plot,p_title,'.png'),combined_plot,width=36,height=24)
    }else{
      ggsave(paste0(output_plot,p_title,'.png'),combined_plot,width=28,height=48)
    }
    
   
    
    ################################################################################
    message('277 # bar plots')
    interleaved_palette <- c(
  "#67000D", "#D8B365", "#6BAED6", "#FDBE85", "#00441B", "#6A51A3",
  "#A50F15", "#FFDA66", "#C6DBEF", "#F16913", "#238B45", "#807DBA",
  "#CB181D", "#FFFFB2", "#2171B5", "#7F2704", "#66C2A4", "#3F007D",
  "#EF3B2C", "#8C510A", "#08306B", "#D94801", "#C7E9C0", "#DADAEB")
    
    interleaved_palette=palette_common
    ################################################################################
    # df2<- df_matrix %>%
    #   xPlode_sample_name()%>%
    #   pivot_longer(cols = any_of(feature_names), names_to = 'feature', values_to = 'counts') %>%
    #   group_by(sample_name)%>%
    #   mutate(sample_counts=sum(counts))%>%
    #   ungroup()%>%
    #   mutate(sample_rel_abun=counts/sample_counts)%>%
    #   group_by(across(any_of(y_facet_group))) %>%
    #   mutate(y_facet_counts=sum(counts))%>%
    #   group_by(feature, across(any_of(y_facet_group))) %>%
    #   mutate(y_facet_feature_counts = sum(counts))%>%
    #   mutate(mean_rel_abun=mean(sample_rel_abun))%>%
    #   ungroup()%>%
    #   filter(y_facet_feature_counts > 0)%>%
    #   mutate(sample_rel_abun=counts/sample_counts)%>%
    #   mutate(y_facet_rel_abun=y_facet_feature_counts/y_facet_counts)%>%
    #   ungroup()
     
    plot_list=list()
    ratio_vector=numeric()
    #names(ratio_vector)=y_facet_group_order
    for(facet in y_facet_group_order){
      
      ###############################################################################
      message('207 # strip backgrounds with outlines')
      ################################################################################
      unique_y_labels <- facet
      
      
      outline_color_y='black'
      
      s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),'"]))')
      print(s)
      eval(parse(t=s))
      #text_color=c("white", rep("black", length(unique_y_labels) - 1))

      y_strip=as.data.frame(unique_y_labels)%>%
       # mutate(text_color=ifelse(unique_y_labels=='Soil','white','black'))%>%
        mutate(text_color='white')%>%
        mutate(face='bold')%>%
        mutate(text_size=strip_text_size)%>%
        ungroup()

      ###############################################################################
    
    df_plot2<-plot_df%>%
      select(any_of(y_facet_group),group_mean_rel_abun,feature,plot_order)%>%
      filter(group_mean_rel_abun>.001)%>%
      distinct()%>%
      mutate(!!sym(y_facet_group):=factor(!!sym(y_facet_group), levels = y_facet_group_order))%>%
    #  mutate(!!sym(y_axis_group):=factor(!!sym(y_axis_group), levels = y_axis_group_order))%>%
     # mutate(!!sym(x_facet_group):=factor(!!sym(x_facet_group), levels = x_facet_group_order))%>%
      arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group))))%>%
      filter(.data[[y_facet_group]] %in% facet)%>%
        arrange(desc(group_mean_rel_abun))%>%
        mutate(plot_order=factor(plot_order,levels=unique(plot_order)))%>%
        ungroup()
  
      # 
      n=nrow(df_plot2)
      # default_palette=scales::hue_pal()(n)
      # step=3
      # idx <- (seq(0, n - 1) * step) %% n + 1
      # shuffled_palette <- default_palette[idx]
      

    plot_title=paste0(taxa_levs, '_taxonomy_bar_',data_set,run_suffix)
    
  p=ggplot(data=df_plot2,aes(x=!!sym(y_facet_group),y=group_mean_rel_abun,fill=plot_order))+
    geom_col(stat = "identity", position = position_stack(reverse = TRUE)) +
    facet_grid2(formula(paste(y_facet_group,'~','.')),
                #facet_grid2(formula(paste(primary_group,'+',subgroup2, "~.")),
                #facet_grid2(formula(paste(primary_group, "~.")),
                strip = strip_themed(
                  background_y = backgrounds_y,
                  text_y = elem_list_text(
                    face = y_strip$face,
                    size = y_strip$text_size,
                    color = y_strip$text_color)),
                scale = "free",space = "free") +
   # scale_fill_manual(values = interleaved_palette)+
    theme(legend.text = element_markdown())+
    scale_fill_manual(
      values = interleaved_palette,  # palette in desired stack order (bottom first)
      guide = guide_legend(reverse = TRUE)  # flip legend only
    )+
    labs(title='',
         fill='',
         x='',
         y='')
  p  
    
  plot_list[[facet]]=p
  ratio_vector[facet]=n
  
    }
    ratio_vector
    type_plots <- wrap_plots(rev(plot_list),ncol=1,heights = rev(ratio_vector))
    #type_plots <- wrap_plots(plot_list)
    type_plots
    
    library(patchwork)
    
    library(patchwork)
    
    # # Extract plots
    # left_plot  <- plot_list[['Soil']]
    # top_right  <- plot_list[['Skin']]
    # bottom_right <- plot_list[['Feces']]
    # 
    # # Stack right-side plots vertically with equal height
    # right_plots <- top_right / bottom_right + 
    #   plot_layout(heights = ratio_vector[3:2])  # equal vertical heights
    #   #plot_layout(heights = c(1, 1))  # equal vertical heights
    # 
    # final_layout <- (left_plot | right_plots) + 
    #   plot_layout(widths = c(1, 1)) +  
    #   plot_annotation(title = plot_title)
    # 
    # 
    # final_layout
    
   
    
    
  
    
    
    #ggsave(paste0(output_plot,plot_title,'.png'),type_plots)
    ggsave(paste0(output_plot,plot_title,'.png'),type_plots,width=10,height=8)
      
      
      # strip.text = element_text(face='bold',size = strip_text_size),
      #           strip.text.y = element_text(face='bold',size = strip_text_size),
      #           plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
      #           plot.background = element_rect(fill = palette_color['general_background']),
      #           legend.background = element_rect(fill = palette_color['general_background']),
      #           # legend.title = element_blank(),
      #           legend.text = element_markdown(size = legend_text_size),
      #           axis.title.y = element_text(face='bold',size = axis_title_text_size,margin = margin(r = 10)),
      #           axis.title.x = element_text(face='bold',size = axis_title_text_size,margin = margin(t = 10)),
      #           #plot.title = element_blank(),
      #           plot.caption = element_text(hjust = 0.5))

    ################################################################################
    message('277 # bar plots')
    ################################################################################
      
    
  }
}

#}