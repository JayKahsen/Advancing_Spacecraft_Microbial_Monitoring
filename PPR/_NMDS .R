# Load global settings and helper functions
source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)), '/_globalStuff.R'))
################################################################################
# Script Title & Output Folders
################################################################################
script_title <- 'ordination'
for(plot_set in plot_set_order){
  
  folder=plot_set
  
  
  output_data <- paste0('output_data/', script_title, '/')
  output_plot <- paste0('output_plot/', script_title, '/')
  output_data <- paste0('output_data/',folder,'/')
  output_plot <- paste0('output_plot/',folder,'/')
  
  
  if(plot_set=='standard'){
    normalization_order=c('fixed cycles')
    data_class_order=c('doubleton')
    normalization_name=paste(normalization_order,collapse='_')
    data_class_name=paste(data_class_order,collapse='_')
    
  }
  if(plot_set=='normalization'){
    normalization_order=c('fixed cycles','targeted fluorescence')
    data_class_order=c('doubleton')
    normalization_name=paste(normalization_order,collapse='_')
    data_class_name=paste(data_class_order,collapse='_')
    
  }
  if(plot_set=='data_class'){
    normalization_order=c('targeted fluorescence')
    data_class_order=c('singleton','doubleton')
    normalization_name=paste(normalization_order,collapse='_')
    data_class_name=paste(data_class_order,collapse='_')
    
  }
  
  ################################################################################
  # Dataset selection
  ################################################################################
  #library(rlang)
  starting_r <- 3  # Control starting dataset index
  ending_r <- 3    # Control ending dataset index
  testing='yes';testing_r=4 # run single dataset index
  
  number_of_permutations=1000
  #number_of_permutations=2
  
  
  
  ################################################################################
  # Data processing settings
  ################################################################################
  #make_new_data_files <- 'yes'
  use_custom_labels <- 'no' # can customize strip labels
  ################################################################################
  # Variable sets for analysis
  ################################################################################
  
  #custom_name=paste0('_',folder)# start with underscore
  
  # defaults
  # test_across = c(color_group,plot_group)
  # testing_group=shape_group
  
  parameter_sets <- list(
    set1 = list(
      filter1_group = 'data_class',
      shape_group = 'temperature',
      color_group = 'normalization',
      plot_group = 'sample_type',
      #   testing_group=c('temperature'),
      number_size = 1
    ),
    set2 = list(
      filter1_group = 'data_class',
      shape_group = 'normalization',
      color_group = 'temperature',
      plot_group = 'sample_type',
      #  testing_group=c('normalization'),
      number_size = 1
    ),
    set3 = list(
      filter1_group = 'normalization',
      shape_group = 'data_class',
      color_group = 'temperature',
      plot_group = 'sample_type',
      #  testing_group=c('data_class'),
      number_size = 1
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
      
      # Adjusting for full loop
      data_set=plot_set
      data_set_order=plot_set_order
      
      # Assign analysis parameters dynamically
      list2env(parameter_sets[[paste0('set', match(data_set, data_set_order))]], envir = .GlobalEnv)
      groups <- c('filter1_group', 'shape_group', 'color_group','plot_group')
      for (var in groups) {
        group_value <- get(var)  # e.g., "data_class"
        assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
      }
      grouping_columns=c(shape_group, color_group,plot_group)
      test_across = c(color_group,plot_group)
      testing_group=shape_group
      # group_within=c(shape_group,plot_group)
      
      run_suffix=''
      if(exists('Run_Group')){run_suffix=paste0('_',lp)}
      
      p_title <- paste0(taxa_levs, '_', script_title,custom_name,'_',normalization_name,'_',data_class_name,'_',data_set,run_suffix)
      
      qPrint(p_title)
      
      ################################################################################
      message (' 124 # Load and Filter Data')
      ################################################################################
      matrix_df1=read.csv(matrix_names[r,'file_path'],check.names=FALSE,row.names=1)%>%
        t()%>%as.data.frame()%>%
        xPlode_sample_name()%>%
        filter(.data[[filter1_group]] %in% filter1_group_order)%>%
        filter(.data[[plot_group]] %in% plot_group_order)%>%
        filter(.data[[shape_group]] %in% shape_group_order)%>%
        filter(.data[[color_group]] %in% color_group_order)
      
      if(exists('Run_Group')){
        matrix_df1=matrix_df1%>%filter(.data[[Loop_Group]] %in% lp)
      }
      
      
      matrix_df=matrix_df1%>%
        imPlode_sample_name()%>%
        mutate_all(as.numeric)%>%
        select(which(colSums(.) > 0))
      
      feature_names=names(matrix_df)
      ################################################################################
      # Statistical Analysis
      ################################################################################
      
      
      ################################################################################
      # matrix_dfx = matrix_df + meta
      ################################################################################
      
      matrix_dfx=matrix_df%>%
        as.data.frame()%>%
        xPlode_sample_name()%>%
        ungroup()
      
      ################################################################################
      print('147 # looping plot_group, separate ordinations')
      ################################################################################
      
      #initialize
      pcs_df=loadings_df=percent_explained_df=data.frame()
      
      permANOVA=permDisp=plot_list <- list()
      
      plot_groups=unique(matrix_dfx[[plot_group]])
      
      p_vector=numeric()
      # names(p_vector)=plot_groups
      test_group_vector=f_vector=disp_p_vector=disp_f_vector=p_vector
      
      
      for (t in plot_groups){
        qPrint(plot_groups)
        qPrint(t)
        
        plot_group_df=matrix_dfx%>%
          filter(.data[[plot_group]] == t) %>%
          imPlode_sample_name()%>%
          select(where(~ sum(.) > 0)) 
        
        meta_df=plot_group_df%>%
          xPlode_sample_name()
        ################################################################################
        print('129 # distance matrix,mds, points, percentage explained, assemble types')
        ################################################################################
        
        pseudo_count= 0.5*(min(plot_group_df[plot_group_df>0]))
        pseudo_count_df=plot_group_df+pseudo_count
        normalized_df=pseudo_count_df/rowSums(pseudo_count_df)
        data=as.data.frame(clr(normalized_df))
        
        # Perform PCA
        results=pca_result <- prcomp(data, scale. = FALSE)
        summary(results)
        
        pca_result_df=as.data.frame(pca_result$x)%>%
          ungroup
        pca_result_df=as.data.frame(pca_result$rotation)%>%
          ungroup
        
        # Extract the principal components
        pcs <- as.data.frame(pca_result$x)%>%
          xPlode_sample_name()%>%
          rename(x='PC1',y='PC2')%>%
          select(any_of(meta_names),x,y)
        
        # percent explained
        percent_explained <- format(round(results$sdev^2 / sum(results$sdev^2) * 100, 1), nsmall = 1, trim = TRUE) %>%
          as.data.frame() %>%
          slice(1:2) %>%
          t() %>%
          as.data.frame() %>%
          `colnames<-`(c('x_lab', 'y_lab')) %>%
          mutate(!!plot_group := t)
        
        
        percent_explained_df=rbind(percent_explained_df,percent_explained)
        # loadings_df=rbind(loadings_df,loadings)
        pcs_df=rbind(pcs_df,pcs)
        
        
        #############################################################################
        print('212 # permANOVA')
        #############################################################################
        
        # meta_df=Type_df%>%
        #   xPlode_sample_name()
        
        qPrint(plot_groups)
        qPrint(t)
        
        # Extract PCA scores (e.g., first two principal components)
        pca_scores <- pca_result$x[, 1:2]
        
        # Calculate distance matrix (Euclidean distance as example)
        dist_matrix <- dist(pca_scores)
        
        # Perform PERMANOVA 
        for (tg in testing_group){
          test_key=paste(t,tg,sep='_')
          test_group_vector[test_key]=tg
          
          if(length(unique(plot_group_df[[tg]]))>2){}
          
          s=paste0('result <- adonis2(dist_matrix ~ ',tg,', data = meta_df, permutations = number_of_permutations)')
          print(s)
          eval(parse(t=s))
          
          print(result)
          p_value=result$'Pr(>F)'[1]
          f_statistic <- result$F[1]
          qPrint(p_value)
          qPrint(f_statistic)
          
          p_vector[test_key]=signif(p_value,2)
          f_vector[test_key]=signif(f_statistic,4)
          
          qPrint(c(t,tg))
          qPrint(p_vector)
          
          # View results
          print(result)
          permANOVA[[test_key]]=result
          #############################################################################
          print('199 #permDISP')
          #############################################################################
          
          # permdisp_res <- betadisper(dist_matrix, group = meta_df[[shape_group]]*meta_df[[color_group]])
          s=paste0('permdisp_res <- betadisper(dist_matrix, group = meta_df$',tg,')')
          print(s)
          eval(parse(t=s))
          # View PERMDISP results
          print(permdisp_res)
          summary(permdisp_res)
          
          perm_test <- permutest(permdisp_res, permutations = 999)
          
          print(perm_test)
          
          perm_p<- perm_test$tab[1, "Pr(>F)"]
          print(perm_p)
          perm_f <- perm_test$tab[1, "F"]
          print(perm_f)
          
          # Plot the results
          plot(permdisp_res)
          
          disp_p_vector[test_key]=signif(perm_p,2)
          disp_f_vector[test_key]=signif(perm_f,4)
          permDisp[[test_key]]=perm_f 
          
        }}
      ################################################################################
      # common Theme and gPlot
      ################################################################################
      annotate_text_size=6
      vjust_margin=1.5
      text_size=18
      title_size=text_size+2
      strip_text_size=title_size
      caption_size=text_size
      
      margin_size=20
      loading_text_size=2.1
      
      common_theme=  theme(
        # plot.title = element_blank(),
        plot.background = element_rect(fill = color_palette['general_background']),
        legend.background = element_rect(fill = color_palette['general_background']),
        #legend.title = element_blank(),
        legend.text = element_text(size = title_size),
        axis.title.x = element_markdown(size = title_size,margin = margin(t = 10)),
        axis.title.y = element_markdown(size = title_size),
        axis.text = element_text(size = text_size),
        plot.caption = element_text(hjust = 0,size = caption_size))
      
      
      
      gPlot <- function(p) {
        p=p+
          scale_color_manual(values = color_palette,labels=label_palette)+
          scale_size_manual(values = size_palette,labels=label_palette)+
          scale_shape_manual(values = shape_palette,labels=label_palette)+
          scale_fill_manual(values = color_palette,labels=label_palette)+
          common_theme 
        p  
      }
      
      
      #############################################################################
      print('231 # plotting begins')
      ##############################################################################
      disp_p_label=disp_f_label=p_label=f_label=list()
      
      for (t in plot_groups){print(t)
        permanova_text='\n'
        for (tg in testing_group){
          test_key=paste(t,tg,sep='_')
          
          if (test_key %in% names(p_vector)){
            
            p_label[test_key]=paste0(tg,'PERMANOVA  p = ' ,p_vector[test_key])
            disp_p_label[test_key]=paste0(tg,'PERMDISP  p = ' ,disp_p_vector[test_key])
            
            f_label[test_key]=paste0(tg,'PERMANOVA F = ' ,f_vector[test_key])
            disp_f_label[test_key]=paste0(tg,'PERMDISP  F = ' ,disp_f_vector[test_key])
            
            perm_combined_label=paste0('\t',tg,' PERMANOVA p= ',p_vector[test_key],', F= ',f_vector[test_key],
                                       '; PERMDISP p= ',disp_p_vector[test_key],'; F= ',disp_f_vector[test_key],'\t')
            
            
            permanova_text=paste(permanova_text,perm_combined_label,'\n')
          }}
        
        x_lab=percent_explained_df%>%
          filter(.data[[plot_group]] == t) %>%
          pull(x_lab)
        
        y_lab=percent_explained_df%>%
          filter(.data[[plot_group]] == t) %>%
          pull(y_lab)
        
        df=pcs_df%>%
          filter(.data[[plot_group]] == t)
        
        #df$Primer=factor(df$Primer,levels=primer_order)
        
        # primary_feature_labels <- feature_labels(plot_group) %>%
        #   filter(.data[[plot_group]] == t) %>%
        #   select(-all_of(plot_group))  
        
        
        #lim_adj=2
        
        ###############################################################################
        print(' 342 # strip backgrounds with outlines')
        ################################################################################
        # unique_x_labels <- unique(df[[plot_group]])
        # 
        # outline_color_x='black'
        # 
        # s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = color_palette["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = color_palette["'),'"]))')
        # print(s)
        # eval(parse(t=s))
        # backgrounds_x
        # 
        # unique_y_labels =rev(unique(df[[plot_group]]))
        # unique_y_labels
        # 
        # outline_color_y='black'
        # 
        # s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
        # eval(parse(t=s))
        # 
        # #colors <- color_palette[match(intersect(primer_order,unique(df$Primer)), names(color_palette))]
        # 
        # # unique(df$primer_label)
        # 
        #############################################################################
        ###############################################################################
        message('207 # strip backgrounds with outlines')
        ################################################################################
        unique_y_labels <- rev(unique(df[[plot_group]]))
        
        
        outline_color_y='black'
        
        s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
        print(s)
        eval(parse(t=s))
        text_color=c("white", rep("black", length(unique_y_labels) - 1))
        
        
        if (use_custom_labels=='yes'){
          custom_y_labels=c("slope","1500TF","x3FC","slope","1500TF","x3FC",'24cycles',rep('Soil',3),rep('Zymo',3),'NTC')
          
          custom_y_labels
          unique_y_labels=custom_y_labels
          s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
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
        unique_x_labels =unique(df[[plot_group]])
        
        
        outline_color_x='black'
        
        s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = color_palette["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = color_palette["'),'"]))')
        print(s)
        eval(parse(t=s))
        text_color=c("white", rep("black", length(unique_x_labels) - 1))
        
        
        if (use_custom_labels=='yes'){
          custom_x_labels=c("slope","1500TF","x3FC","slope","1500TF","x3FC",'24cycles',rep('Soil',3),rep('Zymo',3),'NTC')
          
          custom_y_labels
          unique_y_labels=custom_y_labels
          s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
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
        #if(t=='Skin'){stop()}
        
        # df=df%>%
        #   mutate(Type = case_when(Type %in% names(new_Type) ~ new_Type[Type],TRUE ~ Type))%>%
        #   ungroup()
        
        
        permANOVA
        permDisp[[t]]
        summary_text=paste(capture.output(print(permANOVA[[t]]))[-(1:5)],collapse='\n')
        
        capture.output(print(permANOVA[[t]]))[-(1:5)]
        summary_text
        spacer <- rep(" ", 100) |> paste0(collapse = "")
        permDisp[[t]]$shape_test$tab$`Pr(>F)`[1]
        
        p=ggplot(df, aes(x = x, y = y)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_point(aes(fill = !!sym(color_group),shape = !!sym(shape_group),size = !!sym(shape_group)))+  
          xlab(paste('PC1:\t', x_lab, '% Explained')) +
          ylab(paste('PC2:\t', y_lab, '% Explained')) +
          #  stat_ellipse(aes(group = !!sym(shape_group), color = !!sym(shape_group)), linewidth = 1, level = 0.95)+  
          #  stat_ellipse(aes(group = !!sym(color_group), color = !!sym(color_group)), linewidth = .7, level = 0.95)+  
          facet_wrap2(formula(paste("~", plot_group)),strip.position = "right",
                      strip = strip_themed(
                        
                        background_x = backgrounds_x, 
                        text_x = elem_list_text(
                          face = x_strip$face,
                          size = x_strip$text_size,
                          color = x_strip$text_color), 
                        
                        background_y = backgrounds_y,
                        text_y = elem_list_text(
                          face = y_strip$face,
                          size = y_strip$text_size,
                          color = y_strip$text_color)),
                      
                      scale = 'free') +
          
          
          # theme(legend.position = "bottom",
          #       plot.margin = margin(margin_size, margin_size, margin_size / 2, margin_size)) +
          # guides(color = guide_legend(override.aes = list(shape = 22, linewidth = 8))) +
          
          #  annotate("text", x = -Inf, y = Inf, label = paste0(spacer,p_label), hjust = 0, vjust = vjust_margin * 1, color = 'orangered', size = annotate_text_size) +
          # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,disp_p_label), hjust = 0, vjust = vjust_margin * 2, color = 'orangered', size = annotate_text_size) +
          
          # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,f_label),hjust = 0, vjust = vjust_margin * 3, color = 'blue', size = annotate_text_size) +
          # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,disp_f_label), hjust = 0, vjust = vjust_margin * 4, color = 'blue', size = annotate_text_size) +
          annotate("text", x = 0, y = 0, label = summary_text, hjust = 0, vjust = vjust_margin * 4, color = 'blue') +
          
          labs(caption=paste('\n',permanova_text,'\n',summary_text))+
          
          
          guides(
            size = "none",
            color = guide_legend(title='Ellipses',
                                 override.aes = list(linewidth=1.5)), 
            fill = guide_legend(title='Point Color',
                                override.aes = list(shape=22,size=5)), 
            shape = guide_legend(title='Point Shape',
                                 override.aes = list(size=5)), 
            
          )+
          
          ggtitle(paste(p_title, t))
        
        p
        
        if (shape_group %in% colnames(df) && length(unique(df[[shape_group]])) > 1) {
          p <- p + stat_ellipse(aes(group = !!sym(shape_group), color = !!sym(shape_group)), linewidth = 1, level = 0.95)
        }
        
        if (color_group %in% colnames(df) && length(unique(df[[color_group]])) > 1) {
          p <- p + stat_ellipse(aes(group = !!sym(color_group), color = !!sym(color_group)), linewidth = 1, level = 0.95)
        }
        
        p
        
        
        gPlot(p)
        # df %>% filter(primer_label == 'Tr') %>% select(x, y) %>% summary()
        
        #############################################################################
        # df_loadings$ASV_label
        ################################################################################
        #############################################################################
        #library(grid)
        
        # removing legends
        # if( t%in% c('Soil','Wastewater')){p=p+theme(legend.position = "none",
        #                                             plot.margin = margin(margin_size/2, margin_size, margin_size, margin_size))
        # }
        
        
        plot_list[[t]] <- gPlot(p)
        
      }
      
      
      # ############################################################################
      # print('# save plots')
      # ############################################################################
      
      # grid_plots <- do.call(grid.arrange, c(plot_list))
      # 
      # plot_width=12*max(1,floor(length(plot_list)/2))
      # plot_height=12*length(plot_list)/max(1,floor(length(plot_list)/2))
      # ggsave(paste0(output_plot,p_title,'grid.png'),grid_plots,width=plot_width,height=plot_height)
      
      library(patchwork)
      
      plot_width=12*length(plot_list)/max(1,floor(length(plot_list)/2))
      plot_height=12*max(1,floor(length(plot_list)/2))
      
      patch_plots <- wrap_plots(plot_list)
      
      #type_plots=grid.arrange(plot_list[['Skin']],plot_list[['Soil']],ncol=1,heights = c(1.05, 1))
      
      
      ggsave(paste0(output_plot,p_title,'.png'),patch_plots,width=plot_width,height=plot_height)
      
      #ggsave(paste0(output_plot,p_title,'.png'),type_plots,width=plot_width,height=plot_height)
      
    }
  }
}
################################################################################


################################################################################

