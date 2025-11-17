# Load global settings and helper functions
source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)), '/_globalStuff.R'))
################################################################################
# Script Title & Output Folders
################################################################################
script_title <- 'PCA'

if(!exists('plot_set_order')){plot_set_order='none'}
for(plot_set in plot_set_order){qPrint(plot_set)
  
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
  starting_r <- 4  # Control starting dataset index
  ending_r <- 4    # Control ending dataset index
  testing='yes';testing_r=2 # run single dataset index
  
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
      filter1_group = 'sample_time',
      shape_group = 'treatment',
      color_group = 'location',
      plot_group = 'experiment',
      #   testing_group=c('temperature'),
      number_size = 1
      # ),
      # set2 = list(
      #   filter1_group = 'data_class',
      #   shape_group = 'normalization',
      #   color_group = 'temperature',
      #   plot_group = 'sample_type',
      #   #  testing_group=c('normalization'),
      #   number_size = 1
      # ),
      # set3 = list(
      #   filter1_group = 'normalization',
      #   shape_group = 'data_class',
      #   color_group = 'temperature',
      #   plot_group = 'sample_type',
      #   #  testing_group=c('data_class'),
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
      
      if(plot_set!='none'){
        # # Adjusting for full loop
        data_set=plot_set
        data_set_order=plot_set_order
      }
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
      
      p_title <- paste0(taxa_levs, '_', script_title,custom_name,'_',data_set,run_suffix)
      
      qPrint(p_title)
      
      ################################################################################
      message (' 124 # Load and Filter Data')
      ################################################################################
      matrix_df1=read.csv(matrix_names[r,'raw_path'],check.names=FALSE,row.names=1)%>%
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
      print('185 # pre loop initializing')
      ################################################################################
      plot_groups=unique(matrix_dfx[[plot_group]])
      
      pcs_df=loadings_df=percent_explained_df=data.frame()
      
      df_permdisp=df_permanova=NULL
      
      p_vector=numeric()
      # names(p_vector)=plot_groups
      test_group_vector=f_vector=disp_p_vector=disp_f_vector=p_vector
      ################################################################################
      print('199 # Begin plotting group loop for data')
      ################################################################################
      
      for (pg in plot_groups){
        qPrint(plot_groups)
        qPrint(pg)
        
        plot_group_df=matrix_dfx%>%
          filter(.data[[plot_group]] == pg) %>%
          imPlode_sample_name()%>%
          select(where(~ sum(.) > 0)) 
        
        meta_df=plot_group_df%>%
          xPlode_sample_name()
        ################################################################################
        print('212 # distance matrix,mds, points, percentage explained, assemble types')
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
          mutate(!!plot_group := pg)
        
        percent_explained_df=rbind(percent_explained_df,percent_explained)
        # loadings_df=rbind(loadings_df,loadings)
        pcs_df=rbind(pcs_df,pcs)
        
        
        #############################################################################
        print('250 # permANOVA')
        #############################################################################
        
        # meta_df=Type_df%>%
        #   xPlode_sample_name()
        
        qPrint(plot_groups)
        qPrint(pg)
        
        # Extract PCA scores (e.g., first two principal components)
        pca_scores <- pca_result$x[, 1:2]
        
        # Calculate distance matrix (Euclidean distance as example)
        dist_matrix <- dist(pca_scores)
        
        #############################################################################
        print('266 # permANOVA as NMDS BEGIN')
        #############################################################################
        #   by=NULL
        df_meta=meta_df
        run_adonis <- function(formula, by = NULL) {
          set.seed(123)
          formula_string=paste('dist_matrix~',formula)
          formula_str <- as.formula(formula_string)
          print(formula_str)
          
          # Run adonis2
          res <- adonis2(formula_str, data = df_meta, permutations = 1000,by=by)
          res
          # Tidy results
          if(is.null(by)){by=NA}
          
          as.data.frame(res) %>%
            rownames_to_column(var = "Term") %>%
            rename(
              Df = Df,
              SumOfSqs = SumOfSqs,
              R2 = R2,
              F_value = F,
              p_value = `Pr(>F)`
            ) %>%
            mutate(formula=formula)%>%
            mutate(by=by)%>%
            mutate(test='permaNova')
        }
        
        # Function to run PERMDISP and return tidy result
        run_permdisp <- function(dist_matrix, grouping_var) {
          set.seed(123)
          # grouping_var: string, e.g. "color_group"
          
          # Compute betadisper object
          bd <- betadisper(dist_matrix, group = df_meta[[grouping_var]])
          
          # Run permutation test
          perm_disp <- permutest(bd, permutations = 1000)
          
          # Extract table from result
          df <- as.data.frame(perm_disp$tab) %>%
            rownames_to_column(var = "Term") %>%
            mutate(test='permDisp') %>% 
            mutate(grouping_var = grouping_var)
          
          return(df)
        }
        
        
        # Helper function to check if a grouping variable is valid
        is_valid_group <- function(varname, df_meta) {
          var <- df_meta[[varname]]
          var <- var[!is.na(var)]
          length(unique(var)) > 1
        }
        
        # Then wrap your tests
        panov_color <- if (is_valid_group(color_group, df_meta)) run_adonis(formula=color_group) else NULL
        panov_shape <- if (is_valid_group(shape_group, df_meta)) run_adonis(formula=shape_group) else NULL
        panov_filter1 <- if (is_valid_group(filter1_group, df_meta)) run_adonis(formula=filter1_group) else NULL
        
        panov_full_plus_margin <- if (
          all(sapply(c(color_group, shape_group, filter1_group), is_valid_group, df_meta))
        ) run_adonis(formula=paste(color_group,'+',shape_group,'+',filter1_group), by='margin') else NULL
        
        panov_int1_terms <- if (
          all(sapply(c(color_group, shape_group), is_valid_group, df_meta))
        ) run_adonis(formula=paste(color_group,'*',shape_group), by='terms') else NULL
        
        panov_int2_terms <- if (
          all(sapply(c(shape_group, filter1_group), is_valid_group, df_meta))
        ) run_adonis(formula=paste(shape_group,'*',filter1_group), by='terms') else NULL
        
        panov_int3_terms <- if (
          all(sapply(c(color_group, filter1_group), is_valid_group, df_meta))
        ) run_adonis(formula=paste(color_group,'*',filter1_group), by='terms') else NULL
        
        panov_full_onedf <- if (
          all(sapply(c(color_group, shape_group, filter1_group), is_valid_group, df_meta))
        ) run_adonis(formula=paste(color_group,'+',shape_group,'+',filter1_group), by='onedf') else NULL
        
        # Combine all valid results
        df_permanova1 <- list(
          panov_color, panov_shape, panov_filter1,
          panov_full_plus_margin, panov_int1_terms,
          panov_int2_terms, panov_int3_terms,
          panov_full_onedf
        ) %>%
          compact() %>%  # from purrr; removes NULLs
          bind_rows() %>%
          mutate(!!plot_group := pg)
        
        # Similarly for PERMDISP
        disp_color <- if (is_valid_group(color_group, df_meta)) run_permdisp(dist_matrix, color_group) else NULL
        disp_shape <- if (is_valid_group(shape_group, df_meta)) run_permdisp(dist_matrix, shape_group) else NULL
        disp_filter1 <- if (is_valid_group(filter1_group, df_meta)) run_permdisp(dist_matrix, filter1_group) else NULL
        
        df_permdisp1 <- list(disp_color, disp_shape, disp_filter1) %>%
          compact() %>%
          bind_rows() %>%
          mutate(!!plot_group := pg)
        
        print(df_permdisp)
        print(df_permanova)
        
        df_permanova=bind_rows(df_permanova,df_permanova1)
        df_permdisp=bind_rows(df_permdisp,df_permdisp1)
        #############################################################################
        print('374 # permANOVA as NMDS END')
        #############################################################################
        
      }
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
        plot.background = element_rect(fill = palette_color['general_background']),
        legend.background = element_rect(fill = palette_color['general_background']),
        #legend.title = element_blank(),
        legend.text = element_text(size = title_size),
        axis.title.x = element_markdown(size = title_size,margin = margin(t = 10)),
        axis.title.y = element_markdown(size = title_size),
        axis.text = element_text(size = text_size),
        plot.caption = element_text(hjust = 0,size = caption_size))
      
      
      
      gPlot <- function(p) {
        p=p+
          scale_color_manual(values = palette_color,labels=palette_label)+
          scale_size_manual(values = palette_size,labels=palette_label)+
          scale_shape_manual(values = palette_shape,labels=palette_label)+
          scale_fill_manual(values = palette_color,labels=palette_label)+
          common_theme 
        p  
      }
      
      
      ################################################################################
      print('416 # pre loop initializing')
      ################################################################################
      plot_list =plot_list_test<- list()
      
      ################################################################################
      print('421 # Begin plotting group loop for plotting')
      ################################################################################
      disp_p_label=disp_f_label=p_label=f_label=list()
      
      for (pg in plot_groups){print(pg)
        
        ################################################################################
        print('428 # Annotations as NMDS begin')
        ################################################################################
        col_order=c('test','testing','formula','Term','Df','R2','F_value','p_value')
        
        df_permdisp1=df_permdisp%>%
          rename('SumOfSqs'='Sum Sq','F_value'='F','p_value'='Pr(>F)','formula'='grouping_var')
        
        df_tests <- bind_rows(df_permdisp1,df_permanova)%>% 
          filter(.data[[plot_group]] == pg)%>%
          mutate(rownames = row_number())%>%
          mutate(testing=case_when(
            is.na(by)~'Main Effect',
            by=='terms'~'Interaction',
            by=='margin'~'Margin',
            by=='onedf'~'onedf'
          )) %>% 
          filter(!Term %in% c('Residual','Total','Residuals'))%>% 
          mutate(Term = ifelse(
            Term == "locationCrew Sleeping Quarters (subbed for Laundry area)",
            "locationCrew Sleeping Quarters",
            Term
          )) %>% 
          mutate(F_value=round(F_value,2))%>%
          mutate(R2=round(R2,2))%>%
          mutate(p_value=round(p_value,4))%>%
          mutate(short_name=paste(test,testing,formula,sep='_')) %>% 
          mutate(test_name=paste(short_name,rownames,sep='_')) %>% 
          mutate(across(everything(), as.character)) %>% 
          select(any_of(col_order),everything()) %>% 
          ungroup()
        head(df_tests)
        unique(df_tests$test)
        
        tests_text=df_tests %>% 
          filter(is.na(by) | by != 'onedf') %>% 
          select(any_of(col_order))%>%
          mutate(
            label = ifelse(
              test == 'permDisp',
              paste0(test, ' ',testing,' (',formula,'):   Df=',Df,' F=',F_value,' p=',p_value),
              paste0(test, ' ',testing,' (',formula,') ',Term,':   R2=',R2,' Df=',Df,' F=',F_value,' p=',p_value)
            )
          ) %>% 
          pull(label) %>%
          paste(collapse = "\n")
        
        ################################################################################
        print('# perm test plot')
        ################################################################################
        
        df_permdisp1=df_permdisp%>%
          rename('SumOfSqs'='Sum Sq','F_value'='F','p_value'='Pr(>F)','formula'='grouping_var')
        
        col_order=c('test','testing','formula','Term','Df','R2','F_value','p_value')
        
        # Example data
        df_tests_plot <- df_tests %>% 
          pivot_longer(any_of(col_order), names_to = "col", values_to = "val")%>%
          mutate(widths = case_when(
            col == "formula" ~ 5,
            col == "Term" ~ 5,
            col == "testing" ~ 2,
            col == "test" ~ 2,
            TRUE~1
          ))%>%
          
          mutate(col = factor(col, levels = col_order))%>%
          mutate(test_name = factor(test_name, levels = unique(test_name)))
        
        plot_testing=ggplot(df_tests_plot, aes(x = col, y = rev(test_name))) +
          geom_tile(aes(width = widths,fill=short_name), 
                    #fill = "white",
                    color = "gray") +
          geom_text(aes(label = val)) +
          facet_grid(. ~ col, scales = "free_x", space = "free_x")+
          theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x  = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = 'none')
        plot_testing
        
        plot_list_test[[pg]] <- plot_testing
        
        ################################################################################
        print('387 # Annotations as NMDS end')
        ################################################################################
        
        
        x_lab=percent_explained_df%>%
          filter(.data[[plot_group]] == pg) %>%
          pull(x_lab)
        
        y_lab=percent_explained_df%>%
          filter(.data[[plot_group]] == pg) %>%
          pull(y_lab)
        
        df_plot=pcs_df%>%
          mutate(!!sym(filter1_group):=factor(!!sym(filter1_group), levels = filter1_group_order))%>%
          mutate(!!sym(shape_group):=factor(!!sym(shape_group), levels = shape_group_order))%>%
          mutate(!!sym(color_group):=factor(!!sym(color_group), levels = color_group_order))%>%
          filter(.data[[plot_group]] == pg)
        
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
        # s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = palette_color["'),'"]))')
        # print(s)
        # eval(parse(t=s))
        # backgrounds_x
        # 
        # unique_y_labels =rev(unique(df[[plot_group]]))
        # unique_y_labels
        # 
        # outline_color_y='black'
        # 
        # s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),'"]))')
        # eval(parse(t=s))
        # 
        # #colors <- palette_color[match(intersect(primer_order,unique(df$Primer)), names(palette_color))]
        # 
        # # unique(df$primer_label)
        # 
        #############################################################################
        ###############################################################################
        message('207 # strip backgrounds with outlines')
        ################################################################################
        unique_y_labels <- rev(unique(df_plot[[plot_group]]))
        
        
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
        unique_x_labels =unique(df_plot[[plot_group]])
        
        
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
        ################################################################################
        print('# plot sample_time')
        plot_title=paste(p_title,pg,'sample_time',sep='_')
        plot_title
        ################################################################################
        p=ggplot(df_plot, aes(x = x, y = y)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_point(aes(fill = !!sym(filter1_group),
                         #shape = !!sym(shape_group),size = !!sym(shape_group)
          ),
          shape=21,size=5
          )+  
          geom_path(aes(group=filter1_group),size = .8,alpha=.5) +
          facet_grid2(formula(paste(shape_group,"~",color_group)),
                      labeller = labeller(
                        .rows = function(value) palette_label[as.character(value)],
                        .cols = function(value) palette_label[as.character(value)]
                      ),
                      #strip.position = "right",
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
                      #  scale = 'free'
          ) +
          guides(
            size = "none",
            color = guide_legend(title='Ellipses',override.aes = list(linewidth=1.5)), 
            fill = guide_legend(title='Point Color',override.aes = list(shape=22,size=5)), 
            shape = guide_legend(title='Point Shape',override.aes = list(size=5)), 
          )+
          labs(
            x=paste('PC1:\t', x_lab, '% Explained'),
            y=paste('PC2:\t', y_lab, '% Explained'),
            title=plot_title,
            #caption=tests_text
            )+
          theme(
            # Add border around the panel
            panel.border = element_rect(color = "gray30", fill = NA, size = 1),
            
            # Add border around the strip background
            strip.background = element_rect(color = "gray30", size = 1),
            legend.position = 'bottom'
            
          )
        p
        plot_sample_time=gPlot(p)
        ggsave(paste0(output_plot,plot_title,'.png'),plot_sample_time,width=36,height=10)
        ################################################################################ 
        print('# plot normal plot')
        plot_title=paste(p_title,pg,sep='_')
        plot_title
        ###############################################################################
        
        p=ggplot(df_plot, aes(x = x, y = y)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_point(aes(fill = !!sym(color_group),shape = !!sym(shape_group),size = !!sym(shape_group)))+  

          #  stat_ellipse(aes(group = !!sym(shape_group), color = !!sym(shape_group)), linewidth = 1, level = 0.95)+  
          #  stat_ellipse(aes(group = !!sym(color_group), color = !!sym(color_group)), linewidth = .7, level = 0.95)+  
          # facet_grid2(formula(paste(color_group,"~", filter1_group)),
          #strip.position = "right",
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
                      
                      scale = 'free'
          ) +
          
          
          # theme(legend.position = "bottom",
          #       plot.margin = margin(margin_size, margin_size, margin_size / 2, margin_size)) +
          # guides(color = guide_legend(override.aes = list(shape = 22, linewidth = 8))) +
          
          #  annotate("text", x = -Inf, y = Inf, label = paste0(spacer,p_label), hjust = 0, vjust = vjust_margin * 1, color = 'orangered', size = annotate_text_size) +
          # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,disp_p_label), hjust = 0, vjust = vjust_margin * 2, color = 'orangered', size = annotate_text_size) +
          
          # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,f_label),hjust = 0, vjust = vjust_margin * 3, color = 'blue', size = annotate_text_size) +
          # annotate("text", x = -Inf, y = Inf, label = paste0(spacer,disp_f_label), hjust = 0, vjust = vjust_margin * 4, color = 'blue', size = annotate_text_size) +
          #annotate("text", x = 0, y = 0, label = summary_text, hjust = 0, vjust = vjust_margin * 4, color = 'blue') +
          
          labs(
            x=paste('PC1:\t', x_lab, '% Explained'),
            y=paste('PC2:\t', y_lab, '% Explained'),
            title=plot_title,
            caption=tests_text)+
          
          
          guides(
            size = "none",
            color = guide_legend(title='Ellipses',
                                 override.aes = list(linewidth=1.5)), 
            fill = guide_legend(title='Point Color',
                                override.aes = list(shape=22,size=5)), 
            shape = guide_legend(title='Point Shape',
                                 override.aes = list(size=5)), 
          )
        
        p
        
        if (shape_group %in% colnames(df_plot) && length(unique(df_plot[[shape_group]])) > 1) {
          p <- p + stat_ellipse(aes(group = !!sym(shape_group), color = !!sym(shape_group)), linewidth = 2, level = 0.95)
        }
        
        if (color_group %in% colnames(df_plot) && length(unique(df_plot[[color_group]])) > 1) {
          p <- p + stat_ellipse(aes(group = !!sym(color_group), color = !!sym(color_group)), linewidth = 1, level = 0.95)
        }
        if (filter1_group %in% colnames(df_plot) && length(unique(df_plot[[filter1_group]])) > 1) {
          p <- p + stat_ellipse(aes(group = !!sym(filter1_group), color = !!sym(filter1_group)), linewidth = 1, level = 0.95)
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
        
        
        plot_list[[pg]] <- gPlot(p)
        
      }
      
      
      
      
      
      # ############################################################################
      # print('# save plots')
      # ############################################################################
      
      plot_width=12*length(plot_list)/max(1,floor(length(plot_list)/2))
      plot_height=12*max(1,floor(length(plot_list)/2))
      
      patch_plots <- wrap_plots(plot_list)
      print(patch_plots)
      ggsave(paste0(output_plot,p_title,'.png'),patch_plots,width=plot_width+4,height=plot_height+4)
      
      patch_test_plots <- wrap_plots(plot_list_test)
      print(patch_test_plots)
      ggsave(paste0(output_plot,p_title,'permTests.png'),patch_test_plots,width=plot_width,height=plot_height/2)
      
    } # end plotset
  } # end loop/run group
  
} # end asv levels

################################################################################
################################################################################

################################################################################
################################################################################




################################################################################

