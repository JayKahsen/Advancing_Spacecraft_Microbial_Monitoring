# Load global settings and helper functions
source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)), '/_globalStuff.R'))
################################################################################
# Script Title & Output Folders
################################################################################
script_title <- 'NMDS'

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
  starting_r <- 3  # Control starting dataset index
  ending_r <- 3    # Control ending dataset index
  testing='yes';testing_r=5 # run single dataset index
  
  number_of_permutations=1000
  #number_of_permutations=2
  trymax_value=20
  
  k_value=2
  calculate_new_nmds='yes'
  
  
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
  
  ################################################################################
  # optional Run_Group
  
  # Run_Group='primer_set' # will make separate runs for each value in subgroup_order
  # Run_Group_order='Standard'
  ################################################################################
  # functions 
  ################################################################################
  
  # Extract NMDS Scores (only first 2 dimensions)
  extract_nmds_scores <- function(nmds_model, model_name) {
    scores <- as.data.frame(scores(nmds_model, display = "sites")[, 1:2])
    colnames(scores) <- c(paste0("x_", model_name), paste0("y_", model_name))
    scores$sample_name <- rownames(scores)
    return(scores)
  }
  
  # Extract NMDS Species Scores (only first 2 columns)
  extract_nmds_species_scores <- function(nmds_model, model_name) {
    # Extract site and species scores for the first two NMDS dimensions
    scores <- as.data.frame(scores(nmds_model, display = "sites")[, 1:2])
    loadings <- as.data.frame(scores(nmds_model, display = "species")[, 1:2])
    
    # Calculate combined R value for each species as the Euclidean norm of the first two dimensions
    combined_R <- sqrt(loadings[, 1]^2 + loadings[, 2]^2)
    
    # Normalize combined R to range 0-1 for comparison
    normalized_R <- combined_R / max(combined_R)
    loadings$normalized_R <- normalized_R  # Add normalized R as a third column
    
    # Scaling factor calculation for plotting
    min_plot_dim <- min(max(abs(c(min(scores[, 1]), max(scores[, 1]), 
                                  min(scores[, 2]), max(scores[, 2])))))
    max_loading <- max(abs(c(loadings[, 1], loadings[, 2])))
    scf <- min_plot_dim / max_loading
    
    # Apply scaling factor to loadings for plotting
    loadings[, 1] <- loadings[, 1] * scf
    loadings[, 2] <- loadings[, 2] * scf
    
    # Rename columns
    colnames(loadings)[1:3] <- c(paste0("x_", model_name), paste0("y_", model_name), paste0("R_", model_name))
    loadings$Species <- rownames(loadings)
    
    return(loadings)
  }
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
      df_matrix1=read.csv(matrix_names[r,'file_path'],check.names=FALSE,row.names=1)%>%
        t()%>%as.data.frame()%>%
        xPlode_sample_name()%>%
        filter(.data[[filter1_group]] %in% filter1_group_order)%>%
        filter(.data[[plot_group]] %in% plot_group_order)%>%
        filter(.data[[shape_group]] %in% shape_group_order)%>%
        filter(.data[[color_group]] %in% color_group_order)
      
      if(exists('Run_Group')){
        df_matrix1=df_matrix1%>%filter(.data[[Loop_Group]] %in% lp)
      }
      
      df_matrix=df_matrix1%>%
        imPlode_sample_name()%>%
        mutate_all(as.numeric)%>%
        select(which(colSums(.) > 0))
      
      feature_names=names(df_matrix)
      
      ################################################################################
      # dfx_matrix = df_matrix + meta
      ################################################################################
      
      dfx_matrix=df_matrix%>%
        as.data.frame()%>%
        xPlode_sample_name()%>%
        ungroup()
      
      df_richness=dfx_matrix %>% 
        rowwise()%>%
        mutate(smp_richness = qRichness(c_across(any_of(feature_names)))) %>% 
        ungroup() %>% 
        select(sample_name,smp_richness)
      ################################################################################
      print('185 # pre loop initializing')
      ################################################################################
      plot_groups=unique(dfx_matrix[[plot_group]])
      
      pcs_df=loadings_df=percent_explained_df=data.frame()
      
      df_permdisp=df_permanova=NULL
      
      p_vector=numeric()
      # names(p_vector)=plot_groups
      test_group_vector=f_vector=disp_p_vector=disp_f_vector=p_vector
      ################################################################################
      print('199 # Begin plotting group loop for data')
      for (pg in plot_group_order){
        qPrint(plot_group)
        qPrint(pg)
        ################################################################################
        df_plot_group=dfx_matrix%>%
          filter(.data[[plot_group]] == pg) %>%
          imPlode_sample_name()%>%
          select(where(~ sum(.) > 0)) 
        
        df_meta=df_plot_group%>%
          xPlode_sample_name()%>%
          select(any_of(meta_names))%>%
          mutate(temp=sample_name)%>%
          column_to_rownames(var='sample_name')%>%
          mutate(!!sym(filter1_group):=factor(!!sym(filter1_group), levels = filter1_group_order))%>%
          mutate(!!sym(shape_group):=factor(!!sym(shape_group), levels = shape_group_order))%>%
          mutate(!!sym(color_group):=factor(!!sym(color_group), levels = color_group_order))
        
        df_normalized=df_plot_group/rowSums(df_plot_group)
        ################################################################################
        print('212 # distance matrix,mds, points, percentage explained, assemble types')
        ################################################################################
        #  if(exists('calculate_new_nmds')){
        set.seed(123)
        
        nmds_bray <- metaMDS(df_normalized, k = k_value, 
                             maxit = number_of_permutations, 
                             trymax = trymax_value,
                             autotransform = FALSE, # ask
                             # init = init,
                             distance = "bray",
                             wascores = TRUE)
        
        nmds_bray_scores <- extract_nmds_scores(nmds_bray, "nmds_bray")
        
        nmds_bray_species_scores <- extract_nmds_species_scores(nmds_bray, "nmds_bray")
        
        df_nmds=nmds_bray_scores%>%
          mutate(stress= nmds_bray$stress)%>%
          mutate(!!plot_group:=pg)
        
        df_NMDS_loadings=nmds_bray_species_scores%>%
          mutate(!!plot_group:=pg)
        
        #############################################################################
        print('266 # permANOVA as NMDS BEGIN')
        #############################################################################
        dist_matrix <- vegdist(df_normalized, method = "bray")
        
        if(!exists('df_meta')){df_meta=meta_df}
        run_adonis <- function(formula, method = NULL) {
          set.seed(123)
          formula_string <- paste('dist_matrix ~', formula)
          formula_str <- as.formula(formula_string)
          print(formula_str)
          
          # Run adonis2 with renamed argument
          res <- adonis2(formula_str, data = df_meta, permutations = 1000, by = method)
          
          if (is.null(method)) method <- NA
          
          # Tidy the results
          as.data.frame(res) %>%
            rownames_to_column(var = "Term") %>%
            rename(
              Df = Df,
              SumOfSqs = SumOfSqs,
              R2 = R2,
              F_value = F,
              p_value = `Pr(>F)`
            ) %>%
            mutate(
              formula = formula,
              method = method,
              test = "permANOVA"
            ) %>% 
            select(method,formula,Term,everything())
        }
        
        
        
        
        # Function to run PERMDISP and return tidy result
        run_permdisp <- function(dist_matrix, grouping_var) {
          # Example usage: run_permdisp(dist_matrix, "Treatment")
          
          set.seed(123)
          
          # Compute betadisper object using median
          bd <- betadisper(dist_matrix, group = df_meta[[grouping_var]], type = "median")
          
          # Run permutation test
          perm_disp <- permutest(bd, permutations = 1000)
          
          # Extract distances to median
          dist_to_median <- bd$distances
          
          # Compute group-wise average distances
          group_means <- tapply(dist_to_median, bd$group, mean, na.rm = TRUE)
          group_means <- round(group_means, 2)
          
          # Create formatted string
          average_dist_median <- paste0(names(group_means), "=", group_means, collapse = "___")
          
          # Format output table
          result_df <- as.data.frame(perm_disp$tab) %>%
            rownames_to_column(var = "Term") %>%
            mutate(
              test = "permDISP",
              grouping_var = grouping_var,
              average_dist_median = average_dist_median
            )
          
          return(result_df)
        }
        
        
        
        # Helper function to check if a grouping variable is valid
        is_valid_group <- function(varname, df_meta) {
          var <- df_meta[[varname]]
          var <- var[!is.na(var)]
          length(unique(var)) > 1
        }
        
        perm_group <- unique(c(
          if (is_valid_group(color_group, df_meta)) color_group,
          if (is_valid_group(shape_group, df_meta)) shape_group,
          if (is_valid_group(filter1_group, df_meta)) filter1_group
        ))
        
        
        panov_margin<- if (length(perm_group) > 0)
          run_adonis(formula = paste(perm_group, collapse = " + "),method='margin')
        else
          NULL
        
        
        
        
        panov_terms <- NULL
        
        if (is_valid_group(color_group, df_meta) && is_valid_group(shape_group, df_meta)) {
          panov_terms <- rbind(panov_terms,
                               run_adonis(formula = paste(color_group, "*", shape_group), method = 'terms'))
        }
        
        if (is_valid_group(color_group, df_meta) && is_valid_group(filter1_group, df_meta)) {
          panov_terms <- rbind(panov_terms,
                               run_adonis(formula = paste(color_group, "*", filter1_group), method = 'terms'))
        }
        
        if (is_valid_group(shape_group, df_meta) && is_valid_group(filter1_group, df_meta)) {
          panov_terms <- rbind(panov_terms,
                               run_adonis(formula = paste(shape_group, "*", filter1_group), method = 'terms'))
        }
        
        
        panov_onedf <- if (length(perm_group) > 0)
          run_adonis(formula = paste(perm_group, collapse = " + "), method = 'onedf')
        else
          NULL
        
        df_permanova1 <- dplyr::bind_rows(panov_margin, panov_terms, panov_onedf) %>% 
          mutate(!!plot_group := pg)
        
        
        # Similarly for PERMDISP
        disp_color <- if (is_valid_group(color_group, df_meta)) run_permdisp(dist_matrix, color_group) else NULL
        disp_shape <- if (is_valid_group(shape_group, df_meta)) run_permdisp(dist_matrix, shape_group) else NULL
        disp_filter1 <- if (is_valid_group(filter1_group, df_meta)) run_permdisp(dist_matrix, filter1_group) else NULL
        
        df_permdisp1 <- list(disp_color, disp_shape, disp_filter1) %>%
          compact() %>%
          bind_rows() %>%
          mutate(!!plot_group := pg)
        
        print(df_permdisp1)
        print(df_permanova1)
        
        df_permanova=bind_rows(df_permanova,df_permanova1)
        df_permdisp=bind_rows(df_permdisp,df_permdisp1)
        #############################################################################
        print('396 # permANOVA as NMDS END')
        #############################################################################
        
      }# end plot group data loop
      
      ################################################################################
      print('147 # optional loadings not done yet')
      ################################################################################
      if(exists('use_loadings')){
        loadings_long_df <- df_NMDS_loadings %>%
          rename(Feature = Species) %>%
          pivot_longer(cols = starts_with("x_") | starts_with("y_")| starts_with("R_"),
                       names_to = c(".value", "model","data_type"),
                       names_sep = "_")%>%
          rename(normalized_R='R')
        
        df_feature_labels=feature_labels()
        
        plot_loadings=loadings_long_df%>%
          mutate(ASV=Feature)%>%
          mutate(feature=Feature)%>%
          group_by(model,data_type)%>%
          mutate(x_scale_adjust=max(x)-min(x))%>%
          mutate(y_scale_adjust=max(y)-min(y))%>%
          mutate(scale_adjust=max(x_scale_adjust,y_scale_adjust))%>%
          mutate(scale_adjust=1)%>%
          mutate(x_adj=x*scale_adjust)%>%
          mutate(y_adj=y*scale_adjust)%>%
          left_join(df_feature_labels)%>%
          mutate(data_type=factor(data_type,levels=c('bray','normalized','pseudo','euclidean','clr')))
        
        if(taxa_levs!='Phylum'){  
          
          plot_loadings=loadings_long_df%>%
            mutate(contribution=sqrt((x*x+y*y)))#%>%
          mutate(ASV=Feature)%>%
            mutate(feature=Feature)%>%
            group_by(model,data_type)%>%
            arrange(desc(contribution))%>%
            slice(1:10)%>%
            mutate(x_scale_adjust=max(x)-min(x))%>%
            mutate(y_scale_adjust=max(y)-min(y))%>%
            mutate(scale_adjust=max(x_scale_adjust,y_scale_adjust))%>%
            mutate(scale_adjust=1)%>%
            mutate(x_adj=x*scale_adjust)%>%
            mutate(y_adj=y*scale_adjust)%>%
            left_join(df_feature_labels)#%>%
          mutate(ASV_label4=sapply(str_split(ASV, ";"), function(x) x[6]))%>%
            mutate(ASV_label3=paste0(ASV_label4,'_',rel_abun_label))%>%
            mutate(data_type=factor(data_type,levels=c('bray','normalized','pseudo','euclidean','clr')))
        }   
        
        p=q+geom_segment(aes(x = 0, y = 0, xend = x_adj, yend = y_adj),
                         data=plot_loadings, arrow = arrow(length = unit(0.2, "cm")), linewidth = .5,color='gray30',
                         alpha=.5)+
          geom_text(data=plot_loadings,aes(x=x_adj*1.1,y=y_adj*1.1,label = feature_label3,color=Domain),show.legend = FALSE,alpha=.8,size=2)
        p
        
        
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
      
      theme_common<-theme_global(base_size = 11)+
        theme(      
          axis.title.x = element_markdown(size = rel(1.5),margin = margin(t = 10)),
          axis.title.y = element_markdown(size = rel(1.5),margin = margin(r = 10)),
          panel.grid.major = element_line(color='gray85',size=.25),
          plot.caption = element_text(size = rel(1.2) , face = 'bold',hjust=0),
          plot.margin = margin(margin_size, margin_size, margin_size, margin_size)
        )
      
      
      
      gPlot <- function(p) {
        p=p+
          scale_color_manual(values = palette_color,labels=palette_label)+
          scale_size_manual(values = palette_size,labels=palette_label)+
          scale_shape_manual(values = palette_shape,labels=palette_label)+
          scale_fill_manual(values = palette_color,labels=palette_label)+
          theme_common
        print(p)
        p
      }
      tPlot <- function(p) {
        p=p+
          theme_common+
          theme(
            strip.background = element_rect(fill = "gray40", color = "black"),
            strip.text = element_text(color = "white",face='bold'),
            panel.grid.major = element_blank(),  # remove major gridlines
            panel.grid.minor = element_blank(),   # remove minor gridlines
            panel.spacing = unit(0, "pt"),
            #panel.background = element_rect(fill = "gray40", color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x  = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_blank(),
            legend.position = 'none')+
          scale_fill_manual(breaks=breaks,values = values, guide = "none") +
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0))
        print(p)
        p  
      }
      
      
      ################################################################################
      print('416 # pre loop initializing')
      ################################################################################
      plot_list =plot_list_anova=plot_list_disp<- list()
      
      ################################################################################
      print('421 # Begin plotting group loop for plotting')
      ################################################################################
      disp_p_label=disp_f_label=p_label=f_label=list()
      
      df_permdisp1=df_permdisp%>%
        rename('SumOfSqs'='Sum Sq','F_value'='F','p_value'='Pr(>F)','formula'='grouping_var') %>% 
        mutate(Term=case_when(!Term %in% c('Residual','Total','Residuals')~formula,
                              TRUE~Term))
      
      df_tests <- bind_rows(df_permdisp1,df_permanova) %>% 
        select(experiment,test,method,formula,Term,Df,R2,F_value,p_value,everything())
      df_tests
      
      for (pg in plot_groups){print(pg)
        
        ################################################################################
        print('529 # Annotations as NMDS begin')
        ################################################################################
        
        col_order=c('test','method','formula','Term','average_dist_median','R2','Df','F_value','p_value')
        
        df_tests1=df_tests %>%   
          filter(.data[[plot_group]] == pg)%>%
          mutate(rownames = row_number())%>%
          filter(!Term %in% c('Residual','Total','Residuals'))%>% 
          #ungroup()
          # mutate(Term = ifelse(
          #   Term == "locationCrew Sleeping Quarters (subbed for Laundry area)",
          #   "locationCrew Sleeping Quarters",
          #   Term
          # )) %>% 
          mutate(F_value=round(F_value,2))%>%
          mutate(R2=round(R2,2))%>%
          mutate(p_value=round(p_value,4))%>%
          mutate(short_name=paste(test,method,formula,sep='_')) %>% 
          mutate(test_name=paste(short_name,rownames,sep='_')) %>% 
          mutate(across(everything(), as.character)) %>% 
          select(any_of(col_order),everything()) %>% 
          ungroup()
        head(df_tests1)
        unique(df_tests1$test)
        
        names(df_tests1)
        
        tests_text=df_tests1 %>% 
          #   filter(is.na(by) | by != 'onedf') %>% 
          select(any_of(col_order))%>%
          mutate(
            label = ifelse(
              test == 'permDISP',
              paste0(test, ' ',testing,' (',formula,'):   ADM=(',average_dist_median,') Df=',Df,' F=',F_value,' p=',p_value),
              paste0(test, ' ',testing,' (',formula,') ',Term,':   R2=',R2,' Df=',Df,' F=',F_value,' p=',p_value)
            )
          ) %>% 
          pull(label) %>%
          paste(collapse = "\n")
        ################################################################################
        print('450 # Annotations as NMDS end')
        ################################################################################
        ################################################################################
        print('# 497 perm test plot')
        ################################################################################
        
        # df_permdisp1=df_permdisp%>%
        #   rename('SumOfSqs'='Sum Sq','F_value'='F','p_value'='Pr(>F)','formula'='grouping_var')
        
        names(df_tests1)
        col_order=c('method','formula','Term','average_dist_median','R2','Df','F_value','p_value')
        col_disp=c('formula','Term','average_dist_median','R2','Df','F_value','p_value')
        col_anova=c('method','formula','Term','R2','Df','F_value','p_value')
        
        # Example data
        df_tests_plot1 <- df_tests1 %>% 
          mutate(average_dist_median = str_replace_all(average_dist_median, "___", "\n")) %>% 
          mutate(formula = paste0('(~',formula,')')) %>% 
          pivot_longer(any_of(col_order), names_to = "col", values_to = "val")%>%
          filter(!is.na(val)) %>% 
          mutate(widths = case_when(
            col == "average_dist_median" ~ 10,
            col == "formula" ~ 5,
            col == "Term" ~ 5,
            col == "method" ~ 2,
            col == "test" ~ 2,
            TRUE~2
          ))%>%
          
          
          mutate(col = factor(col, levels = col_order))%>%
          mutate(test_name = factor(test_name, levels = unique(test_name))) %>% 
          group_by(short_name) %>%
          mutate(group_color = as.integer(as.factor(short_name)) %% 2 + 1)
        
        unique(df_tests_plot1$col)
        names(df_tests_plot1)
        names(df_tests)
        
        alt_colors <- c("gray60", "gray80")
        alt_colors <- rep(c("gray60", "gray80"), 10)
        ################################################################################
        print('540 # first plot  permDISP')
        ################################################################################
        
        df_disp_plot=df_tests_plot1 %>% 
          filter(test=='permDISP') %>% 
          filter(col %in% col_disp) %>% 
          ungroup() %>% 
          mutate(short_name=factor(short_name,levels=unique(short_name)))
        
        
        breaks=as.character(unique(df_disp_plot$short_name))
        values=alt_colors[1:length(breaks)]
        names(values)=breaks
        
        p=ggplot(df_disp_plot, aes(x = col, y = short_name,fill=short_name)) +
          geom_tile(aes(width = widths
                        #  ,height=heights
          ), 
          #fill = "white",
          color = "black") +
          geom_text(aes(label = paste(val))) +
          facet_grid2(test ~ col, scales = "free", space = "free"
                      #,independent = "all"
          )
        
        plot_disp=tPlot(p)
        
        plot_list_disp[[pg]] <- plot_disp
        print(plot_disp)
        ggsave('test.png',plot_disp,width=12,height=8)
        
        ################################################################################
        print('540 # repeat for permANOVA')
        ################################################################################
        df_anova_plot=df_tests_plot1 %>% 
          filter(test=='permANOVA') %>% 
          filter(col %in% col_anova) %>% 
          #filter(!col %in% c('method','Term')) %>% 
          ungroup() %>% 
          mutate(short_name=factor(short_name,levels=unique(short_name)))
        
        
        breaks=as.character(unique(df_anova_plot$short_name))
        values=alt_colors[1:length(breaks)]
        names(values)=breaks
        
        breaks
        values
        
        
        p=ggplot(df_anova_plot, aes(x = col, y = rev(test_name),fill=short_name)) +
          geom_tile(aes(width = widths
                        #  ,height=heights
          ), 
          #fill = "white",
          color = "black") +
          geom_text(aes(label = paste(val))) +
          facet_grid2(test ~ col, scales = "free", space = "free"
                      #,independent = "all"
          )
        plot_anova=tPlot(p)
        
        plot_list_anova[[pg]] <- plot_anova
        print(plot_anova)
        ################################################################################
        print('540 # Annotations as NMDS end')
        ################################################################################
        ################################################################################
        print('387 # df_plot')
        ################################################################################
        df=df_nmds%>%
          left_join(meta)%>%
          left_join(df_richness)%>%
          mutate(richness_size=case_when(
            smp_richness>80~'g80',
            smp_richness>40~'g40',
            smp_richness>20~'g20',
            smp_richness>10~'g10',
            smp_richness>0~'g0',
          )) 
        
        
        stress=df_nmds%>%
          filter(.data[[plot_group]] == pg) %>%
          distinct(stress)%>%
          pull(stress)
        
        df_plot <- df%>%
          select(sample_name,experiment,location,treatment,sample_time,x_nmds_bray,y_nmds_bray,Bacterial_dPCR,Fungal_dPCR,smp_richness,richness_size) %>% 
          mutate(log10_Bacterial_dPCR=log10(Bacterial_dPCR)) %>% 
          mutate(log10_Bacterial_dPCR_size=case_when(
            log10_Bacterial_dPCR>10~'log10',
            log10_Bacterial_dPCR>8~'log8',
            log10_Bacterial_dPCR>6~'log6',
            log10_Bacterial_dPCR>4~'log4',
            log10_Bacterial_dPCR>0~'log0',
          )) %>% 
          pivot_longer(cols = c('richness_size','log10_Bacterial_dPCR_size'),names_to ='method' ,values_to = 'size') %>% 
          # mutate(log10_Bacterial_dPCR_size=factor(log10_Bacterial_dPCR_size, levels = c('log0','log4','log6','log8','log10')))%>%
          mutate(!!sym(filter1_group):=factor(!!sym(filter1_group), levels = filter1_group_order))%>%
          mutate(!!sym(shape_group):=factor(!!sym(shape_group), levels = shape_group_order))%>%
          mutate(!!sym(color_group):=factor(!!sym(color_group), levels = color_group_order))%>%
          rename(x=x_nmds_bray,y=y_nmds_bray )
        names(df_plot)
        
        ################################################################################
        ################################################################################
        
        
        
        
        ###############################################################################
        message('607 # strip backgrounds with outlines')
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
          # scale_size(range = c(1, 12))+  
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
            # size = "none",
            color = guide_legend(title='Ellipses',override.aes = list(linewidth=1.5)), 
            #fill = guide_legend(title='Point Color',override.aes = list(shape=22,size=5)), 
            fill = guide_legend(title='',override.aes = list(shape=22,size=5)), 
            shape = guide_legend(title='Point Shape',override.aes = list(size=5)), 
            
          )+
          labs(
            y=paste0("NMDS2; <span style='color:transparent'>_</span> k = ",k_value),
            x=paste0("NMDS1; <span style='color:transparent'>_</span> stress = ",round(stress,2)),
            title = plot_title,
            fill='',
            color='',
          )+
          theme(
            # Add border around the panel
            panel.border = element_rect(color = "gray30", fill = NA, size = 1),
            
            # Add border around the strip background
            strip.background = element_rect(color = "gray30", size = 1),
            legend.position = 'bottom'
            
          )
        p
        plot_sample_time=gPlot(p)+theme(legend.position = 'bottom')
        plot_sample_time
        
        # ggsave(paste0(output_plot,plot_title,'.png'),plot_sample_time,width=36,height=10)
        ggsave(paste0(output_plot,'Figure_4_',plot_title,'_ILMAH_',current_date,'.png'),plot_sample_time,width=36,height=10)
        ################################################################################ 
        ################################################################################
        print('# plot sample_time_richness')
        plot_title=paste(p_title,pg,'sample_time_richness',sep='_')
        plot_title
        ################################################################################
        p=ggplot(df_plot %>% filter(method=='richness_size'), aes(x = x, y = y)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_point(aes(fill = !!sym(filter1_group),size=size
                         #shape = !!sym(shape_group),size = !!sym(shape_group)
          ),
          shape=21
          )+  
          #scale_size(range = c(1, 12))+  
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
            # size = "none",
            color = guide_legend(title='Ellipses',override.aes = list(linewidth=1.5)), 
            fill = guide_legend(title='',override.aes = list(shape=22,size=5),order=1,nrow=1), 
            shape = guide_legend(title='Point Shape',override.aes = list(size=5)), 
            size = guide_legend(title='Richness >',order=2,nrow=1), 
          )+
          labs(
            y=paste0("NMDS2; <span style='color:transparent'>_</span> k = ",k_value),
            x=paste0("NMDS1; <span style='color:transparent'>_</span> stress = ",round(stress,2)),
            title = plot_title,
            fill='',
            color='',
          )+
          theme(
            # Add border around the panel
            panel.border = element_rect(color = "gray30", fill = NA, size = 1),
            
            # Add border around the strip background
            strip.background = element_rect(color = "gray30", size = 1),
            legend.position = 'bottom',
            legend.box = "vertical"
            
          )
        p
        gPlot(p)
        plot_sample_time=gPlot(p)+theme(legend.position = 'bottom',legend.box = "vertical")
        plot_sample_time
        # ggsave(paste0(output_plot,plot_title,'.png'),plot_sample_time,width=36,height=10)
        ggsave(paste0(output_plot,'Figure_4_',plot_title,'_ILMAH_',current_date,'.png'),plot_sample_time,width=36,height=10)
        ################################################################################ 
        ################################################################################
        print('# plot sample_time_dPCR')
        plot_title=paste(p_title,pg,'sample_time_dPCR',sep='_')
        plot_title
        ################################################################################
        p=ggplot(df_plot %>% filter(method=='log10_Bacterial_dPCR_size'), aes(x = x, y = y)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_point(aes(fill = !!sym(filter1_group),size=size,
                         #shape = !!sym(shape_group),size = !!sym(shape_group)
          ),
          shape=21
          )+  
          #scale_size(range = c(1, 12))+  
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
            # size = "none",
            color = guide_legend(title='Ellipses',override.aes = list(linewidth=1.5)), 
            fill = guide_legend(title='',override.aes = list(shape=22,size=5),order=1,nrow=1), 
            shape = guide_legend(title='Point Shape',override.aes = list(size=5)), 
            size = guide_legend(title='Bacterial dPCR >',order=2,nrow=1), 
          )+
          labs(
            y=paste0("NMDS2; <span style='color:transparent'>_</span> k = ",k_value),
            x=paste0("NMDS1; <span style='color:transparent'>_</span> stress = ",round(stress,2)),
            title = plot_title,
            fill='',
            color='',
          )+
          theme(
            # Add border around the panel
            panel.border = element_rect(color = "gray30", fill = NA, size = 1),
            
            # Add border around the strip background
            strip.background = element_rect(color = "gray30", size = 1),
            legend.position = 'bottom'
            
          )
        p
        plot_sample_time=gPlot(p)+theme(legend.position = 'bottom')
        # ggsave(paste0(output_plot,plot_title,'.png'),plot_sample_time,width=36,height=10)
        ggsave(paste0(output_plot,'Figure_4_',plot_title,'_ILMAH_',current_date,'.png'),plot_sample_time,width=36,height=10)
        ################################################################################ 
        print('# plot normal plot')
        plot_title=paste(p_title,pg,sep='_')
        plot_title
        ###############################################################################
        
        
        p=ggplot(df_plot, aes(x = x, y = y)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_point(aes(fill = !!sym(color_group),shape = !!sym(shape_group),size = !!sym(shape_group)))+  
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
          labs(
            y=paste0("NMDS2; <span style='color:transparent'>_</span> k = ",k_value),
            x=paste0("NMDS1; <span style='color:transparent'>_</span> stress = ",round(stress,2)),
            #caption=tests_text,
            title=plot_title
          )+
          guides(
            size = "none",
            color = guide_legend(title='Ellipses',override.aes = list(linewidth=1.5)), 
            fill = guide_legend(title='Point Color',override.aes = list(shape=22,size=5)), 
            shape = guide_legend(title='Point Shape',override.aes = list(size=5)), 
          )+
          ggtitle(plot_title)
        
        # p
        # if (shape_group %in% colnames(df) && length(unique(df[[shape_group]])) > 1) {
        #   p <- p + stat_ellipse(aes(group = !!sym(shape_group), color = !!sym(shape_group)), linewidth = 2, level = 0.95)
        # }
        # 
        # if (color_group %in% colnames(df) && length(unique(df[[color_group]])) > 1) {
        #   p <- p + stat_ellipse(aes(group = !!sym(color_group), color = !!sym(color_group)), linewidth = 1, level = 0.95)
        # }
        # if (filter1_group %in% colnames(df) && length(unique(df[[filter1_group]])) > 1) {
        #   p <- p + stat_ellipse(aes(group = !!sym(filter1_group), color = !!sym(filter1_group)), linewidth = 1, level = 0.95)
        # }
        
        p
        gPlot(p)
        plot_list[[pg]] <- gPlot(p)
        
      }
      
      
      # ############################################################################
      # print('# save plots')
      # ############################################################################
      
      # Create a list of combined plots (each combined plot is one object)
      combined_plot_list <- Map(function(top,middle,bottom) {
        wrap_elements(top) / wrap_elements(middle)/ wrap_elements(bottom) +
          plot_layout(heights = c(5,3,4))
      }, plot_list, plot_list_disp, plot_list_anova)
      
      combined_plot_list
      
      
      combined_test_list <- Map(function(top,bottom) {
        wrap_elements(top) / wrap_elements(bottom) +
          plot_layout(heights = c(1,2))
      },plot_list_disp, plot_list_anova)
      
      combined_test_list
      
      plot_width=12*length(combined_plot_list)/max(1,floor(length(combined_plot_list)/2))
      plot_height=12*max(1,floor(length(combined_plot_list)/2))
      
      patch_combined_plots <- wrap_plots(combined_plot_list)
      print(patch_plots)
      ggsave(paste0(output_plot,p_title,'.png'),patch_plots,width=plot_width+4,height=plot_height+8)
      
      
      patch_test_plots <- wrap_plots(combined_test_list)
      print(patch_test_plots)
      
      patch_plots <- wrap_plots(plot_list)
      print(patch_plots)
      
      ggsave(paste0(output_plot,'Figure_4_',p_title,'_combined_ILMAH_',current_date,'.png'),patch_combined_plots,width=plot_width,height=plot_height*2)
      
      ggsave(paste0(output_plot,'Figure_4_',p_title,'_TESTS_ILMAH_',current_date,'.png'),patch_test_plots,width=plot_width,height=plot_height)
      
      ggsave(paste0(output_plot,'Figure_4_',plot_title,'_ILMAH_',current_date,'.png'),patch_plots,width=plot_width+4,height=plot_height)
      
      write.csv(df_plot,paste0(output_plot,'Figure_4_',plot_title,'_ILMAH_',current_date,'.csv'),row.names=FALSE)
      write.csv(df_tests,paste0(output_plot,'Figure_4_',plot_title,'_TESTS_ILMAH_',current_date,'.csv'),row.names=FALSE)
    }# plot group  r/ASV level
  } # r loop ASV levels # run group
  
} # end asv levels # plot set order
#}# end Run_Group/loop group 

#plot_list_disp   df_tests      print(plot_anova)
################################################################################
################################################################################

################################################################################
################################################################################




################################################################################

