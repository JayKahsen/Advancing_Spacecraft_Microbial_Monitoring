source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
###############################################################################
# Description Start line 80ish
################################################################################
script_title='Figure_5'
options(scipen = 999) # don;t use scientific notation

figure_number='5'
################################################################################
# Running info and switches
################################################################################

# output_data=paste0('output_data/',script_title,'/')
# output_plot=paste0('output_plot/',script_title,'/')

################################################################################
# Extra Script Variables and functions
################################################################################
################################################################################
# Theme and gPlot
################################################################################
theme_common<-theme_global(base_size = 11)+
  theme(      legend.position='none' 
    
  )
gPlot <- function(p) {
  p=p+
    theme_common+
    #scale_fill_manual(values = palette_color,labels=legend_labels) +
    scale_fill_manual(values = palette_color,labels=palette_label) +
    scale_color_manual(values = palette_color,labels=palette_label)+
    scale_linetype_manual(values = palette_linetype,labels=palette_label)+
    scale_shape_manual(values = palette_shape,labels=palette_label)+
    scale_y_continuous(labels = log10_labels_bold,expand = expansion(mult = c(0, 0.1)))+
    scale_x_discrete(labels = palette_label)

  
  print(p)
  p  
}

################################################################################
# Extra Script Variables and functions
################################################################################
################################################################################
# plotting settings
################################################################################

################################################################################
# load files
################################################################################
file_path=paste0('data_tables/',script_title,'/')
library(stringr)

df=meta %>% 
  filter(str_detect(Figure, figure_number)) %>% 
  select(unique_replicate_id,sample_name,sample_id,copies_original_qpcr_log10,copies_original_dpcr_log10,copies_original_qpcr,copies_original_dpcr,sampling_device,experiment_type,figure,Figure,sample_reads)


figure_5_selection=df %>% 
  mutate(figure_5='5') %>% 
  mutate(figure_5=case_when(
    sampling_device %in% c('cotton','macrofoam')~'5A',
    sampling_device %in% c('SALSA','Polyester Wipe')~'5B',
    TRUE~'5'
  ))%>%
  select(sample_id,figure_5) %>% 
  write.csv('data_tables/figure_5_selection.csv',row.names = FALSE)

df_non_log10 <- df %>%
  mutate(copies_original_dpcr=as.numeric(copies_original_dpcr)) %>% 
  mutate(copies_original_qpcr=as.numeric(copies_original_qpcr)) %>% 
  pivot_longer(
    cols = c(copies_original_qpcr, copies_original_dpcr),
    # cols = c(copies_original_qpcr, copies_original_dpcr),
    names_to = "method",
    values_to = "value_non_log10"
  ) %>%
  mutate(method = str_replace(method, "copies_original_", "")) %>% 
  #mutate(method = str_replace(method, "_log10", "")) %>%
  select(sample_name,method,value_non_log10,) %>% 
  ungroup() %>% 
  mutate(method=factor(method,levels=method_order)) %>% 
  distinct()



df_plot <- df %>%
  pivot_longer(
   cols = c(copies_original_qpcr_log10, copies_original_dpcr_log10),
   # cols = c(copies_original_qpcr, copies_original_dpcr),
    names_to = "method",
    values_to = "value"
  ) %>%
  mutate(method = str_replace(method, "copies_original_", "")) %>% 
  mutate(method = str_replace(method, "_log10", "")) %>% 
  distinct() %>% 
  left_join(df_non_log10,by=c('sample_name','method')) %>% 
  ungroup() %>% 
  mutate(method=factor(method,levels=method_order)) %>% 
  distinct()
################################################################################
# Setup
################################################################################


plot_title <- 'cotton_comparison'

################################################################################
# Prepare Data
################################################################################
################################################################################
# Libraries
################################################################################


################################################################################
plot_title='cotton1'
################################################################################
df_plot99=df_plot %>% 
  filter(figure %in% c('5A','5B'))

results <- summary_brackets(
  df = df_plot99,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value_non_log10",
  facet_columns = c('figure','method'),
  testing_col= "sampling_device",
  test_type = NULL,
  adjust_method = "BH",
  # dodge_width = 0.8,
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE,
)
df_tests3   <- results$tests





df_plot1=df_plot %>% 
  filter(figure=='5A')
################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value_non_log10",
  facet_columns = 'method',
  testing_col= "sampling_device",
  test_type = NULL,
  adjust_method = "BH",
  # dodge_width = 0.8,
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE,
)
#df_summary  <- results$summary
#df_segments <- results$segments
df_tests1    <- results$tests
################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = 'method',
  testing_col= "sampling_device",
  test_type = NULL,
  adjust_method = "BH",
  # dodge_width = 0.8,
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE,
)
df_summary  <- results$summary
df_segments <- results$segments
df_tests=df_tests1 %>% 
  left_join(df_segments %>% 
              filter(part=='top') %>% 
              select(method,y)) %>% 
  mutate(label_value=y)

df_tests=df_tests3 %>% 
  filter(figure=='5A') %>% 
  select(-figure) %>% 
  left_join(df_segments %>% 
              filter(part=='top') %>% 
              select(method,y)) %>% 
  mutate(label_value=y) %>% 
  mutate(adj.significance=mann_log_sig) %>% 
  mutate(adj.p.value=mann_log_val) 
  
################################################################################
################################################################################
################################################################################
library(ggplot2)
library(ggsignif)

p=ggplot(df_plot1, aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_jitter(width = 0.2, size = point_size, alpha = 0.8, shape = 21) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(method~.,scale='free_x',
             labeller = labeller(.cols = palette_label, .rows = palette_label),
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = adj.significance
    ),
    inherit.aes = FALSE,
    size = 4,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(adj.p.value)
    ),
    inherit.aes = FALSE,
    size = 4,vjust=1.5
  )+
  geom_errorbar(
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_minus_sd,
      ymax = mean_plus_sd,
    ),
    width = 0.1, size = 1, inherit.aes = FALSE, color=bar_color
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_value,
      ymax = mean_value,
    ),
    width = 0.4, size = 1.5, inherit.aes = FALSE, color=bar_color
  )+
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
 coord_cartesian(ylim = zoom)+
  labs(
    fill='',
    color = '',
    title = 'A) JPL-Spacecraft Assembly Facility<br>
    <span style="color:white">A) </span>(25 cm<sup>2</sup>)',
    y = "log10(16S rRNA gene copies) / 25 cm<sup>2</sup>",
    x = ''
  ) 

p
plotA=gPlot(p)
plotA


################################################################################
plot_title='SALSA1'
################################################################################

df_plot2=df_plot %>% 
  filter(figure=='5B')

################################################################################
results <- summary_brackets(
  df = df_plot2,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value_non_log10",
  facet_columns = 'method',
  testing_col= "sampling_device",
  test_type = NULL,
  adjust_method = "BH",
  # dodge_width = 0.8,
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE,
)
#df_summary  <- results$summary
#df_segments <- results$segments
df_tests1    <- results$tests
################################################################################
results <- summary_brackets(
  df = df_plot2,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = 'method',
  testing_col= "sampling_device",
  test_type = NULL,
  adjust_method = "BH",
  # dodge_width = 0.8,
  width = 0.8,
  log = FALSE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE,
)
df_summary  <- results$summary
df_segments <- results$segments
df_tests=df_tests1 %>% 
  left_join(df_segments %>% 
              filter(part=='top') %>% 
              select(method,y)) %>% 
  mutate(label_value=y)

df_tests=df_tests3 %>% 
  filter(figure=='5B') %>% 
  select(-figure) %>% 
  left_join(df_segments %>% 
              filter(part=='top') %>% 
              select(method,y)) %>% 
  mutate(label_value=y)  %>% 
  mutate(adj.significance=welch_log_sig) %>% 
  mutate(adj.p.value=welch_log_val) 
################################################################################
library(ggplot2)
library(ggsignif)

p=ggplot(df_plot2, aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_jitter(width = 0.2, size = point_size, alpha = 0.8, shape = 21) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(method~.,scale='free_x',
             labeller = labeller(.cols = palette_label, .rows = palette_label),
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = adj.significance
    ),
    inherit.aes = FALSE,
    size = 4,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(adj.p.value)
    ),
    inherit.aes = FALSE,
    size = 4,vjust=1.5
  )+
  geom_errorbar(
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_minus_sd,
      ymax = mean_plus_sd,
    ),
    width = 0.1, size = 1, inherit.aes = FALSE, color=bar_color
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_value,
      ymax = mean_value,
    ),
    width = 0.4, size = 1.5, inherit.aes = FALSE, color=bar_color
  )+
 coord_cartesian(ylim = zoom)+
  scale_y_log10(
    labels = log_labels_bold,
    breaks = log_breaks,
    limits = log_limits,
    expand = expansion(mult = c(0, 0.05))
  )+
  labs(
    #fill = 'Sampling Device:',
    fill='',
    color = '',
    title = 'B) JPL-Spacecraft Assembly Facility<br>
    <span style="color:white">B) </span>(1000 cm<sup>2</sup>)',
    y = "log10(16S rRNA gene copies) / 1000 cm<sup>2</sup>",
    x = ''
  ) 

p

plotC=gPlot(p)
plotC

################################################################################
################################################################################
plot_title='SALSA2'
################################################################################
################################################################################
plotA

# =plotA+  
#   stat_count(aes(x = sampling_device,
#                  label = after_stat(paste0("n = ", count))
#   ),
#   inherit.aes = FALSE,
#   geom = "text",
#   vjust = 0,
#   fontface = "bold",
#   size = 3.5
#   )

plotC

# =plotC+  
#   stat_count(aes(x = sampling_device,
#                  label = after_stat(paste0("n = ", count))
#   ),
#   inherit.aes = FALSE,
#   geom = "text",
#   vjust = 0,
#   fontface = "bold",
#   size = 3.5
#   )


#final_plot <- (plotA / plotB) | (plotC / plotD)

final_plot <- (plotA |plotC)
final_plot


#final_plot <- (plotA / plotB) | (plotC / plotD) +
final_plot <- (plotA |plotC)


final_plot


plot_title=script_title
ggsave(paste0(output_plot,plot_title,'_PPR_',current_date,'.pdf'),final_plot,width=12,height=8)
ggsave(paste0(output_plot,plot_title,'A_PPR_BH_corrected_',current_date,'.pdf'),plotA,width=6,height=8)
ggsave(paste0(output_plot,plot_title,'B_PPR_BH_corrected_',current_date,'.pdf'),plotC,width=6,height=8)
saveRDS(plotA, paste0(output_plot, plot_title, 'A_PPR_BH_corrected_', current_date, '.rds'))
saveRDS(plotC, paste0(output_plot, plot_title, 'B_PPR_BH_corrected_', current_date, '.rds'))

#write.csv(df,paste0(output_plot,plot_title,'_PPR_',current_date,'.csv'))
