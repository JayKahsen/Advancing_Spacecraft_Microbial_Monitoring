source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
###############################################################################
# Description Start line 80ish
################################################################################
script_title='figure_4'
options(scipen = 999) # don;t use scientific notation
################################################################################
# Running info and switches
################################################################################

# output_data=paste0('output_data/',script_title,'/')
# output_plot=paste0('output_plot/',script_title,'/')

################################################################################
# Extra Script Variables and functions
################################################################################

################################################################################
# Extra Script Variables and functions
################################################################################
new_labels=c('cotton'='Cotton Swab','macrofoam'='Macrofoam Swab')
palette_label=c(new_labels,palette_label)




################################################################################
# plotting settings
################################################################################
strip_text_size=12
x_text_size=10
margin_size=10

pool_height=10
pool_width=7
ind_height=10
ind_width=8

p_value_size=3
geom_text_size=3
ind_text_size=3
legend_text_size=20
axis_text_size=axis_title_size=legend_text_size+1
################################################################################
# Theme and gPlot
################################################################################
theme_common<-theme_global(base_size = 11)+
  theme(       strip.text = element_markdown(size = rel(1.2), face = 'bold'),
               legend.position = 'bottom'
    
  )
gPlot <- function(p) {
  p=p+
    theme_common+
    #scale_fill_manual(values = palette_color,labels=legend_labels) +
    #scale_fill_manual(values = palette_color,labels=palette_label) +
    scale_fill_manual(
      breaks = sampling_device_order,
      values = palette_color,
      labels = filtered_labels,
      na.translate = FALSE
    ) +
    scale_color_manual(values = palette_color,labels=palette_label)+
    scale_linetype_manual(values = palette_linetype,labels=palette_label)+
    scale_shape_manual(values = palette_shape,labels=palette_label)
  print(p)
  p  
}



################################################################################
# load files
################################################################################

figure_number='4'

experiment_type_order=c( 'metal_deposition','swab_head_retention' )

df=meta %>% 
  filter(experiment_type%in%c( 'metal_deposition','swab_head_retention')) %>% 
  filter(!is.na(percent_recovery_dpcr)) %>% 
  filter(sampling_device %in% sampling_device_order) %>% 
  select(sample_id,percent_recovery_qpcr,percent_recovery_dpcr,sampling_device,experiment_type,figure,Figure,bar_fill)

figure_4_selection=df %>% 
  mutate(figure_4='4') %>% 
  mutate(figure_4=case_when(
    experiment_type=='swab_head_retention'~'4A',
    experiment_type=='metal_deposition' & sampling_device %in% c('macrofoam','cotton')~'4B',
    sampling_device %in% c('SALSA','Polyester Wipe')~'4C',
    TRUE~NA
  ))%>%
  select(sample_id,figure_4) %>% 
  write.csv('data_tables/figure_4_selection.csv',row.names = FALSE)

df_plot <- df %>%
  filter(str_detect(Figure, figure_number)) %>% 
  pivot_longer(
    cols = c(percent_recovery_dpcr, percent_recovery_qpcr),
    names_to = "method",
    values_to = "value"
  ) %>%
  mutate(method = str_replace(method, "percent_recovery_", "")) %>% 
  # mutate(bar_type = factor(bar_type, levels = bar_type_order))%>% 
  # mutate(bar_name = factor(bar_name, levels = bar_name_order))%>% 
  # mutate(bar_color=factor(bar_color,levels=bar_color_order)) %>% 
   mutate(method=factor(method,levels=method_order)) %>% 
  ungroup()


# sampling_device_order
# vPrint(unique(df$experiment_type))

################################################################################
# make meta data
################################################################################

################################################################################
# Setup
################################################################################

################################################################################
# Prepare Data
################################################################################
################################################################################
# Libraries
################################################################################


################################################################################
plot_title='Figure_4'
################################################################################
df_plot1= df_plot%>%
#  filter(figure=='4A')%>%
  ungroup()

################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value",
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
df_summary11  <- results$summary %>% 
  left_join(meta %>% distinct(bar_fill,sampling_device))
df_segments11 <- results$segments %>% 
  mutate(x=ifelse(group1 %in% c('cotton','macrofoam'),(x-2),x)) %>% 
  mutate(xend=ifelse(group1 %in% c('cotton','macrofoam'),(xend-2),xend)) %>%
  ungroup()
df_tests11    <- results$tests %>% 
  mutate(label_position=ifelse(group1 %in% c('cotton','macrofoam'),(label_position-2),label_position)) %>%
  ungroup()

df_summary1  <- results$summary%>% 
  left_join(meta %>% distinct(bar_fill,sampling_device))
df_segments1 <- results$segments 
df_tests1    <- results$tests %>% 
  mutate(p.value=welch_raw_val) %>% 
  mutate(significance=welch_raw_sig)


filtered_labels <- palette_label[sampling_device_order]

################################################################################
message('plot_A_13')
################################################################################

df_summary=df_summary1 %>% filter(figure %in% c('4A','4B','4C'))
df_segments=df_segments1%>% filter(figure %in% c('4A','4B','4C'))
df_tests=df_tests1  %>% filter(figure %in% c('4A','4B','4C')) 

p=ggplot(df_plot1 , aes(x = sampling_device, y = value, fill = sampling_device)) +
 geom_jitter(width = 0.2, size = point_size+1, alpha = 0, shape = 21) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(method~figure,scale='free_x',
             labeller = labeller(.cols = palette_label, .rows = palette_label),
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=1.5
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_value,
      ymax = mean_value,
      color=sampling_device
    ),
    width = 0.6, size = 1.5, inherit.aes = FALSE
  )+
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width  = 0.8
    ),
    size = point_size + 1,
    alpha = 0.8,
    shape = 21
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)))+
  labs(
    fill = 'Sampling Device',
    color = '',
title = '',
    y = "Percent Recovery (%)",
    x = ''
  ) +
  guides(color = "none")

p


plot_A_13=gPlot(p)+theme(legend.position='right')
plot_A_13
################################################################################
################################################################################
message('plot_A_1')
################################################################################
df_summary=df_summary1 %>% filter(figure=='4A')
df_segments=df_segments1%>% filter(figure=='4A')
df_tests=df_tests1  %>% filter(figure=='4A')

p=ggplot(df_plot1%>% filter(figure=='4A'), aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_jitter(width = 0.2, size = point_size+1, alpha = 0, shape = 21) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(method~figure,scale='free_x',
             labeller = labeller(.cols = palette_label, .rows = palette_label),
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=1.5
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_value,
      ymax = mean_value,
      color=sampling_device
    ),
    width = 0.6, size = 1.5, inherit.aes = FALSE
  )+
  geom_jitter(width = 0.2, size = point_size+1, alpha = 0.8, shape = 21) +
  geom_errorbar(
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_minus_sd,
      ymax = mean_plus_sd,
    ),
    width = 0.1, size = 1, inherit.aes = FALSE, color=bar_color
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)))+
  labs(
    fill = 'Sampling Device',
    color = '',
    title = '',
    y = "Percent Recovery (%)",
    x = ''
  ) +
  guides(color = "none")

p


plot_A_1=gPlot(p)
################################################################################
################################################################################
message('plot_A_2')
################################################################################
df_summary=df_summary1 %>% filter(figure=='4B')
df_segments=df_segments1%>% filter(figure=='4B')
df_tests=df_tests1  %>% filter(figure=='4B')

p=ggplot(df_plot1%>% filter(figure=='4B'), aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_jitter(
    width = 0.2,      # horizontal jitter
    height = 0,       # NO vertical jitter
    size = point_size + 1,
    alpha = 0,
    shape = 21
  )+

  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(method~figure,scale='free_x',
             labeller = labeller(.cols = palette_label, .rows = palette_label),
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=1.5
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_value,
      ymax = mean_value,
      color=sampling_device
    ),
    width = 0.6, size = 1.5, inherit.aes = FALSE
  )+
  geom_point(
    position = position_jitter(
      width = 0.2,   # horizontal jitter
      height = 0     # NO vertical jitter
    ),
    size = point_size + 1,
    alpha = 0,
    shape = 21
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)))+
  labs(
    fill = 'Sampling Device',
    color = '',
    title = '',
    y = "Percent Recovery (%)",
    x = ''
  ) +
  guides(color = "none")

p


plot_A_2=gPlot(p)
################################################################################
################################################################################
message('plot_A_3')
################################################################################
df_summary=df_summary1 %>% filter(figure=='4C')
df_segments=df_segments1%>% filter(figure=='4C')
df_tests=df_tests1  %>% filter(figure=='4C')

p=ggplot(df_plot1%>% filter(figure=='4C'), aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_jitter(width = 0.2, size = point_size+1, alpha = 0, shape = 21) +
  geom_point(
    position = position_jitter(
      width = 0.2,   # horizontal jitter
      height = 0     # NO vertical jitter
    ),
    size = point_size + 1,
    alpha = 0,
    shape = 21
  )+
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_grid(method~figure,scale='free_x',
             labeller = labeller(.cols = palette_label, .rows = palette_label),
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=1.5
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_value,
      ymax = mean_value,
      color=sampling_device
    ),
    width = 0.6, size = 1.5, inherit.aes = FALSE
  )+
  geom_jitter(width = 0.2, size = point_size+1, alpha = 0.8, shape = 21) +
  geom_errorbar(
    data = df_summary,
    aes(
      x = sampling_device,
      ymin = mean_minus_sd,
      ymax = mean_plus_sd,
    ),
    width = 0.1, size = 1, inherit.aes = FALSE, color=bar_color
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)))+
  labs(
    fill = 'Sampling Device',
    color = '',
    title = '',
    y = "Percent Recovery (%)",
    x = ''
  ) +
  guides(color = "none")

p


plot_A_3=gPlot(p)
################################################################################

write.csv(df,paste0(output_plot,plot_title,'_PPR_',current_date,'.csv'))
#ggsave(paste0(output_plot,plot_title,'_PPR_',current_date,'.png'),width=15,height=10)

################################################################################
skip='yes'
if(skip=='no'){
################################################################################
plot_title='Percent Recovery 2'
################################################################################
df_plot1= df_plot%>%
  filter(figure=='4B')

################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = NA,
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
df_tests    <- results$tests

################################################################################
p=ggplot(df_plot1, aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_col(data = df_summary,aes(x = sampling_device,y= mean_value))+
  geom_jitter(width = 0.2, size = point_size, alpha = 0.8, shape = 21) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=1.5
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
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)))+
  
  
  # coord_cartesian(ylim = zoom)+
  #  scale_y_log10(
  #    labels = log_labels_bold,
  #    breaks = log_breaks,
  #    limits = c(0.1, 1e6)
  #  )+
  
  labs(
    fill = 'Sampling Device',
    color = '',
    title = 'Microbial DNA Recovery Efficiency<br>
    of Swabs on Spiked Coupons (25 cm<sup>2</sup>)',
    y = "Percent Recovery (%)",
    x = ''
  ) 

p

plotB=gPlot(p)
plotB
################################################################################
################################################################################
plot_title='Percent Recovery 3'
################################################################################
df_plot1= df_plot%>%
  filter(figure=='4C')

################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "sampling_device",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = NA,
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
df_tests    <- results$tests

################################################################################
p=ggplot(df_plot1, aes(x = sampling_device, y = value, fill = sampling_device)) +
  geom_col(data = df_summary,aes(x = sampling_device,y= mean_value))+
  geom_jitter(width = 0.2, size = point_size, alpha = 0.8, shape = 21) +
  geom_segment(
    data = df_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    # linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3.2,vjust=1.5
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
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)))+

  
  
  # coord_cartesian(ylim = zoom)+
  #  scale_y_log10(
  #    labels = log_labels_bold,
  #    breaks = log_breaks,
  #    limits = c(0.1, 1e6)
  #  )+
  
  labs(
    fill = 'Sampling Device',
    color = '',
    title = 'Microbial DNA Recovery Efficiency<br>
    of Swabs on Spiked Coupons (25 cm<sup>2</sup>)',
    y = "Percent Recovery (%)",
    x = ''
  ) 


p

plotC=gPlot(p)
plotC
################################################################################
final_plot <- (plotA  | plotB | plotC) +
  #plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = 'A',
    tag_suffix = '.',
theme_plot
  )


plot_title=script_title

final_plot <- (plotA  | plotB | plotC) +
 # plot_layout(guides = "collect") +
  plot_annotation(
    title = plot_title,
    tag_levels = 'A',
    tag_suffix = '.',
    theme = theme_plot
   #   legend.position='bottom'
)
final_plot


ggsave(paste0(output_plot,plot_title,'_PPR_',current_date,'.png'),final_plot,width=15,height=10)
write.csv(df,paste0(output_plot,plot_title,'_PPR_',current_date,'.csv'))
}# end skip
################################################################################
################################################################################
################################################################################
sc=3
st_plot=readRDS(paste0(output_plot,'Figure_',figure_number,'_st_plot.rds'))

st_plot_alt1=readRDS(paste0(output_plot,'Figure_',figure_number,'_st_plot1.rds'))
st_plot_alt2=readRDS(paste0(output_plot,'Figure_',figure_number,'_st_plot2.rds'))
st_legend1=readRDS(paste0(output_plot,'Figure_',figure_number,'_legend_plot.rds'))
st_right1=readRDS(paste0(output_plot,'Figure_',figure_number,'_legend_plot_right.rds'))
#mn_plot1=readRDS(paste0(output_plot,'Figure-',figure_number,'_mn_plot.rds'))

theme_custom=theme(
  axis.text.y = element_markdown(size = rel(1.2), face = 'bold'),
  legend.text = element_text(size =rel(1.5)),
  legend.title = element_markdown(size =rel(1.5)),
  strip.text = element_markdown(size = rel(1.4), face = 'bold'),
  axis.title.x = element_markdown(size = rel(1.2), face = 'bold',color='black'),
  axis.text = element_markdown(size = rel(1.2), face = 'bold')
  
)

st_legend=st_legend1+
  theme(plot.margin=margin(t=0,b=0))

st_right=st_right1

st_plot_alt1
st_plot1=st_plot_alt1+
 # ggtitle('B) Percent Relative Abundance')+
  theme_custom+
  theme(
    plot.margin=margin(l=10),
    #legend.key.size = unit(1.2, "cm"),
   #   axis.text.y = element_blank(),
  #   axis.ticks.y = element_blank(),
     legend.position='bottom'
    #axis.text.y = element_markdown(size = rel(1.2), face = 'bold')
  #  plot.margin=margin(l=55)
  )
st_plot1

st_plot_alt2
st_plot2=st_plot_alt2+
  # ggtitle('B) Percent Relative Abundance')+
  theme_custom+
  theme(
    plot.margin=margin(l=10),
    #legend.key.size = unit(1.2, "cm"),
    #  axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    #axis.text.y = element_markdown(size = rel(1.2), face = 'bold')
   # plot.margin=margin(l=55)
  )+
  guides(color='none',
         fill='none')
st_plot2
################################################################################
library(cowplot)
plot_A_13=plot_A_13+theme_custom
# Extract the legend from plot_A_13
legend_A13 <- get_legend(plot_A_13 + theme(legend.position = "right"))
legend_A13

# Wrap as a plot
legend_plot <- ggdraw(legend_A13)

legend_plot 


################################################################################
plot_A_1

plotA1=plot_A_1 +
  ggtitle('')+
  labs(y='% Released')+
  scale_x_discrete(labels = palette_label)+
  theme_custom+
  theme(legend.position = 'none',
        #legend.margin=margin(,l=105,r=100)
        )
plotA1
################################################################################
################################################################################
plot_A_2
  
  plotA2=plot_A_2 +
  ggtitle('')+
    labs(y='% Recovered')+
  scale_x_discrete(labels = palette_label)+
  theme_custom+
  theme(legend.position = 'none')
plotA2
################################################################################
################################################################################
plot_A_3

plotA3=plot_A_3 +
  ggtitle('')+
  labs(y='% Recovered')+
  scale_x_discrete(labels = palette_label)+
  theme_custom+
  theme(legend.position = 'none',
        
        )
plotA3
################################################################################

row1=plotA1 | legend_plot
row1

row2=plotA2 | plotA3
row2



col1=row1/row2
col1


col2 =st_plot1 / st_plot2 
col2
# +
#   plot_layout(guides = "collect") 
# #col2=wrap_elements(col2)
# col2
# col2 =(st_plot1 | st_plot2)/st_legend +
#   plot_layout(heights = c(30, 1))
# col2


combined <- (col1 | col2)+
  plot_layout(widths = c(1.2, 1))

combined
plot_height=12
plot_width=20

ggsave(paste0(output_plot,'Figure_',figure_number,'_PPR_combined_',current_date,'.pdf'),combined,width=plot_width,height=plot_height,dpi=dpi)

ggsave(paste0(output_plot,'Figure_',figure_number,'_PPR_combined_',current_date,'.png'),combined,width=plot_width,height=plot_height,dpi=dpi)

plot_A_13=plot_A_13+
  theme(plot.margin = margin())
# +
#   theme(legend.position = 'bottom')
#+ labs(color='',fill='')

st_plot=st_plot+
  theme(legend.position = 'none',
        plot.margin=margin(l=50,r=80))

row2=wrap_elements(st_plot)

row1=wrap_elements(plot_A_13)

combined_alt=row1/row2
combined_alt
ggsave(paste0(output_plot,'Figure_',figure_number,'_PPR_combined_alt_',current_date,'.pdf'),combined_alt,width=plot_width,height=plot_height,dpi=dpi)

ggsave(paste0(output_plot,'Figure_',figure_number,'_PPR_combined_alt_',current_date,'.png'),combined_alt,width=plot_width,height=plot_height,dpi=dpi)






