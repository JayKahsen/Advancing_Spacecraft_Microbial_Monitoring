source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
###############################################################################
# Description Start line 80ish
################################################################################

options(scipen = 999) # don;t use scientific notation

figure_number='2'

script_title=paste0('Figure_',figure_number)
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
################################################################################
# plotting settings
################################################################################
# strip_text_size=12
# x_text_size=10
# margin_size=10
# 
# pool_height=10
# pool_width=7
# ind_height=10
# ind_width=8

theme_common<-theme_global(base_size = 11)+
  theme(       
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, 
                               #size = axis_text_size,
                               hjust=1,vjust=1),
    legend.position='bottom'
  )

# plot_width=10
# plot_height=10
# margin_buffer=30


#You can start X-axis with Timepoint 4 at 0.
#Can you please use the scale limits for the Y-axis to be from 0 – 4.


gPlot <- function(p) {
  p=p+
    scale_color_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
    scale_size_manual(values = palette_size,labels=palette_label)+
    labs(
      fill='',
      size='Sequence Reads'
      )+
    theme_common
  print(p)
}

df_sample_counts=read.csv(paste(output_plot,'df_sample_counts.csv')) %>% 
  mutate(experiment=factor(experiment, levels = experiment_order))%>%
  mutate(bar_color=factor(bar_color, levels = rev(bar_color_order)))%>%
  mutate(bar_type_2=factor(bar_type_2, levels = bar_type_2_order))%>%
  distinct()

names(meta)
################################################################################
# load files
###############################################################################
names(meta)

df=meta %>% 
  filter(str_detect(Figure, figure_number)) %>%  
  select(unique_replicate_id,sample_id,Figure,figure,bar_name,copies_rxn_dpcr,copies_rxn_qpcr,bar_color,bar_fill,bar_type_2,bar_axis_2,sample_reads)
names(meta)

figure_2_selection=df %>% 
  mutate(figure_2='2') %>% 
  select(sample_id,figure_2) %>% 
  write.csv('data_tables/figure_2_selection.csv',row.names = FALSE)

sample_read_size_order=c('ng100','g100','g1K','g10K','g100K')

df_plot <- df %>%
  mutate(sample_read_size = case_when(
    is.na(sample_reads)      ~ "ng100",
    sample_reads > 100000    ~ "g100K",
    sample_reads > 10000     ~ "g10K",
    sample_reads > 1000      ~ "g1K",
    sample_reads > 100       ~ "g100",
    TRUE                     ~ "ng100",
  )) %>% 
  mutate(
    copies_rxn_dpcr = as.numeric(copies_rxn_dpcr),
    copies_rxn_qpcr = as.numeric(copies_rxn_qpcr)
  ) %>%
  pivot_longer(
    cols = c(copies_rxn_dpcr, copies_rxn_qpcr),
    names_to = "method",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>% 
  mutate(method = str_replace(method, "copies_rxn_", "")) %>% 
  mutate(bar_type_2 = factor(bar_type_2, levels = bar_type_2_order))%>% 
  mutate(bar_name = factor(bar_name, levels = bar_name_order))%>% 
  mutate(bar_axis_2 = factor(bar_axis_2, levels = bar_axis_2_order))%>% 
  mutate(bar_color=factor(bar_color,levels=bar_color_order)) %>% 
  ################################################################################
  mutate(method=factor(method,levels=method_order)) %>% 
  ################################################################################
  mutate(sample_read_size=factor(sample_read_size,levels=sample_read_size_order))

df_summary <- df_plot%>%
  group_by(bar_name,bar_type_2,bar_color,bar_fill,bar_axis_2,method) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = sd / sqrt(n)
  ) %>%
  mutate(bar_type_2 = factor(bar_type_2, levels = bar_type_2_order))%>% 
  mutate(bar_name = factor(bar_name, levels = bar_name_order))%>% 
  mutate(bar_axis_2 = factor(bar_axis_2, levels = bar_axis_2_order))%>% 
  ################################################################################
mutate(method=factor(method,levels=method_order)) %>% 
  ################################################################################
  mutate(bar_color=factor(bar_color,levels=bar_color_order))

filtered_colors <- palette_color[bar_color]
filtered_labels <- palette_label[bar_color]

plot_title='Figure_2'

#name_order=c('Cotton Device Control' ,'Macrofoam Device Control','Cotton Environmental Control','Macrofoam Environemental Control' ,'Cotton Dry + Metal','Macrofoam Dry + Metal','Cotton Wet + Metal','Macrofoam Wet + Metal' ,'Water Only','NTC')
#################
p=ggplot(df_plot,aes(x=bar_axis_2)) +
  geom_point(aes(x=bar_axis_2,y = value,color=bar_color),inherit.aes = FALSE,shape=15,alpha=0) +
 geom_errorbar(data=df_summary,aes(ymin = mean_value, ymax = mean_value,color=bar_color), width = 0.6,linewidth=1) +
  geom_point(aes(x=bar_axis_2,y = value,fill=bar_color,size=sample_read_size), position = position_jitter(width = 0.1), inherit.aes = FALSE,color='gray18',shape=21,alpha=.8
  ) +
  geom_errorbar(data=df_summary,aes(ymin = mean_value - sd, ymax = mean_value + sd), width = 0.2,color='gray18',linewidth=1) +
  
  theme(axis.text.x = element_text(
    angle = 45, 
    #size = axis_text_size,
    hjust=1,vjust=1))+
  scale_x_discrete(labels = palette_label)+
  scale_fill_manual(
    breaks = bar_color,
    values = palette_color,
    labels = filtered_labels,
    #  labels = palette_label,
    na.translate = FALSE
  ) +
  scale_color_manual(
    #     breaks = bar_color,
    values = palette_color,
    labels = palette_label,
    na.translate = FALSE
  ) +
  scale_size_manual(
    #     breaks = bar_color,
    values = palette_size,
    labels = palette_label,
    na.translate = FALSE
  ) +
  facet_grid(method~bar_type_2,scale='free_x',
             labeller = labeller(.cols = palette_label, .rows = palette_label),
  )+
  guides(
    color = guide_legend(
      override.aes = list(shape=15,linetype=NA,alpha=1,size=8)),
    fill ='none' 
    #  guide_legend(
    #   override.aes = list(shape=21,linetype=NA,alpha=1,size=2))
  )+
  labs(title='Recovery of 16S rRNA Gene Copies from Negative Controls',
       color='Control Type',
       fill='',
       size='# of Sequence Reads',
       # caption='can shorten names',
       y='Log10 (16S dPCR copies) / reaction',
       x='Negative Controls')
p



p1=p

################################################################################
range_summary <- df_plot %>%
  group_by(bar_axis_2,method,bar_type_2) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE)
  )

p_range <- p1 +
  geom_text(
    data = range_summary,
    aes(
      x = bar_axis_2,
      y = 100,
      label = paste0(
        "min=", round(min_value, 2),
        "\nmax=", round(max_value, 2),
        "\nmean=", round(mean_value, 2)
      )
    ),
    size = 3,
    color = "black",
    lineheight = 0.9,
    vjust = 0
  )
p_range


################################################################################
write.csv(range_summary,paste0(output_plot,'Figure_',figure_number,'_ranges_',current_date,'.csv'),row.names=FALSE)
ggsave(paste0(output_plot,'Figure_',figure_number,'_ranges_',current_date,'.png'),p_range,width=12,height=12,dpi=dpi)


#ggsave(paste0(output_plot,plot_title,'_PPR_',current_date,'.png'),width=10,height=6.6,dpi=dpi)
#write.csv(df_plot,paste0(output_plot,plot_title,'_PPR_',current_date,'.csv'))
################################################################################
################################################################################
################################################################################
################################################################################
# df_sample_counts=df_plot%>%
#   filter(!is.na(sample_reads)) %>% 
#   select(sample_id,bar_color,bar_type_2,sample_reads) %>% 
#   mutate(bar_color=factor(bar_color, levels = rev(bar_color_order)))%>%
#   mutate(bar_type_2=factor(bar_type_2, levels = bar_type_2_order))%>%
#   distinct()
df_sample_counts=df_plot %>% 
select(sample_id,bar_color,bar_type_2,sample_reads) %>% 
  distinct()

df_summary = df_sample_counts %>%
  mutate(value=sample_reads) %>% 
  mutate(log10=log10(value)) %>% 
  group_by(bar_type_2, bar_color) %>%
  summarize(
    log10_mean=mean(log10),
    log10_sum = sum(log10, na.rm = TRUE),
    log10_mean = mean(log10, na.rm = TRUE),
    log10_median = median(log10, na.rm = TRUE),
    log10_min  = min(log10, na.rm = TRUE),
    log10_max  = max(log10, na.rm = TRUE),
    log10_sd   = sd(log10,  na.rm = TRUE),
    sum = sum(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    min  = min(value, na.rm = TRUE),
    max  = max(value, na.rm = TRUE),
    sd   = sd(value,  na.rm = TRUE),
    n    = sum(!is.na(value))   # or just n()
  )

sample_bar1=ggplot(df_sample_counts, aes(x = bar_color, y = sample_reads)) +
  facet_grid2(~bar_type_2,
              # ,  labeller = labeller(
              #   .rows = as_labeller(palette_label),
              #   .cols = as_labeller(palette_label)
              #)
              ,scale = "free"
              #,nrow=1
              # ,space = "free"
  )+

  stat_summary(
    aes(fill = bar_color),
    fun = mean,
    geom = "bar"
  )+
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),   # mult = 1 → 1 SD, mult = 2 → 2 SD
    geom = "errorbar",
    width = 0.2
  ) +
  #geom_col(data=df_summary,aes(x=bar_color,y=10^log10_mean,fill = bar_color),width=.3,alpha=.3)+
  scale_y_log10(labels=log_labels_bold)  +
  
  labs(x='',y='',
       title='Sequence Reads')+
  scale_color_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
  scale_fill_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
  theme_global(base_size = 11)+
  theme(       
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # axis.text.x = element_markdown(angle = 45, 
    #                            #size = axis_text_size,
    #                            hjust=1,vjust=1),
    legend.position='bottom'
  )
  # stat_count(
  #   aes(x=bar_color,
  #       
  #       label = after_stat(paste0("n = ", count))
  #   ),
  #   inherit.aes = FALSE,
  #   geom = "text",
  #   vjust = 0,
  #   fontface = "bold",
  #   size = 3.5
  # )
sample_bar1

ggsave(paste0(output_plot,'Figure-',figure_number,'_sample_bar.png'),sample_bar1)

################################################################################
################################################################################

sc=4

theme_custom=theme(
  legend.key.size = unit(.9, "cm"),
  plot.title = element_markdown(size = rel(1.6), face = 'bold',hjust=0),
  axis.text.y = element_markdown(size = rel(1.2), face = 'bold',margin = margin(r = 10)),
  strip.text = element_text(size = rel(1.35), face = 'bold'),
 legend.text = element_text(size = rel(1.35)),
 legend.title = element_markdown(size = rel(1.35), face = 'bold'),
  axis.text = element_markdown(size = rel(1.2), face = 'bold')
  
)

st_plot1=readRDS(paste0(output_plot,'Figure-',figure_number,'_st_plot.rds'))
st_plot1
mn_plot1=readRDS(paste0(output_plot,'Figure-',figure_number,'_mn_plot.rds'))
sample_plot1=readRDS(paste0(output_plot,'Figure-',figure_number,'_sample_plot.rds'))
sample_plot1


mn_plot1
mn_plot=mn_plot1+
  theme_custom+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(b = 0,l=5,r=20)
  )+
  guides(color='none',fill='none')+
  ggtitle('C) Mean Differential Abundance')
mn_plot

sample_plot=sample_plot1+
  # guides(color='none',fill='none')+
  ggtitle('D) Sequence Reads')+
  theme_custom+
  theme(    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(size = rel(1.2)),
    legend.position='right',
    legend.margin = margin(t=0, r=26.5, b=0, l=28.5),
    plot.margin=margin(
      #l=310,
      l=60)
  )+
  labs(y='')
sample_plot



st_plot1
st_plot=st_plot1+
  theme_custom+
  theme(
#legend.margin=margin(),
    plot.margin=margin(t=10,l=60)
  )+
  guides(color='none',fill='none')+
  ggtitle('B) Percent Relative Abundance')+
  scale_x_log10(breaks=c(1,10,100),labels=log10_labels_percent,limits = c(1, 100))

st_plot

ggplot_build(st_plot)$layout$panel_params[[1]]$x$breaks



sample_plot1

p=gPlot(p1)+theme_custom+
  # theme(    plot.title = element_markdown(size = rel(1.6), face = 'bold',hjust=0),
  #                    plot.caption = element_text(hjust = 0.5))
#  guides(color='none',fill='none')+
  ggtitle('A) Recovery of 16S rRNA Gene Copies and Sequence Reads from Negative Controls')+
  
  theme(legend.position ='right',
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x  = element_blank(),
        plot.margin  = margin(b = 0),
        legend.margin = margin(t=20, r=90, b=0, l=25)
        )
p

sample_bar=sample_bar1+
  theme_custom+
  labs(y='Sequence Reads')+
  theme(
    plot.title = element_blank(),
    strip.text = element_blank(),
    plot.margin = margin(t = 0),
    legend.position='none'
  )


    row1=(p/sample_bar)+
      plot_layout(heights = c(2, 1), guides = "collect")
    row1
    
   row1=wrap_elements(row1)
    
    row2=st_plot
  row2=wrap_elements(row2)
   
    
       #row3=wrap_elements(sample_plot)
    #row1=p
   # row2=free(st_plot)
  #  row3=sample_plot
      
      
    combined <- (row1 / row2) +
      plot_layout(heights = c(5,5))# &
    #  theme(
    #    legend.box.margin = margin(l=-240,b=550,r=115,t=0)
        
     # )
    
  

#combined

ggsave(paste0(output_plot,'Figure_',figure_number,'_Neg_Control_',current_date,'.pdf'),combined,width=16,height=16,dpi=dpi)





################################################################################





################################################################################
################################################################################