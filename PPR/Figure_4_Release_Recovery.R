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
names(meta)

figure_number='4'

experiment_type_order=c( 'metal_deposition','swab_head_retention' )

df=meta %>% 
  filter(experiment_type%in%c( 'metal_deposition','swab_head_retention')) %>% 
  filter(!is.na(percent_recovery_dpcr)) %>% 
  filter(sampling_device %in% sampling_device_order) %>% 
  select(sample_id,percent_recovery_qpcr,percent_recovery_dpcr,sampling_device,experiment_type,figure,Figure,bar_fill,sample_reads)

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

df_summary=df_summary1 %>% filter(figure %in% c('4A','4B','4C'))%>% 
  mutate(method=factor(method,levels=method_order))
df_segments=df_segments1%>% filter(figure %in% c('4A','4B','4C'))%>% 
  mutate(method=factor(method,levels=method_order))
df_tests=df_tests1  %>% filter(figure %in% c('4A','4B','4C')) %>% 
  mutate(method=factor(method,levels=method_order))

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
  guides(color = "none")+
  scale_x_discrete(labels = palette_label)

p


plot_A_13=gPlot(p)+theme(legend.position='right')
plot_A_13
################################################################################
message('reads_plot')
################################################################################
df_plot2 <- df %>%
  filter(str_detect(Figure, figure_number)) %>% 
  mutate(sampling_device=factor(sampling_device,levels=sampling_device_order))

################################################################################
#──────────────────────────────────────────────
# Function: plot_reads_by_device
# Example: plot_reads_by_device(df, "4C")
#──────────────────────────────────────────────
plot_reads_by_figure <- function(df, figure_id) {
  require(ggplot2)
  require(dplyr)
  require(ggh4x)
  require(stringr)
  
  df_plot_data <- df %>%
    filter(str_detect(Figure, figure_id)) %>%
    mutate(
      sampling_device = factor(
        sampling_device,
        levels = sampling_device_order
      )
    )
  
  p <- ggplot(df_plot_data, aes(x = sampling_device, y = sample_reads, fill = sampling_device)) +
    facet_grid2(~figure, scales = "free_x", axes = "all", remove_labels = "none") +
    scale_y_log10(labels = log_labels_bold) +
    stat_summary(fun = mean, geom = "col", color = "black", width = 0.5) +
    geom_jitter(width = 0.1, alpha = 0.7, size = 2,shape=21)+
    labs(
      y = "Mean Sequence Reads",
      x = "",
      color = "",
      fill = ""
    )
  p
  # +
  #   coord_cartesian(ylim = c(10^5.5, 10^7.5))+
  #   coord_cartesian(clip = "off")

  
  p_out <- gPlot(p) +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      strip.text = element_blank()
    ) +

     coord_cartesian(ylim = c(10^4.5, 10^7.65))+
    ggtitle('')
  #   coord_cartesian(clip = "off")
  
  return(p_out)
}

################################################################################
df_plot2$sampling_device
reads_plot_4A=plot_reads_by_figure(df=df_plot2,'4A')
reads_plot_4A

reads_plot_4B=plot_reads_by_figure(df=df_plot2,'4B')
reads_plot_4B

reads_plot_4C=plot_reads_by_figure(df=df_plot2,'4C')
reads_plot_4C

################################################################################

################################################################################






  ################################################################################
sc=3
st_plot=readRDS(paste0(output_plot,'Figure_',figure_number,'_st_plot.rds'))+
  
  scale_x_log10_shift(
    shift = 2,                                   # shift right by 10^2 = ×100
    breaks = c(0.01, 0.1, 1, 10, 100),           # logical/original values
    labels = log10_labels_percent,
    limits = c(0.01, 100)
  )

  
  
  
  
 # scale_x_log10(breaks=c(1,10,100),labels=log10_labels_percent,limits = c(1, 100))

st_plot


theme_custom=theme(
  axis.text.y = element_markdown(size = rel(1.2), face = 'bold'),
  legend.text = element_text(size =rel(1.5)),
  legend.title = element_markdown(size =rel(1.5)),
  strip.text = element_markdown(size = rel(1.4), face = 'bold'),
  axis.title.x = element_markdown(size = rel(1.2), face = 'bold',color='black'),
  axis.text = element_markdown(size = rel(1.2), face = 'bold')
  
)


################################################################################

plot_height=12
plot_width=20


plot_A_13=plot_A_13+
  theme_custom+
  # stat_count(aes(x = sampling_device,
  #                label = after_stat(paste0("n = ", count))
  # ),
  # inherit.aes = FALSE,
  # geom = "text",
  # vjust = 0,
  # fontface = "bold",
  # size = 3.5
  # )+
  theme(legend.box.margin = margin(l=-20))+
  ggtitle('')
plot_A_13



st_plot=st_plot+
  theme_custom+
  theme(
#  axis.title.y = element_text(margin = margin(r = 115)),
   # axis.text.y = element_markdown(size = rel(1.2), face = 'bold',margin = margin(l = 0)),
    legend.position = 'none',
        plot.margin=margin(l=54.5,r=43.5)
        )
st_plot

row1=plot_A_13

row1=wrap_elements(row1)
row1

row2=st_plot
row2

row2=wrap_elements(row2)

combined1=row1/row2+
  plot_layout(heights = c(10, 12))
  
combined1
top1=.55
bottom1=.088
width1=.08
width2=.274

left2_adj=.00

left1=.208; right1=left1+width1
left2=left1+width2-left2_adj; right2=left2+width1
left3=left2+width2+left2_adj; right3=left3+width1



combined <- combined1 + 
  inset_element(reads_plot_4A,left = left1, bottom = bottom1,right = left1+width1, top = top1)+
  inset_element(reads_plot_4B,left = left2, bottom = bottom1,right = left2+width1, top = top1)+
  inset_element(reads_plot_4C,left = left3, bottom = bottom1,right = left3+width1, top = top1)
combined

ggsave(paste0(output_plot,'Figure_',figure_number,'_Release_Recovery_',current_date,'.png'),combined,width=plot_width,height=plot_height,dpi=dpi)

ggsave(paste0(output_plot,'Figure_',figure_number,'_Release_Recovery_BH_corrected_',current_date,'.pdf'),combined,width=plot_width,height=plot_height,dpi=dpi)

######################################################################################################
################################################################################
# RANGE SUMMARY FOR SECOND DATASET
################################################################################

range_summary2 <- df_plot1 %>%
  group_by(method, figure, sampling_device) %>%
  summarise(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

range_summary2
write.csv(range_summary2,paste0(output_plot,'fig_4_ranges.csv'),row.names = FALSE)

################################################################################
# ADD RANGE LABELS TO SECOND PLOT
################################################################################

p_range2 <- p +
  geom_text(
    data = range_summary2,
    aes(
      x = sampling_device,
      y = 50,   # or: max(df_plot1$value, na.rm=TRUE) * 1.05
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

p_range2


ggsave(paste0(output_plot,'Figure_',figure_number,'_Release_Recovery_ranges.png'),p_range2,width=12,height=6)


