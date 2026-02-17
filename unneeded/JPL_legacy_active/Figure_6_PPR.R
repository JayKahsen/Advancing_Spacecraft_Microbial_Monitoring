source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
###############################################################################
# Description Start line 80ish
################################################################################
script_title='Figure_7_PPR'
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


theme_common<-theme_global(base_size = 11)+
  theme(       
    
  )


gPlot <- function(p) {
  p=p+ theme_common+
    labs(
      fill='Treatment Type',
      shape='Location'
    )+
    scale_y_continuous(labels = log10_labels_bold)+
    scale_x_discrete(labels = palette_label)+
    #   scale_fill_manual(values = palette_color) +
    scale_color_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
    scale_fill_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
    scale_shape_manual(values = palette_shape,labels=palette_label)
    # guides(color = "none")+
    #scale_x_log10(limits = c(1, sc),labels = adjust_labels,expand = expansion(mult = c(0,.3)))+ 
    #scale_x_log10()+ 
    #  scale_x_continuous(trans = "log") +
    #scale_x_continuous(trans = "log", labels = adjust_labels) +
    # scale_x_continuous(trans = 'log2', limits = c(1, sc), labels = adjust_labels)+
    #  scale_x_log10(labels = log_labels,breaks = log_breaks)+
    # labs(
    #   title = taxa_levs,
    #   x = "Mean Reads",
    #   y = "",
    #   caption = 'mean differential abunance'  )+
    
   # scale_x_discrete(labels = palette_label)+
  
   
  print(p)
  return(p)
}
################################################################################
# Setup
################################################################################
################################################################################
# Prepare Data
################################################################################
################################################################################
# Libraries
################################################################################

figure_number='7'
################################################################################
plot_title='bacterial left 1'
################################################################################
# df_plot1= df_bacterial%>%
#   pivot_longer(cols = c('No PMA',PMA),names_to = 'Treatment',values_to = 'value')
df_plot1=meta %>% 
  mutate(value_non_log10=value_bacterial) %>% 
  mutate(
    value = ifelse(value_non_log10 < 0, NA,
                   log10(ifelse(value_non_log10 == 0, 1, value_non_log10)))
  ) %>% 
  filter(!is.na(value)) %>% 
  filter(treatment %in% treatment_order) %>% 
  filter(location %in% location_order) %>% 
  ungroup()
names(df_plot1)

################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "treatment",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = 'location',
  testing_col= "treatment",
  test_type = NULL,
  adjust_method = "BH",
  # dodge_width = 0.8,
  width = 0.8,
  log = TRUE,
  space_adj = 1,
  tip_adj = 1,
  run_tests = TRUE,
)

df_summary  <- results$summary
df_segments <- results$segments
df_tests    <- results$tests
################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "treatment",
  #  dodge_column = "extraction_method",
  value_column = "value_non_log10",
  facet_columns = 'location',
  testing_col= "treatment",
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
  axis_column = "treatment",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = 'location',
  testing_col= "treatment",
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
              select(location,y)) %>% 
  mutate(label_value=y)
################################################################################
################################################################################
library(ggplot2)
library(ggsignif)

p=ggplot(df_plot1, aes(x = treatment, y = value, fill = treatment)) +
  geom_jitter(width = 0.2, size=point_size,alpha=0,fill='white'
              #shape = 21
              ) +
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
    size = 3,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3,vjust=1.5
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = treatment,
      ymin = mean_value,
      ymax = mean_value,
      color=treatment
    ),
    width = 0.6, size = 1.5, inherit.aes = FALSE
    #, color=bar_color
  )+
  geom_jitter(aes(shape=location),width = 0.2, size=point_size,
              #shape = 21
  ) +

  labs(
    fill = 'Treatment Type',
    Shape = 'Location',
 #   color = 'Treatment Type',
 title = 'A) KSC-PHSF facility surface- (m<sup>2</sup>) - Bacteria',
    y = "log10(16S rRNA gene copies) / m<sup>2</sup>",
    #x = 'Treatment Type'
    x=''

  ) +
  guides(  color='none',
           alpha='none',
           fill = guide_legend(override.aes = list(shape = 21)))+
  facet_grid(~location)

p

plotA=gPlot(p)
plotA

################################################################################
# ################################################################################
# plot_title='bacterial left2'
# ################################################################################
# 
# df_summary <- df_bacterial %>%
#   distinct(sample_name, No_PMA_minus_PMA) %>%
#   summarise(
#     mean_value = mean(No_PMA_minus_PMA),
#     sd_value = sd(No_PMA_minus_PMA)
#   )%>%
#   mutate(
#     ymin = mean_value - sd_value,
#     ymax = mean_value + sd_value
#   )
# 
# # Prepare distinct difference data
# df_diff <- df_plot1 %>% distinct(sample_name, No_PMA_minus_PMA)
# 
# #log_breaks <- c(0,10^seq(0, 12, by = 1))
# 
# # Left panel: log10-transformed values for cotton and macrofoam
# p_left <- ggplot(df_plot1, aes(x = Treatment, y = value, group = sample_name)) +
#   geom_line(color=bar_color) +
#   geom_point(aes(fill = Treatment), size=point_size, shape = 21) +
#  # coord_cartesian(ylim = zoom)+
#   scale_y_log10(
#     labels = log_labels_bold,
#     breaks = log_breaks,
#     limits = log_limits
#   )+
#   theme_bw() +
#   common_theme +
#   theme(
#     #axis.title.x = element_blank(),
#     legend.position = "none",
#     plot.margin = margin(5, 0, 5, 5)  # no space on right side
#   )+
#   labs(
#     fill='',
#     color='',
#     title='',
#     y="log10(16S rRNA gene copies) / 25 cm<sup>2</sup>",
#     x=''
#   )
# p_left=gPlot(p_left)
# # Right panel: linear difference with axis on the right
# p_right <- ggplot(df_diff, aes(x = 'No PMA - PMA', y = No_PMA_minus_PMA)) +
#   geom_point(shape = 21, fill = "gray30", size=point_size) +
#   geom_segment(data = df_summary, aes(
#     x = 1 - 0.2,
#     xend = 1 + 0.2,
#     y = mean_value,
#     yend = mean_value
#   ),size = 1.5, inherit.aes = FALSE,color=bar_color)+
#   
#   geom_errorbar(data = df_summary, aes(x = 'No PMA - PMA', ymin = ymin, ymax = ymax), 
#                 width = 0.1, size = 1, inherit.aes = FALSE,color=bar_color)+
#   # scale_y_continuous(name = "Mean of differences", position = "right") +
#   scale_y_continuous( position = "right") +
#   theme_bw() +
#   common_theme +
#   theme(
#     # axis.title.x = element_blank(),
#     plot.margin = margin(5, 5, 5, 0),  # no space on left side
#     legend.position = "none"
#   )+
#   labs(
#     fill='',
#     color='',
#     title='',
#     y="Mean of differences",
#     x=''
#   )
# p_right
# p_right=gPlot(p_right)
# 
# # Combine the two plots side by side with touching panels
# combined_plot <- p_left + p_right +
#   plot_layout(
#     widths = c(2.2, 1),
#     ##  guides = "collect"
#   ) &
#   theme(
#     #  legend.position = "bottom",
#     plot.margin = margin(0, 0, 0, 0),  # remove external spacing
#     legend.position = "none")
# 
# # Print the final combined plot
# plotB=combined_plot
# plotB
# 

################################################################################
plot_title='fungal left 1'
################################################################################

################################################################################
# df_plot1= df_bacterial%>%
#   pivot_longer(cols = c('No PMA',PMA),names_to = 'Treatment',values_to = 'value')
df_plot2=meta %>% 
  mutate(value_non_log10=value_fungal) %>% 
  mutate(
    value = ifelse(value_non_log10 < 0, NA,
                   log10(ifelse(value_non_log10 == 0, 1, value_non_log10)))
  ) %>% 
  filter(!is.na(value)) %>% 
  filter(treatment %in% treatment_order) %>% 
  filter(location %in% location_order) %>% 
  ungroup()

names(df_plot2)

################################################################################
results <- summary_brackets(
  df = df_plot2,
  axis_column = "treatment",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = 'location',
  testing_col= "treatment",
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
results <- summary_brackets(
  df = df_plot2,
  axis_column = "treatment",
  #  dodge_column = "extraction_method",
  value_column = "value_non_log10",
  facet_columns = 'location',
  testing_col= "treatment",
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
  axis_column = "treatment",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = 'location',
  testing_col= "treatment",
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
              select(location,y)) %>% 
  mutate(label_value=y)
################################################################################
################################################################################
library(ggplot2)
library(ggsignif)

p=ggplot(df_plot2, aes(x = treatment, y = value, fill = treatment)) +
  geom_jitter(width = 0.2, size=point_size,alpha=0,fill='white'
              #shape = 21
  ) +
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
    size = 3,vjust=-.5
  )+
  geom_text(
    data = df_tests,
    aes(
      x = label_position,
      y = label_value ,
      label = scales::scientific(p.value)
    ),
    inherit.aes = FALSE,
    size = 3,vjust=1.5
  )+
  # geom_errorbar(
  #   data = df_summary,
  #   aes(
  #     x = treatment,
  #     ymin = mean_minus_sd,
  #     ymax = mean_plus_sd,
  #   ),
  #   width = 0.1, size = 1, inherit.aes = FALSE, color=bar_color
  # )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = treatment,
      ymin = mean_value,
      ymax = mean_value,
      color=treatment
    ),
    width = 0.4, size = 1.5, inherit.aes = FALSE
    #, color=bar_color
  )+
  geom_jitter(aes(shape=location),width = 0.2, size=point_size
              #shape = 21
  ) +
  
  
  ## coord_cartesian(ylim = zoom)+
  scale_y_log10(
    labels = log_labels_bold,
    breaks = log_breaks,
    limits = log_limits
  )+
  
  labs(
    shape='Location',
    fill = 'Treatment Type',
    #   color = 'Treatment Type',
    title = 'B) KSC-PHSF facility surface (m<sup>2</sup>) - Fungi',
    y = "log10(16S rRNA gene copies) / m<sup>2</sup>",
    #x = 'Treatment Type'
    x=''
    
  ) +
  guides(
  color='none',
  alpha='none',
  fill = guide_legend(override.aes = list(shape = 21)))+
  facet_grid(~location)

p

plotB=gPlot(p)
plotB
################################################################################

# final_plot <- ((plotA |plotB)) +
#   plot_annotation(
#     tag_levels = 'A',
#     tag_suffix = '.',
#     theme = theme(
#       plot.tag = element_text(size = 16, face = "bold")
#     )
#   )
# 
# library(patchwork)
# 
# final_plot <- (
#   (plotA / wrap_elements(plotB)) |
#     (plotC / wrap_elements(plotD))
# ) +
#   plot_annotation(
#     tag_levels = 'A',
#     tag_suffix = '.',
#     theme = theme(
#       plot.tag = element_text(size = 16, face = "bold")
#     )
#   )
# library(patchwork)
# library(ggplot2)
# 
# plot_title=paste(script_title)
# 
# # Add right margin to A and C
# plotA_adj <- plotA + theme(plot.margin = margin(r=55))  # top, right, bottom, left
# plotC_adj <- plotC + theme(plot.margin = margin(r=55))
# 
# # Wrap composite plots B and D
# plotB_wrap <- wrap_elements(plotB)
# plotD_wrap <- wrap_elements(plotD)

# Final combined layout

plotB=plotB+theme(legend.position='none')

plot_title=script_title
final_plot_right <- (plotA | plotB)+
plot_layout(guides = "collect")

final_plot_right

ggsave(paste0(output_plot,plot_title,current_date,'_right.pdf'),final_plot_right,width=19,height=8,dpi=dpi)


final_plot <- (plotA | plotB) +
  plot_annotation(
    caption = "legend in the middle as it pertains to both? or to move it back to the right\nBacteria, and Fungi, or Bacterial and Fungal, figure 5 is the latter",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = 14)
    )
  )

#   (wrap_elements(plotA_adj) / plotB_wrap) |
#     (wrap_elements(plotC_adj) / plotD_wrap)
# ) +
  # +
#  plot_annotation(
#   title = '',
  #  tag_levels = 'A',
#  #  tag_suffix = '.',
#    theme = theme(
 #     plot.tag = element_text(size = 16, face = "bold")
 #   )
#  )

final_plot

ggsave(paste0(output_plot,plot_title,current_date,'_middle_question.pdf'),final_plot,width=19,height=8,dpi=dpi)

plotD=plotB+labs(color='',shape='',fill='')+
  stat_count(aes(x = treatment,
                 label = after_stat(paste0("n = ", count))
  ),
  inherit.aes = FALSE,
  geom = "text",
  vjust = 0,
  fontface = "bold",
  size = 3.5
  )
plotD
plotC=plotA+labs(color='',shape='',fill='')+
  stat_count(aes(x = treatment,
                 label = after_stat(paste0("n = ", count))
  ),
  inherit.aes = FALSE,
  geom = "text",
  vjust = 0,
  fontface = "bold",
  size = 3.5
  )
plotC

final_plot_bottom <- (plotC | plotD) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  ) 

final_plot_bottom

ggsave(paste0(output_plot,plot_title,current_date,'.pdf'),final_plot_bottom,width=19,height=8,dpi=dpi)





write.csv(df_plot1,paste0(output_plot,plot_title,current_date,'_bacterial.csv'),row.names = FALSE)
write.csv(df_plot2,paste0(output_plot,plot_title,current_date,'_fungal.csv'),row.names = FALSE)
################################################################################
