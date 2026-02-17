source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
###############################################################################
# Description Start line 80ish
################################################################################
script_title='ILMAH'
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


# common_theme <- theme(
#   
#   # legend.justification = "center",
#   # legend.margin = margin(t = -10, r = 10, b = 10, l = 10),
#   #legend.box.margin = margin(,b=20), 
#   panel.border = element_rect(color = "black", fill = NA, size = 1),
#   legend.background = element_rect(color = "black"),
#   #plot.background = element_rect(color = "black"),
#   legend.text = element_text(face = "bold", size = legend_text_size), 
#   axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",color='black',size=axis_text_size),
#   axis.text.y = element_text(face = "bold",color='black',size=axis_text_size),
#   axis.ticks = element_blank(),
#   plot.title = element_blank(),
#   axis.title = element_text(face = "bold", size = axis_title_size), 
#   plot.caption = element_text(hjust = 0.5)
# )
# 
# plot_width=10
# plot_height=10
# margin_buffer=30


#You can start X-axis with Timepoint 4 at 0.
#Can you please use the scale limits for the Y-axis to be from 0 – 4.


# gPlot <- function(p) {
#   p=p+
#     scale_color_manual(values = palette_color,labels=palette_label)+
#     scale_fill_manual(values = palette_color,labels=palette_label)+
#     theme_bw() +
#     common_theme
#     # labs(x = "", y = "Mean Shannon", color = "", fill = "",title='',
#     #      caption='assuming normal distribution instead of individual t distributions')+
#     
#   #  scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) 
#   print(p)
#   return(p)
# }

################################################################################
# load files
################################################################################
file_path=paste0('data_tables/',script_title,'/')
library(stringr)

df_bacterial <- read_excel("data_tables/RUSH_woof_tidy FIGURE 6 CLIPPER.xlsx", sheet = 2)%>%
  rename(value=`Bacterial Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`)%>%
  select(`Subject ID`,Treatment,Location,value)%>%
  rename(sample_name=`Subject ID`)%>%
  mutate(value=as.numeric(value))%>%
  filter(!is.na(value))%>%
  mutate(sample_name = sub("-[^-]*$", "", sample_name))%>%
  pivot_wider(names_from = Treatment, values_from = value) %>%
  mutate(No_PMA_minus_PMA = `No PMA` - PMA)%>%

ungroup()
  

df_fungal <- read_excel("data_tables/RUSH_woof_tidy FIGURE 6 CLIPPER.xlsx", sheet = 3)%>%
  rename(value=`Fungall Digital PCR conversion (1uL to 1m2) dPCR*50*[Column L/0.375]`)%>%
  select(`Subject ID`,Treatment,Location,value)%>%
  rename(sample_name=`Subject ID`)%>%
  mutate(value=as.numeric(value))%>%
  filter(!is.na(value))%>%
  mutate(sample_name = sub("-[^-]*$", "", sample_name))%>%
  pivot_wider(names_from = Treatment, values_from = value) %>%
  mutate(No_PMA_minus_PMA = `No PMA` - PMA)%>%
  
  ungroup()


# %>%
#   rename(total_copies='Total Copies in Orginal Sample (50uL) [dPCR copies*dil*50]')%>%
#   select(unique_replicate_id,Treatment,total_copies)%>%
#   mutate(sample_name = str_remove_all(unique_replicate_id, "[mspc]"))%>%
#   select(-unique_replicate_id)%>%
#   pivot_wider(names_from = Treatment,values_from = total_copies)%>%
#   mutate(
#     No_PMA_minus_PMA = No PMA - PMA,
#     polyester_wipe_minus_salsa = `Polyester Wipe`-SALSA 
#   )%>%
#   ungroup()
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
plot_title='bacterial left 1'
################################################################################
df_plot1= df_bacterial%>%
  pivot_longer(cols = c('No PMA',PMA),names_to = 'Treatment',values_to = 'value')
  

################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "Treatment",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = NA,
  testing_col= "Treatment",
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
library(ggplot2)
library(ggsignif)

p=ggplot(df_plot1, aes(x = Treatment, y = value, fill = Treatment)) +
  geom_jitter(width = 0.2, size=point_size, alpha = 0.8, shape = 21) +
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
  geom_errorbar(
    data = df_summary,
    aes(
      x = Treatment,
      ymin = mean_minus_sd,
      ymax = mean_plus_sd,
    ),
    width = 0.1, size = 1, inherit.aes = FALSE, color=bar_color
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = Treatment,
      ymin = mean_value,
      ymax = mean_value,
    ),
    width = 0.4, size = 1.5, inherit.aes = FALSE, color=bar_color
  )+
  

 ## coord_cartesian(ylim = zoom)+
  scale_y_log10(
    labels = log_labels_bold,
    breaks = log_breaks,
    limits = log_limits
  )+
  
  labs(
    fill = 'Treatment Type',
    color = 'Treatment Type',
    title = 'JPL - Spacecraft Assembly Facilty 25 (cm<sup>2</sup>) - Bacterial',
    y = "log10(16S rRNA gene copies) / 25 cm<sup>2</sup>",
    #x = 'Treatment Type'
    x=''

  ) 

p

plotA=gPlot(p)
plotA

################################################################################
################################################################################
plot_title='bacterial left2'
################################################################################

df_summary <- df_bacterial %>%
  distinct(sample_name, No_PMA_minus_PMA) %>%
  summarise(
    mean_value = mean(No_PMA_minus_PMA),
    sd_value = sd(No_PMA_minus_PMA)
  )%>%
  mutate(
    ymin = mean_value - sd_value,
    ymax = mean_value + sd_value
  )

# Prepare distinct difference data
df_diff <- df_plot1 %>% distinct(sample_name, No_PMA_minus_PMA)

#log_breaks <- c(0,10^seq(0, 12, by = 1))

# Left panel: log10-transformed values for cotton and macrofoam
p_left <- ggplot(df_plot1, aes(x = Treatment, y = value, group = sample_name)) +
  geom_line(color=bar_color) +
  geom_point(aes(fill = Treatment), size=point_size, shape = 21) +
 # coord_cartesian(ylim = zoom)+
  scale_y_log10(
    labels = log_labels_bold,
    breaks = log_breaks,
    limits = log_limits
  )+
  theme_bw() +
  common_theme +
  theme(
    #axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 0, 5, 5)  # no space on right side
  )+
  labs(
    fill='',
    color='',
    title='',
    y="log10(16S rRNA gene copies) / 25 cm<sup>2</sup>",
    x=''
  )
p_left=gPlot(p_left)
# Right panel: linear difference with axis on the right
p_right <- ggplot(df_diff, aes(x = 'No PMA - PMA', y = No_PMA_minus_PMA)) +
  geom_point(shape = 21, fill = "gray30", size=point_size) +
  geom_segment(data = df_summary, aes(
    x = 1 - 0.2,
    xend = 1 + 0.2,
    y = mean_value,
    yend = mean_value
  ),size = 1.5, inherit.aes = FALSE,color=bar_color)+
  
  geom_errorbar(data = df_summary, aes(x = 'No PMA - PMA', ymin = ymin, ymax = ymax), 
                width = 0.1, size = 1, inherit.aes = FALSE,color=bar_color)+
  # scale_y_continuous(name = "Mean of differences", position = "right") +
  scale_y_continuous( position = "right") +
  theme_bw() +
  common_theme +
  theme(
    # axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 5, 0),  # no space on left side
    legend.position = "none"
  )+
  labs(
    fill='',
    color='',
    title='',
    y="Mean of differences",
    x=''
  )
p_right
p_right=gPlot(p_right)

# Combine the two plots side by side with touching panels
combined_plot <- p_left + p_right +
  plot_layout(
    widths = c(2.2, 1),
    ##  guides = "collect"
  ) &
  theme(
    #  legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0),  # remove external spacing
    legend.position = "none")

# Print the final combined plot
plotB=combined_plot
plotB


################################################################################
plot_title='fungal left 1'
################################################################################
df_plot1= df_fungal%>%
  pivot_longer(cols = c('No PMA',PMA),names_to = 'Treatment',values_to = 'value')


################################################################################
results <- summary_brackets(
  df = df_plot1,
  axis_column = "Treatment",
  #  dodge_column = "extraction_method",
  value_column = "value",
  facet_columns = NA,
  testing_col= "Treatment",
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
library(ggplot2)
library(ggsignif)

p=ggplot(df_plot1, aes(x = Treatment, y = value, fill = Treatment)) +
  geom_jitter(width = 0.2, size=point_size, alpha = 0.8, shape = 21) +
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
  geom_errorbar(
    data = df_summary,
    aes(
      x = Treatment,
      ymin = mean_minus_sd,
      ymax = mean_plus_sd,
    ),
    width = 0.1, size = 1, inherit.aes = FALSE, color=bar_color
  )+
  geom_errorbar(  # Mean segment
    data = df_summary,
    aes(
      x = Treatment,
      ymin = mean_value,
      ymax = mean_value,
    ),
    width = 0.4, size = 1.5, inherit.aes = FALSE, color=bar_color
  )+
  
  
 # coord_cartesian(ylim = zoom)+
  scale_y_log10(
    labels = log_labels_bold,
    breaks = log_breaks,
    limits = log_limits
  )+
  
  labs(
    fill = 'Treatment Type',
    color = 'Treatment Type',
    title = 'JPL - Spacecraft Assembly Facilty 1000 (cm<sup>2</sup>) - Fungal',
    y = "log10(16S rRNA gene copies) / 1000 cm<sup>2</sup>",
    #x = 'Treatment Type'
    x = ''
  ) 

p

plotC=gPlot(p)
plotC

################################################################################
################################################################################
plot_title='fungal left2'
################################################################################

df_summary <- df_fungal %>%
  distinct(sample_name, No_PMA_minus_PMA) %>%
  summarise(
    mean_value = mean(No_PMA_minus_PMA),
    sd_value = sd(No_PMA_minus_PMA)
  )%>%
  mutate(
    ymin = mean_value - sd_value,
    ymax = mean_value + sd_value
  )

# Prepare distinct difference data
df_diff <- df_plot1 %>% distinct(sample_name, No_PMA_minus_PMA)

#log_breaks <- c(0,10^seq(0, 12, by = 1))

# Left panel: log10-transformed values for cotton and macrofoam
p_left <- ggplot(df_plot1, aes(x = Treatment, y = value, group = sample_name)) +
  geom_line(color=bar_color) +
  geom_point(aes(fill = Treatment), size=point_size, shape = 21) +
 # coord_cartesian(ylim = zoom)+
  scale_y_log10(
    labels = log_labels_bold,
    breaks = log_breaks,
    limits = log_limits
  )+
  theme_bw() +
  common_theme +
  theme(
    #axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 0, 5, 5)  # no space on right side
  )+
  labs(
    fill='',
    color='',
    title='',
    y="log10(16S rRNA gene copies) / 1000 cm<sup>2</sup>",
    x=''
  )
p_left=gPlot(p_left)
# Right panel: linear difference with axis on the right
p_right <- ggplot(df_diff, aes(x = 'No PMA - PMA', y = No_PMA_minus_PMA)) +
  geom_point(shape = 21, fill = "gray30", size=point_size) +
  geom_segment(data = df_summary, aes(
    x = 1 - 0.2,
    xend = 1 + 0.2,
    y = mean_value,
    yend = mean_value
  ),size = 1.5, inherit.aes = FALSE,color=bar_color)+
  
  geom_errorbar(data = df_summary, aes(x = 'No PMA - PMA', ymin = ymin, ymax = ymax), 
                width = 0.1, size = 1, inherit.aes = FALSE,color=bar_color)+
  # scale_y_continuous(name = "Mean of differences", position = "right") +
  scale_y_continuous( position = "right") +
  theme_bw() +
  common_theme +
  theme(
    # axis.title.x = element_blank(),
    plot.margin = margin(5, 5, 5, 0),  # no space on left side
    legend.position = "none"
  )+
  labs(
    fill='',
    color='',
    title='',
    y="Mean of differences",
    x=''
  )
p_right
p_right=gPlot(p_right)

# Combine the two plots side by side with touching panels
combined_plot <- p_left + p_right +
  plot_layout(
    widths = c(2.2, 1),
    ##  guides = "collect"
  ) &
  theme(
    #  legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0),  # remove external spacing
    legend.position = "none" )

# Print the final combined plot
plotD=combined_plot
plotD

################################################################################

final_plot <- ((plotA / plotB) | (plotC / plotD)) +
  plot_annotation(
    tag_levels = 'A',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )

library(patchwork)

final_plot <- (
  (plotA / wrap_elements(plotB)) |
    (plotC / wrap_elements(plotD))
) +
  plot_annotation(
    tag_levels = 'A',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )
library(patchwork)
library(ggplot2)

# Add right margin to A and C
plotA_adj <- plotA + theme(plot.margin = margin(r=55))  # top, right, bottom, left
plotC_adj <- plotC + theme(plot.margin = margin(r=55))

# Wrap composite plots B and D
plotB_wrap <- wrap_elements(plotB)
plotD_wrap <- wrap_elements(plotD)

# Final combined layout
final_plot <- (
  (wrap_elements(plotA_adj) / plotB_wrap) |
    (wrap_elements(plotC_adj) / plotD_wrap)
) +
  plot_annotation(
    tag_levels = 'A',
    tag_suffix = '.',
    theme = theme(
      plot.tag = element_text(size = 16, face = "bold")
    )
  )





final_plot

ggsave(paste0(output_plot,'figure6 PPR.png'),final_plot,width=18,height=12)
################################################################################
