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
log_limits=c(NA, NA)
zoom=c(1e3, 1e11)
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


# if (!require(remotes)) {
#   install.packages("remotes")
# }
# remotes::install_github('njudd/ggrain')

#library(gghalves)

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

plot_width=10
plot_height=10
margin_buffer=30


#You can start X-axis with Timepoint 4 at 0.
#Can you please use the scale limits for the Y-axis to be from 0 – 4.


gPlot <- function(p) {
  p=p+
    scale_color_manual(values = color_palette)+
    scale_pattern_manual(values = pattern_palette)+
    scale_linetype_manual(values = linetype_palette)+
    labs(
      fill='',
      color='',
      pattern='',
      shape='',
      linetype=''
    )
    # labs(x = "", y = "Mean Shannon", color = "", fill = "",title='',
    #      caption='assuming normal distribution instead of individual t distributions')+
    
   # scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) 
  print(p)
}

################################################################################
# load files
################################################################################
file_path=paste0('data_tables/',script_title,'/')

df <- read.csv("data_tables/data_May22.csv")

pasteC(unique(df$Location))

location=location_order=c('Entrance','Kitchen','Bathroom','Crew Quarters','Exercise Module',
                          'Medical Bay','Plant Production','EVA Module')
metric=metric_order=c('ATP_RLU','Bacterial_CFU','Fungal_CFU','Bacterial_dPCR','Fungal_dPCR')
#timepoint_order=c('Day 0','Day 5','Day 10','Day 10B','Day 15','Day 20')



df_plot=df%>%
  filter(Location !='')%>%
  mutate(across(c(Bacterial_CFU, Fungal_CFU, ATP_RLU),
                ~ ifelse(. == "BDL", "0.0.E+00", .)))%>%
  mutate(across(c(Bacterial_CFU, Fungal_CFU, ATP_RLU),
                ~ as.numeric(gsub("\\.E", "E", .))))%>%
  
  pivot_longer(cols = any_of(metric),names_to = 'metric',values_to = 'values')%>%
  rename(timepoint=Sampling.time)%>%
 # mutate(timepoint=ifelse(timepoint %in% c('Day 10B','Day 10A'),'Day 10',timepoint))%>%



  # filter(!Location %in% c('Field Control','Neg Control'))%>%

  mutate(domain = case_when(
    str_starts(metric, "Bacterial") ~ "Bacterial",
    str_starts(metric, "Fungal") ~ "Fungal",
    str_starts(metric, "ATP") ~ "Both",
    TRUE ~ "Other"  # Default case if no condition matches
  ))%>%

 # mutate(Location = ifelse(Location == 'EVA Module', 'EVA Module', Location))%>%

  #mutate(timepoint=ifelse(timepoint=='Day 10' & Location=='EVA Module','Day 10B',timepoint))%>%
  mutate(plot_group=metric)%>%
  mutate(plot_group=factor(plot_group,levels=metric_order))%>%
  #mutate(timepoint=factor(timepoint,levels=timepoint_order))%>%
  
  mutate(Location=ifelse(Location=='Eva Module','EVA Module',Location))%>%
  mutate(metric=factor(metric,levels=metric_order))%>%
  mutate(Location=factor(Location,levels=location_order))%>%
  filter(Location %in% location_order)%>%
  ungroup()

unique(df_plot$timepoint)
unique(df$Location)



################################################################################
################################################################################
df_box=df_plot%>%
  mutate(timepoint='Day 10')

df_ATP =df_plot%>%
  select(-values)%>%
  left_join(df_plot%>%filter(metric=='ATP_RLU')%>%select(timepoint,Location,values), by = c("timepoint", "Location"))%>%
  filter(metric!='ATP_RLU')%>%
  #mutate(timepoint=factor(timepoint,levels=timepoint_order))%>%
  mutate(metric=factor(metric,levels=metric_order))%>%
  mutate(plot_group='ATP_RLU')%>%
  mutate(plot_group=factor(plot_group,levels=metric_order))

df_ATP_Fungal=df_plot%>%
  filter(metric=='ATP_RLU')%>%
  mutate(domain='Fungal')

df_ATP_Bacterial=df_plot%>%
  filter(metric=='ATP_RLU')%>%
  mutate(domain='Bacterial')

df_ATP_domain=rbind(df_ATP_Bacterial,df_ATP_Fungal)%>%
 # mutate(timepoint=factor(timepoint,levels=timepoint_order))%>%
  mutate(metric=factor(metric,levels=metric_order))%>%
  mutate(plot_group='ATP_RLU')%>%
  mutate(plot_group=factor(plot_group,levels=metric_order))

unique(df_plot$metric)


################################################################################

#plot_title='Microbial Populations of Various ILMAH Surfaces'
plot_title='Microbial Load Over Time'
plot_list=list()
df_plot2=df_plot
for (loco in location_order){
  df_plot=df_plot2%>%
    filter(Location==loco)
  
  palette_pattern_colour=c(
    'Fungal_CFU'='gray80',
    'Bacterial_CFU'='black',
    'ATP_RLU'='gray80',
    'Fungal_dPCR'='gray80',
    'Bacterial_dPCR'='gray80',
    'default'='gray80'
  )
  palette_pattern_size=c(
    'Fungal_CFU'=1,
    'Bacterial_CFU'=1,
    'ATP_RLU'=1,
    'Fungal_dPCR'=1,
    'Bacterial_dPCR'=1,
    'default'=1
  )

  df_plot$pattern_colour <- palette_pattern_colour[match(df_plot$plot_group, names(palette_pattern_colour))]
  df_plot$pattern_size <- palette_pattern_size[match(df_plot$plot_group, names(palette_pattern_size))]
  
  column_width=3
    dodge_width=column_width+.3
  ################################################################################

    
    ################################################################################
 p= ggplot(df_plot,aes(x=timepoint,y=values,group=metric,color=metric,fill=metric,pattern=metric))+
      geom_point(aes(shape=metric),alpha=0)+
 
   geom_col(data=df_plot%>%filter(metric %in% c('Bacterial_CFU','Fungal_CFU')),position = position_dodge(width = dodge_width),width=column_width)+
   geom_col(data=df_plot%>%filter(metric %in% c('Bacterial_CFU','Fungal_CFU')),position = position_dodge(width = dodge_width),color='gray30',width=column_width)+

         geom_col_pattern(data=df_plot%>%filter(metric %in% c('Bacterial_CFU','Fungal_CFU')),
                    position = position_dodge(width = dodge_width),width=column_width,
                    # aes(pattern_size  = pattern_size ,
                    #     pattern_colour = pattern_colour),
                    pattern_key_scale_factor= .2,
                 #   pattern_shape=2,
                    # pattern_spacing=.04,#.05
                    #
                    # pattern_density= 0.4,#.2
                    #   pattern_fill = "black",
                    #   pattern_density = 0.5,
                    #   pattern_spacing = 0.04,
                       color = "black"
                    )+

     #  guides( pattern_colour = "none", pattern_size = "none")+
      geom_line(data=df_plot%>%filter(metric %in% c('ATP_RLU')),aes(color = metric,linetype=metric),linewidth=1 ) +
      geom_point(data=df_plot%>%filter(metric %in% c('ATP_RLU')),aes(shape=metric,size=metric))+
      
    geom_line(data=df_plot%>%filter(metric %in% c('Bacterial_dPCR','Fungal_dPCR')),aes(linetype=metric),position = position_dodge(width = dodge_width),linewidth=1)+
    geom_point(data=df_plot%>%filter(metric %in% c('Bacterial_dPCR','Fungal_dPCR')),aes(shape=metric,size=metric),position = position_dodge(width = dodge_width),fill='white')+
   
  facet_wrap(~Location,nrow=2)+
   # scale_y_log10(labels = log_labels_bold,breaks = log_breaks)+
      coord_cartesian(ylim = zoom)+
      scale_y_log10(
        labels = log_labels_bold,
        breaks = log_breaks,
        limits = log_limits
      )+
    scale_color_manual(values = palette_color,labels=palette_label) + 
    scale_fill_manual(values = palette_color,labels=palette_label) + 
      scale_shape_manual(values = palette_shape,labels=palette_label) + 
   scale_pattern_manual(values = palette_pattern,labels=palette_label) + 
      scale_linetype_manual(values = palette_linetype,labels=palette_label) + 
      scale_size_manual(values = palette_size,labels=palette_label) + 
    labs(
      fill='',
      color='',
      pattern='',
      shape='',
      linetype='',
      size='',
         # y="Microbial Load (log scale) dPCR per m<sup>2</sup>",
         # x='Sampling Day'
         y="",
         x=''
         
    )+
    #ggtitle(paste('Microbial Load over Time - ',loco))+
    theme_bw()+
    common_theme
 p

  plot_list[[loco]]=p
}

library(patchwork)
plot_layout <- (
  wrap_plots(plot_list, nrow = 2, guides = "collect") &
    theme(legend.position = "bottom")
) +
  plot_annotation(
    title = plot_title,
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,        # Center the title
        size = 20,          # Increase font size
        face = "bold"       # Make it bold
      )
    )
  )
plot_layout
################################################################################
################################################################################
library(patchwork)
plot_title='Temporal Dynamics of Microbial Load Across Spacecraft Modules Using dPCR, CFU, and ATP Assays'
plot_title=''
y_axis_title='Microbial Load (log scale)'
x_axis_title='Sampling Days'
plot_layout <- (
  wrap_plots(plot_list, nrow = 2, guides = "collect") &
    theme(legend.position = "bottom",
          #legend.margin = margin(t = 25),
          legend.box.spacing = unit(2, "lines"))
)

final_plot <- plot_layout +
  plot_annotation(
    title = plot_title,
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold")
    )
  ) &
  theme(
    plot.margin = margin(20, 20, 40, 60)  # give room for labels
  )

# Add overall x and y axis labels using patchwork::inset_element
library(grid)
library(patchwork)
x_width=1;y_width=1;y_adj=-4.2;x_adj=-.18;y_center=.655;x_center=-2.05

final_plot_with_labels <- final_plot +
  inset_element(grid::textGrob(x_axis_title, gp = gpar(fontsize = 18)), 
                left = x_center, bottom = x_adj, right = x_center+x_width, top = x_adj+.1) +
 inset_element(grid::textGrob(y_axis_title, rot = 90, gp = gpar(fontsize = 18,fontface='bold')), 
               left = y_adj, bottom = y_center, right = y_adj+.1, top = y_center+y_width)

print(final_plot_with_labels)


################################################################################
plot_title='Microbial Load'
################################################################################
  ggsave(paste0(output_plot,plot_title,'_ILMAH.png'),final_plot_with_labels,width=24,height=12)


 