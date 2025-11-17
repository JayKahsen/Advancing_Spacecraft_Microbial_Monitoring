################################################################################
# clean work space, set a seed, avoid scientific nptation
################################################################################
rm(list = ls())
set.seed(20240725) # Date
options(scipen = 999)
################################################################################
# Description
################################################################################
# Global Script used as a starting point for other scripts
################################################################################

################################################################################
# set working directory from script location
# load helper file/functions
# load libraries
################################################################################
wd <-dirname(normalizePath(rstudioapi::getSourceEditorContext()$path))
project_folder <- basename(wd)
setwd(wd)

source(file.path(dirname(wd), "helperJ.R"))

library(tidyverse)
# included in tidyverse
# ggplot2,dplyr,tidyr,readr,purr,tibble,stringr,forcats 
library(seqinr)
#library(phyloseq)
library(vegan)
#library(glue)
library(gridExtra)
library(compositions)
#library(robCompositions)
library("ggh4x")
#library(coin)
#library(rlang)
#library(UpSetR)
library(scales)
library(grid)
library(ggtext)
library(patchwork)
conflicts()

conflicts()

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

################################################################################
# Load in data/ objects if exist
################################################################################

qLoad("data_tables/original/Contaminant_Id_Taxa.xlsx")

qLoad('data_tables/meta.csv') # meta data
qLoad('data_tables/matrix_names.csv') # used to select data sets to run scripts on
qLoad('data_tables/ASV_taxa.csv') # taxonomic info for each ASV may include sequence
if(exists('meta')){  
  meta_names=names(meta) # names of all meta columns
  qPrint(meta_names)
  
}

################################################################################
# taxa levels
# common data set descriptor incorporated into matrix_names 
# used to determine which levels scripts are run on
################################################################################
taxa_levels=c('Phylum','Order','Class','Family','Genus','Species','ASV','Contaminated')
taxa_plural=c('Phyla','Orders','Classes','Families','Genera','Species','ASVs','Contaminated')

taxa_levels=c('Phylum','Genus','ASV','Contaminated')
taxa_plural=c('Phyla','Genera','ASVs','Contaminated')

taxa_levels=c('Phylum','Genus','ASV')
taxa_plural=c('Phyla','Genera','ASVs')

if(exists('matrix_names')){  
  data_set_order=data_sets=unique(matrix_names$data_sets)
}

################################################################################
# defining general output folders, refine in scripts if necessary
################################################################################
output_data=paste0('output_data/')
output_plot=paste0('output_plot/')
################################################################################
# plotting orders
# factor levels to aid in plotting order
# data set dependent, but could be set with a generic script based on meta
################################################################################
meta=meta %>% 
  rename(location=Location) %>% 
  mutate(sample_time=Sampling.Time,
         treatment=Treatment,
         experiment='ILMAH'
  ) %>% 
  ungroup()

meta_names=names(meta)

vPrint(unique(meta$sample_time))
vPrint(unique(meta$location))
vPrint(unique(meta$treatment))

sample_time_order=c( 'Day 0','Day 5','Day 10','Day 15','Day 20')
treatment_order=c( 'PMA','No PMA')

location_order=c( 'Entrance','Kitchen','Bathroom Wall','Crew Sleeping Quarters (subbed for Laundry area)',
                  'Exercise Module','Medical Bay','Plant Production Module','EVA Module','Field Control',
                  'Negative Control','No PMA','PMA','NA','Zymo community standard','Extraction Blank PMA',
                  'Extraction Blank No PMA','no template control' )

selected_locations=c( 'Entrance','Kitchen','Bathroom Wall','Crew Sleeping Quarters (subbed for Laundry area)',
                      'Exercise Module','Medical Bay','Plant Production Module','EVA Module','Field Control')

location_order=selected_locations
experiment_order='ILMAH'


# cycle <- 1:24
# plot_cycle_order <- paste0('cycle_', str_pad(cycle, width = 2, side = "left", pad = "0"))

custom_name='' # initializing string that is incorporated into titles and saved names
################################################################################
# Groups
# Scripts run using different groups such as x_axis,y_axis, facets, plot, run
# could possible set them up with generic names that could be defined from meta
################################################################################
# primary group = separate plots for
# primary_group='experiment'
# primary_group_order<-experiment_order
# selected_groups=primary_group_order
# 
# primary_group='sample_type'
# primary_group_order<-sample_type_order
# selected_groups=primary_group_order
# 
# subgroup1=c('primer_set')
# subgroup1_order=primer_order
# 
# subgroup2=c('temperature')
# subgroup2_order=temperature_order
# 
# subgroup3=c('mastermix')
# subgroup3_order=mastermix_order
# #grouping_columns=c(primary_group,subgroup1)
# 
# grouping_columns=c(primary_group,subgroup1,subgroup2)

#data_set_filter_group='sample_type'
# data_sets=data_set_order=c('full_plates','half_plates')
# full_plates=c('G4','G6','G7')
# half_plates=c('G9','G10')

################################################################################
# palette
# defining palettes to control plotting aspects
# first entry takes precedence
################################################################################
palette_color <- c(    
  'Entrance'=palette_common[15],
  'Kitchen'=palette_common[16],
  'Bathroom Wall'=palette_common[3],
  'Crew Sleeping Quarters (subbed for Laundry area)'=palette_common[4],
  'Exercise Module'=palette_common[5],
  'Medical Bay'=palette_common[6],
  'Plant Production Module'=palette_common[7],
  'EVA Module'=palette_common[8],
  'Field Control'=palette_common[9],
  'Day 0'=palette_common[10],
  'Day 5'=palette_common[11],
  'Day 10'=palette_common[12],
  'Day 15'=palette_common[13],
  'Day 20'=palette_common[14],
  'PMA'=palette_common[1],
  'No PMA'=palette_common[2],
  'ILMAH'=palette_common[17],
  'NTC'='gray75',
  'unique_features'=palette_common[18],
  'chimeras'=palette_common[19],
  
  
  
  'ns'='gray50',
  'sig'='red',
  
  
  '1'='blue3',
  '2'='darkgreen',
  '3'='goldenrod',
  '4'='darkred',
  '5'='skyblue1',
  "6"="palegreen",
  "7" = "gold",
  '8'= "hotpink",
  '9'= "plum",  
  "Standard" = "red", "Truncated53" = "goldenrod1", "Truncated50" = "blue",
  "Archaea" = "#FF4040",
  "Bacteria" = "#000080",
  "Eukaryote"="darkgreen",
  "unknown"="gray30",
  'QB'='red3','AMPED'='green4',
  'slope'= "lightskyblue1",
  '1500TF'= "deepskyblue1",
  'x3FC'= "dodgerblue3",
  '24cycles'= "gray30",
  'slope'= "gray90",
  '1500TF'= "gray70",
  'x3FC'= "gray50",
  '24cycles'= "gray30",
  '1' = 'red',
  '2' = 'blue',
  '3' = 'green',
  '4' = 'purple',
  '5' = 'orange',
  '6' = 'pink',
  '7' = 'cyan',
  '8' = 'yellow',
  '9' = 'brown',
  '10' = 'magenta',
  '11' = 'limegreen',
  '12' = 'navy',
  '13' = 'darkolivegreen',
  '14' = 'darkcyan',
  '15' = 'gold',
  '16' = 'violet',
  '17' = 'coral',
  '18' = 'turquoise',
  '19' = 'salmon',
  '20' = 'maroon',
  '21' = 'orchid',
  '22' = 'darkgreen',
  '23' = 'steelblue',
  '24' = 'midnightblue',
  'cycle_01' = 'red',
  'cycle_02' = 'blue',
  'cycle_03' = 'green',
  'cycle_04' = 'purple',
  'cycle_05' = 'orange',
  'cycle_06' = 'pink',
  'cycle_07' = 'cyan',
  'cycle_08' = 'yellow',
  'cycle_09' = 'brown',
  'cycle_10' = 'magenta',
  'cycle_11' = 'limegreen',
  'cycle_12' = 'navy',
  'cycle_13' = 'darkolivegreen',
  'cycle_14' = 'darkcyan',
  'cycle_15' = 'gold',
  'cycle_16' = 'violet',
  'cycle_17' = 'coral',
  'cycle_18' = 'turquoise',
  'cycle_19' = 'salmon',
  'cycle_20' = 'maroon',
  'cycle_21' = 'orchid',
  'cycle_22' = 'darkgreen',
  'cycle_23' = 'steelblue',
  'cycle_24' = 'midnightblue',
  'white'='white')
################################################################################

qnvPrint(unique(meta$sample_time))
qnvPrint(unique(meta$location))
qnvPrint(location_order)
qnvPrint(unique(meta$treatment))



palette_label <- c(   
  'ILMAH'='ILMAH',
  'Day 0'='Day 0',
  'Day 5'='Day 5',
  'Day 10'='Day 10',
  'Day 15'='Day 15',
  'Day 20'='Day 20',
  'NA'='NA',
  'Entrance'='Entrance',
  'Kitchen'='Kitchen',
  'Bathroom Wall'='Bathroom Wall',
  'Crew Sleeping Quarters (subbed for Laundry area)'='Crew Sleeping Quarters',
  'Exercise Module'='Exercise Module',
  'Medical Bay'='Medical Bay',
  'Plant Production Module'='Plant Production Module',
  'EVA Module'='EVA Module',
  'Field Control'='Field Control', 
  'Negative Control'='Negative Control',
  'No PMA'='No PMA',
  'PMA'='PMA',
  'NA'='NA',
  'Zymo community standard'='Zymo community standard',
  'Extraction Blank PMA'='Extraction Blank PMA',
  'Extraction Blank No PMA'='Extraction Blank No PMA',
  'no template control'='no template control', 
  'PMA'='PMA',
  'No PMA'='No PMA',
  'NA'='NA', 
  'unique_features'='Unique Features',
  'chimeras'='Chimeras',
  'both'='Shared',
  "A"="A",
  "C"="C",
  "G"="G",
  "T"="T",
  'highlight_on'='Real data',
  'highlight_off'='Primer',
  'g__Cutibacterium'='Cutibacterium',
  'no'='not detected',
  'G11'='G11',
  'Standard'='Standard',
  'Truncated50'=  'Truncated50',
  'Truncated53'='Truncated53',
  '42°C'='42°C',
  '47°C'='47°C',
  '52°C'='52°C',
  '57°C'='57°C',
  'Skin'='Skin',
  'Feces'='Feces',
  'Soil'='Soil',
  "NTC" = "NTC",
  'AMPED'='AMPED',
  'QB'='QB',
  'Slope'='Slope','Targeted Fluorescence'='Targeted Fluorescence','Fold Change'='Fold Change','Fixed Cycles'='Fixed Cycles',
  'slope'='Slope','targeted fluorescence'='Targeted Fluorescence','fold change'='Fold Change','fixed cycles'='Fixed Cycles',
  "Soil" = "Soil", "Zymo" = "Zymo", "NTC" = "NTC",
  'QB'='QB','AMPED'='AMPED',
  'Slope'='Slope','Targeted Fluorescence 1500'='Targeted','Fold Change x3'='Fold',
  "Standard" = "Standard", "Truncated53" = "Truncated53", "Truncated50" = "Truncated50",
  'slope'= "Slope",
  '1500TF'= "Targeted",
  'x3FC'= "Fold",
  '24cycles'= "Fixed",
  'percent_Muri'='% Muribaculaceae',
  'percent_Lachno'='% Lachnospiraceae',
  'percent_Propi'='% Cutibacterium',
  'percent_Archaea'='% Archaea',
  'percent_good'='% Good',
  'richness'='Richness',
  'evenness'='Evenness',
  'shannon'='Shannon',
  'n'='Reads', 
    'white'='white')
################################################################################
palette_shape <- c(
  
  'PMA'=23,
  'No PMA'=21,
  
  'Slope'=2,'Targeted Fluorescence 1500'=3,'Fold Change x3'=0,
  'Standard'=22,
  'Truncated50'=  24,
  'Truncated53'=23,
  '42°C'=22,
  '47°C'=22,
  '52°C'=22,
  '57°C'=22,
  # "AMPED" = "A",
  # "QB" = "Q",
  'slope'= 23,
  '1500TF'= 21,
  'x3FC'= 24,
  '24cycles'= 22,
  "AMPED" = 17,
  "QB" = 15
)
################################################################################
palette_size <- c(  
  
  'PMA'=4,
  'No PMA'=4,
  'Standard'=3,
  'Truncated50'=3,
  'Truncated53'=3,
  '42°C'=3,
  '47°C'=3,
  '52°C'=3,
  '57°C'=3,
  "AMPED" = 3,
  "QB" = 3,
  'slope'= 5,
  '1500TF'= 5,
  'x3FC'= 5,
  '24cycles'= 5,
  'Standard'= 3,
  'Truncated'= 3,
  'TruncatedR'= 3,
  'TruncatedF'= 3,
  "V1" = 3,
  "V2" = 3,
  "V3" = 3,
  "V4" = 3,
  "V5" = 3,
  "V6" = 3,
  "V7" = 3,
  "V8" = 3,
  'StandardO'= 3,
  'TruncatedFO'= 3,
  'Standard'= 3,
  'TruncatedR'= 3,
  'TruncatedF'= 3,
  'TruncatedA'= 3,
  'All'= 3,
  'Pool'= 3,
  '515F.806R'=3,
  'CasF.806R'=3,
  '515F.CasR'=3,
  'CasF.CasR'=3,
  "shotgun" = 5,
  "515F" = 5,
  "Cascading" = 5,
  "Sampled" = 3,
  "V1.2-12.CC" = 3,
  "V2.2-13.CA" = 3,
  "V3.3-14.TC" = 3,
  "V4.4-15.TA" = 3,
  "V5.5-16.xA" = 3,
  "V6.6-15.xC" = 3,
  "V7.7-17.xA" = 3,
  "V8.8-17.xC" = 3,
  "shotgun_group" = 5,
  "Archaea"=5,
  "Bacteria"=3,
  "unknown"=3,
  "515F_group"=5,
  "Cascading_group"=5,
  "Sampled_group"=3,
  "Ind_group"=3)
################################################################################
palette_linetype <- c(  
  'Standard'='solid',
  'Truncated50'=  'dotted',
  'Truncated53'='dashed',
  '42°C'='42°C',
  '47°C'='47°C',
  '52°C'='52°C',
  '57°C'='57°C',
  'Skin'='Skin',
  'Feces'='Feces',
  'Soil'='Soil',
  'AMPED'='AMPED',
  "NTC" = "NTC",
  'Soil' = 'solid', 'Zymo' = 'dashed', 'NTC' = 'dotted'
)
################################################################################
# Zymo taxonomy for comparison if using Zymo 
################################################################################
Zymo_names=c(
  'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Listeriaceae;g__Listeria',
  'd__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus',
  'd__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus',
  'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus',
  'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus',
  'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia-Shigella',
  'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;__',
  'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas'
)


################################################################################
# Functions
################################################################################
################################################################################
# common Theme and gPlot
# plotting function to add palettes, comoon theme, and ?
# refine more in individual scripts
################################################################################
# library(RColorBrewer)
# brewer_colors <- brewer.pal(8, "Dark2")
# scale_color_manual(values = brewer_colors)

magnify=1
text_size=10
# Define text size variables
axis_title_size <- (2+text_size)*magnify
axis_text_size <- (0+text_size)*magnify
strip_text_size <- (4+text_size)*magnify
legend_text_size <- (2+text_size)*magnify
caption_text_size <- (0+text_size)*magnify
plot_title_size <- (4+text_size)*magnify
subtitle_size <- (2+text_size)*magnify

# Define color variables
background_color <- palette_color['general_background']
panel_background_color <- "white"
grid_major_color <- "grey90"
grid_minor_color <- "grey95"
strip_background_color <- "grey90"
axis_line_color <- "black"
axis_tick_color <- "black"
legend_key_color <- "white"

common_theme <- theme(
  # Plot layout & background
  #  plot.margin = margin(t=5.5, r=5.5, b=45, l=5.5),
  #  plot.margin = margin(5.5, 5.5, 45, 5.5),
  #  plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
  plot.background = element_rect(fill = background_color),
  #  panel.background = element_rect(fill = panel_background_color),
  #  panel.grid.major = element_line(color = grid_major_color),
  #  panel.grid.major.x = element_blank(),
  
  #  panel.grid.minor = element_line(color = grid_minor_color),
  # panel.grid.minor.x = element_blank(),
  #  panel.spacing = unit(0.5, "cm"),
  
  # Titles & captions
  plot.title = element_text(size = plot_title_size, face = "bold", hjust = 0.5),
  #  plot.title = element_blank(),
  plot.subtitle = element_text(size = subtitle_size, hjust = 0.5),
  #plot.caption = element_text(size = caption_text_size, hjust = 0.5),
  
  # Legend
  legend.background = element_rect(fill = background_color),
  legend.key = element_rect(fill = legend_key_color),
  #  legend.position = "bottom",
  #  legend.direction = "vertical",
  #  legend.position = c(0.5, -.09),
  legend.title = element_text(size = legend_text_size, face = "bold"),
  #  legend.title = element_blank().
  legend.text = element_text(size = legend_text_size),
  # legend.margin = margin(5, 5, 5, 5),
  
  # Axes customization
  #  axis.title.y = element_text(face = 'bold', size = axis_title_size, margin = margin(r = 10)),
  axis.title.y = element_text(face = 'bold', size = axis_title_size),
  #  axis.title.x = element_text(face = 'bold', size = axis_title_size, margin = margin(t = 10)),
  axis.title.x = element_text(face = 'bold', size = axis_title_size),
  axis.text.x = element_text(angle = 90, size = axis_text_size,vjust=1,hjust=.25),  
  axis.text.y = element_text(size = axis_text_size),
  #  axis.text.x = element_blank(),
  axis.ticks = element_line(color = axis_tick_color),
  #  axis.ticks.x = element_blank(),
  #  axis.line = element_line(color = axis_line_color),
  
  # Facet labels (strip text)
  strip.text = element_text(face = 'bold', size = strip_text_size),
  #  strip.background = element_rect(fill = strip_background_color),
  plot.caption = element_text(hjust = 0.5)
  
)
################################################################################
# colored x axis labels
################################################################################
# scale_x_discrete(labels = function(x) {
#   color <- color_palette[x]  
#   label <- label_palette[x]
#   #paste0("<span style='color:", color, "'>", x, "</span>")
#   paste0("<span style='color:", color, "'>", label, "</span>")
# }) +
#   theme(
#     axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5,face='bold',size=x_text_size)  # Rotate and enable markdown for x-axis labels
#   )+
################################################################################
gPlot <- function(p) {
  p=p+
    #   scale_fill_manual(values = color_palette) +
    scale_color_manual(values = palette_color,labels=palette_label)+
    scale_fill_manual(values = palette_color,labels=palette_label)+
    scale_size_manual(values = palette_size,labels=palette_label)+
    scale_shape_manual(values = palette_shape,labels=palette_label)+
    scale_linetype_manual(values = palette_linetype,labels=palette_label)+
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
  common_theme
  print(p)
  return(p)
}
################################################################################
# axis modifying functions
################################################################################
# log_breaks <- 10^seq(0, 12, by = 2)
# 
# log_labels <- function(x) {
#   sapply(x, function(val) {
#     if (is.na(val)) {
#       return(NA)  # Return NA if the value is missing
#     } else if (val == 1) {
#       return(expression(10^0))
#     } else {
#       return(as.expression(bquote(10^.(round(log10(val))))))
#     }
#   })
# }
################################################################################
# strip backgrounds with outlines
# using ggh4x and facet_wrap2/ facet_grid2 can control facet strip colors and text,
# method is a bit clunky, may develop function at some point, currently modified in scripts
################################################################################
# unique_x_labels=unique_Primers <- intersect(primer_order, unique(df$Primer))
# 
# outline_color_x='black'
# 
# s=paste0('backgrounds_x <- list(element_rect(color=outline_color_x,fill = color_palette["',paste0(unique_x_labels,collapse='"]),element_rect(color=outline_color_x,fill = color_palette["'),'"]))')
# print(s)
# eval(parse(t=s))
# backgrounds_x   
# 
# unique_y_labels <- intersect(Type_order, unique(df$Type))
# 
# outline_color_y='black'
# 
# s=paste0('backgrounds_y <- list(element_rect(color=outline_color_y,fill = color_palette["',paste0(unique_y_labels,collapse='"]),element_rect(color=outline_color_y,fill = color_palette["'),'"]))')
# eval(parse(t=s))

################################################################################
# facet_grid/wrap2 function
################################################################################
# facet_grid2(Type ~ Primer,
#             labeller = labeller(Primer = function(value) df$primer_Label[match(value, df$Primer)]),
#             strip = strip_themed(
#               background_x = backgrounds_x, background_y =backgrounds_y,
#               text_x=elem_list_text(
#                 color = c('white',rep('black',(length(unique_x)-1))),
#                face=c(rep('bold',length(unique_x))),
#                 size=c(rep(strip_text_size,length(unique_x)))),
#               text_y=elem_list_text(
#                 color = c('white',rep('black',(length(unique_y)-1))),
#                 face=c(rep('bold',length(unique_y))),
#                 size=c(rep(strip_text_size,length(unique_y))))),
#             scale = 'free_y') +
################################################################################
# removes row names, adds meta data, adjusts data frame
################################################################################
xPlode_sample_name <- function(df) {
  df=df %>%
    as.data.frame() %>%
    rownames_to_column(var='sample_name')%>%
    left_join(meta)%>%
    select(any_of(meta_names),everything())%>%
    # mutate(primer_Label = factor(primer_Label, levels = unique(meta$primer_Label)))%>%
    # mutate(Primer = factor(Primer, levels = unique(meta$Primer)))%>%
    
    ungroup()
  
}

################################################################################
# adds row names,removes meta data
################################################################################
imPlode_sample_name <- function(df) {
  df %>%
    as.data.frame() %>%
    column_to_rownames(var='sample_name')%>%
    select(-any_of(meta_names))
}

################################################################################
# script has been refined and moved into helperJ feature_labels

# makes shortened labels from the taxonomy ($taxa)
# some label sets are colored by Domain
# Domain is extracted from taxonomy
# relative abundance is calculated across ASV_labels(group) and enetered as a percentage into some labels

################################################################################
# 
# ASV_labels<-function(group='ASV',df=matrix_df){
#   # ASV_labels_df=ASV_labels('group',df)
#   ASV_names=names(matrix_df)
#   if(!('ASV' %in% group)){group=c(group,'ASV')}
#   
#   df=df%>%
#     as.data.frame()%>%
#     xPlode_sample_name()%>%
#     pivot_longer(col=all_of(ASV_names),names_to='ASV',values_to='counts')%>%
#     mutate(taxa=ASV)%>%
#     mutate(Domain = sapply(taxa, extract_domain))%>%
#     mutate(ASV_label = sub(".*?d__.*?;(.*)", "\\1", ASV))%>%
#     mutate(ASV_label=if_else(ASV_label=='__','unclassified',ASV_label))%>%
#     mutate(ASV_label = sub("^p__","", ASV_label))%>%
#     group_by(sample_name)%>%
#     mutate(sample_counts=sum(counts))%>%
#     ungroup()%>%
#     mutate(rel_abun=counts/sample_counts)%>%
#     group_by(ASV,ASV_label,Domain)%>%
#     mutate(group_rel_abun=rel_abun)%>%
#     group_by(across(all_of(group)),ASV_label,Domain)%>%
#     summarize(group_rel_abun=mean(rel_abun))%>%
#     ungroup()%>%
#     mutate(rel_abun_label = paste0(signif(100*group_rel_abun,2),'%'))%>%
#     mutate(ASV_label2 = paste0("<span style='color:", color_palette[Domain], "'><i>", ASV_label, "</i>"," ",rel_abun_label, "</span>"))%>%
#     mutate(ASV_label3 = paste0(ASV_label,' ',signif(100*group_rel_abun,2),'%'))%>%
#     #distinct()%>%
#     arrange(desc(group_rel_abun)) %>%
#     mutate(plot_order = factor(ASV_label3, levels = rev(unique(ASV_label3))))
#   #color <- color_palette[domain]
#   #   paste0("<span style='color:", color, "'><i>", asv_label, "</i>"," ",rel_abund_label, "</span>")
# }
################################################################################

################################################################################
# END
################################################################################
################################################################################
# get citations, might have a function to curate this info
################################################################################

# Citation for R
citation()

# Citation for compositions
citation("compositions")

# Citation for vegan
citation("vegan")

# Citation for tidyverse
citation("tidyverse")

# Citation for ggplot2
citation("ggplot2")

# Citation for ggtext
citation("ggtext")