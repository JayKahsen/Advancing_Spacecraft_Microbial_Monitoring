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
#source(file.path(wd, "helperJ.R"))
library(patchwork)
qPrint(project_folder)

library(tidyverse)
#library(seqinr)
#library(phyloseq)
#library(tibble)
library(vegan)

#library(glue)
library(gridExtra)
#library(dplyr)
#library(tidyr)
library(stringr)
#library(patchwork)
#library(purrr)
library(compositions)
#library(robCompositions)
library("ggh4x")
#library(coin)
#library(rlang)
#library(UpSetR)
library(scales)
#library(ggplot2)
library(grid)
library(ggtext)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(broom)
library(ggpattern)
library(magick)
library(ggsignif)
conflicts()
version
# install.packages("rlang")
# install.packages("rlang", type = "source")
# 
# Sys.which("make")



# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

################################################################################
# Load in data/ objects if exist
################################################################################
current_date <- format(Sys.Date(), "%Y%m%d")

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

taxa_levels=c('Phylum','Family','Genus','Species')
taxa_plural=c('Phyla','Families','Genera','Species')

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
# sampling_device_order=c( 'SALSA','Polyester Wipe','cotton','macrofoam','Positive Control' )
# sample_type_order=c('Positive Control','Negative Control','Reagent Control',
#                     '2.00E+2','2.00E+3','2.00E+4','2.00E+5','2.00E+6','NTC' )
#extraction_method_order=c('magnetic_bead','spin_column')

dpi=300

treatment_order=c('No PMA','PMA')
extraction_method_order=c( 'mag_beads','spin_column' )
location_order=c("PHSF Highbay floor","Field Control","PHSF Airlock floor")
experiment_order=c('PPR','Clipper')
experiment_order=c('PPR')
sample_type_order=c( 'reagent_control','2.00E+1','2.00E+2','2.00E+3','2.00E+4','2.00E+5','2.00E+6','2.00E+7','water_control','positive_control','no_template_control','100000','NTC','NA' )
vPrint(unique(meta$sample_type))

vPrint(unique(meta$bar_color))
vPrint(unique(meta$bar_name))
bar_name_order=c('Cotton Device Control' ,'Macrofoam Device Control','Cotton Environmental Control','Macrofoam Environmental Control' ,'Cotton Dry + Metal','Macrofoam Dry + Metal','Cotton Wet + Metal','Macrofoam Wet + Metal' ,'Water Only','NTC')
bar_color_order=c('Cotton','Macrofoam','Water','NTC' )
bar_axis_order=c('Cotton Wet','Macrofoam Wet','Cotton Dry','Macrofoam Dry','Cotton','Macrofoam','Water','NTC' )
bar_axis_2_order=c('Cotton','Macrofoam','Water','NTC' )

sampling_device_order=c( 'cotton','macrofoam','Polyester Wipe','SALSA')
bar_type_order=c('Mastermix','Environmental','Device','Metal' )
bar_type_2_order=c('Mastermix','Environmental','Device','Metal Dry','Metal Wet')
method_order=c('qpcr','dpcr')
experiment_type_order=c( 'metal_deposition','swab_head_retention' )
figure_order=c('2','3A','3B','4A','4B','4C','5')


order_vector=unique(c(treatment_order,experiment_order,location_order,bar_color_order,extraction_method_order,sampling_device_order,bar_type_order,method_order,experiment_type_order))

vPrint(unique(meta$bar_type))
vPrint(unique(meta$sampling_device))
vPrint(unique(meta$bar_color))
vPrint(unique(meta$extraction_method))
vPrint(unique(meta$location))
vPrint(unique(meta$treatment))
vPrint(unique(meta$experiment))
vPrint(unique(meta$sample_type))
vPrint(unique(meta$experiment_type))
#vPrint(unique(df$sampling_device))
custom_name='' # initializing string that is incorporated into titles and saved names
################################################################################
# Groups
# Scripts run using different groups such as x_axis,y_axis, facets, plot, run
# could possible set them up with generic names that could be defined from meta
################################################################################

################################################################################
# palette
# defining palettes to control plotting aspects
# first entry takes precedence
################################################################################

################################################################################
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

atcc_genera=genera <- c(
  "Bacillus",
  "Bifidobacterium",
  "Clostridium",
  "Deinococcus",
  "Enterococcus",
  "Escherichia",
  "Lactobacillus",
  "Cereibacter",
  "Staphylococcus",
  "Streptococcus"
  # "g_Bacillus",
  # "g_Bifidobacterium",
  # "g_Clostridium",
  # "g_Deinococcus",
  # "g_Enterococcus",
  # "g_Escherichia",
  # "g_Lactobacillus",
  # "g_Cereibacter",
  # "g_Staphylococcus",
  # "g_Streptococcus"
)

################################################################################
# Functions
################################################################################

################################################################################
library(ggplot2)
library(ggpattern)

################################################################################
################################################################################
################################################################################
################################################################################
log_limits=c(NA, NA)
zoom=c(NA, NA)
point_size=2.5
bar_color='gray30'
################################################################################



#qnvPrint(sample_type_order)

# qnvPrint(unique(df_plot$extraction_method))

#qnvPrint(unique(df_plot$sample_type))


qnvPrint(order_vector)
# qnvPrint(metric_order)

palette_label <- c(   
  '<i>Escherichia</i> '='<i>Escherichia</i> ',
  '<i>Cereibacter</i> '='<i>Cereibacter</i> ',
  '<i>Bacillus</i> '='<i>Bacillus</i> ',
  '<i>Lactobacillus</i> '='<i>Lactobacillus</i> ',
  '<i>Clostridium</i> '='<i>Clostridium</i> ',
  '<i>Streptococcus</i> '='<i>Streptococcus</i> ',
  '<i>Enterococcus</i> '='<i>Enterococcus</i> ',
  '<i>Bifidobacterium</i> '='<i>Bifidobacterium</i> ',
  '<i>Staphylococcus</i> '='<i>Staphylococcus</i> ',
  '<i>Deinococcus</i> '='<i>Deinococcus</i> ',
  'Other'='Other', 

  'reagent_control'='Reagent<br>Control',
  'mag_beads'='Automatic Magnetic Beads',
  'spin_column'='Manual Spin Column',
  '5A'='titleA?',
  '5B'='titleB?',
  'ng100'='<100',
  'g100'='>100',
  'g1K'='>1K',
  'g10K'='>10K',
  'g100K'='>100K',
  
  'metal_deposition'='Recovery from Metal Surface',
  'swab_head_retention'='Swab Head Released' ,
  'unique_features'='unique_features',
  'chimeras'='chimeras',
  'percent_Muri'='% Muri',
  'percent_Lachno'='% Lachno',
  'percent_Propi'='% Propi',
  'percent_Archaea'='% Archaea',
  'percent_good'='% Good',
  'ideal_score'='ideal_score',
  'richness'='Richness',
  'evenness'='Evenness',
  'shannon'='Shannon',
  'n'='Reads',
  'less_than_0'='0',
  'log0'='1',
  'log4'='e4',
  'log6'='e6',
  'log8'='e8',
  'log10'='e10',
  '0P'='0%',
  '20P'='10%',
  '40P'='20%',
  '60P'='40%',
  '80P'='80%',
  'g0'='0',
  'g10'='10',
  'g20'='20',
  'g40'='40',
  'g80'='80',
  
  '4A'='A\\) DNA Released from Swab Heads',
  '4B'='B\\) DNA Recovered by Swabs',
  '4C'='C\\) DNA Rec. by Polyester Wipes and SALSA',
  '4B'='Microbial DNA Recovery Efficiency<br>
    of Swabs on Spiked Coupons (cm<sup>2</sup>)',
  '4C'='Microbial DNA Recovery Efficiency<br>
    on Spiked Coupons (25 cm<sup>2</sup>)',
  'No PMA'='No PMA',
  'PMA'='PMA',
  'PPR'='PPR',
  'Clipper'='Clipper',
  'PHSF Highbay floor'='PHSF Highbay floor',
  'Field Control'='Field Control',
  'PHSF Airlock floor'='PHSF Airlock floor',
  'Cotton Wet'='Cotton Wet',
  'Macrofoam Wet'='Macrofoam Wet',
  'Cotton Dry'='Cotton Dry',
  'Macrofoam Dry'='Macrofoam Dry',
  'Cotton'='Cotton',
  'Macrofoam'='Macrofoam',
  'Water'='Water',
  'NTC'='NTC',
  'positive_control'='Positive Control',
  'negative_control'='Negative Control',

  'water_control'='Water Control',
  'non_template_control'='NTC',
  'no_template_control'='NTC',

  #'cotton'='Cotton Swab',
  'cotton'='Cotton',
  #'macrofoam'='Macrofoam Swab',
  'macrofoam'='Macrofoam',
  'Polyester Wipe'='Polyester Wipe',
  'SALSA'='SALSA',
  'Mastermix'='Mastermix',
  'Environmental'='Environmental',
  'Device'='Swab',
  'Metal Dry'='Metal Dry',
  'Metal Wet'='Metal Wet',
  'Metal'='Metal',
  'dpcr'='dPCR',
  'qpcr'='qPCR',
  
  
  
  
  
  
  'Clipper'='Clipper',
  'PHSF Highbay floor'='PHSF Highbay floor',
  'Field Control'='Field Control',
  'PHSF Airlock floor'='PHSF Airlock floor',
  'ATCC'='ATCC',
  'No PMA'='No PMA',
  'PMA'='PMA', 
  
  'Cotton'='Cotton',
  'Water'='Water',
  'Macrofoam'='Macrofoam',
  'Cotton Device Control'='Cotton Device Control',
  'Macrofoam Device Control'='Macrofoam Device Control',
  'Cotton Environmental Control'='Cotton Environmental Control',
  'Macrofoam Environmental Control'='Macrofoam Environmental Control',
  'Cotton Dry + Metal'='Cotton Dry + Metal',
  'Macrofoam Dry + Metal'='Macrofoam Dry + Metal',
  'Cotton Wet + Metal'='Cotton Wet + Metal',
  'Macrofoam Wet + Metal'='Macrofoam Wet + Metal',
  'Water Only'='Water Only',
  'NTC'='NTC', 
  'PMA'='PMA',
  'No PMA'='No PMA',
  na.value='NTC',
  'Positive Control'='Positive Control',
  'negative_control'='Negative Control',
  'positive_control'='Positive Control',
  'Negative Control'='Negative Control',
  'Reagent Control'='Reagent\nControl',

  '2.00E+2'='2.00E+2',
  '2.00E+3'='2.00E+3',
  '2.00E+4'='2.00E+4',
  '2.00E+5'='2.00E+5',
  '2.00E+6'='2.00E+6',
  'NTC'='NTC', 
  'SALSA'='SALSA',
  'Polyester Wipe'='Polyester Wipe',
 # 'cotton'='Cotton Swab',
 # 'macrofoam'='Macrofoam Swab',
  'ATP_RLU'='ATP RLU',
  'Fungal_CFU'='Fungal CFU',
  'Fungal_dPCR'='Fungal dPCR',
  'Bacterial_CFU'='Bacterial CFU',
  'Bacterial_dPCR'='Bacterial dPCR',
  'default'='default',
  'NA'='NTC',
  'subject_id'='Subject ID',
  'gmcf_number'='GMCF Number',
  'sample_number'='Sample Number',
  'treatment'='Treatment',
  'location'='Location',
  'timepoints'='Timepoints (day)',
  'dna_concentration)'='DNA  Conc. (ng/uL)',
  'miniseq_pf_clusters'='Miniseq PF Clusters',
  'miniseq_p_of_lane'='Miniseq  % of Lane',
  'novaseq_pf_clusters'='NovaSeq X PF Clusters',
  'novaseq_x_p_of_lane'='NovaSeq X % of Lane',
  'Weight_of_liquid'='Weight of liquid after InnovaPrep concentration (1m2 of surface area collected)',
  'bacterial_concentration'='Bacterial Conc. [cp/¬µL] (dPCR reaction)',
  'bacterial_partitions_valid'='Bacterial  Partitions (Valid)',
  'bacterial_partitions_positive'='Bacterial  Partitions (Positive)',
  'bacterial_copies'='Bacterial  Copies in Orginal Sample (2uL)',
  'bacterial_final_dilution_factor'='FINAL Dilution Factor...17',
  'bacterial_copies_in_orginal_sample'='Copies in Orginal Sample (50uL)[dPCR copies*dil*50*3]...18',
  'bacterial_digital_pcr_conversion'='Bacterial Digital PCR conversion (2uL to 1m2) dPCR*50*dil*3*[Column L/0.375] or [Column L/0.3775]',
  'blank'='...20',
  'fungal_concentration'='Fungal Conc. [cp/¬µL] (dPCR reaction)',
  'fungal_partitions_valid'='Fungal  Partitions (Valid)',
  'fungal_partitions_positive'='Fungal Partitions (Positive)',
  'fungal_copies'='Fungal Copies in Orginal Sample (2uL)',
  'fungal_final_dilution_factor'='FINAL Dilution Factor...25',
  'fungal_copies_in_orginal_sample'='Copies in Orginal Sample (50uL)[dPCR copies*dil*50*3]...26',
  'fungal_digital_pcr_conversion'='Fungall Digital PCR conversion (2uL to 1m2) dPCR*50*[Column L/0.375] or [Column L/0.3775]',
  'Entrance'='Entrance',
  'Kitchen'='Kitchen',
  'Bathroom'='Bathroom',
  'Crew Quarters'='Crew Quarters',
  'Exercise Module'='Exercise Module',
  'Medical Bay'='Medical Bay',
  'Plant Production'='Plant Production Module',
  'EVA Module'='EVA Module',
  'Field Control'='Field Control',
  'Neg Control'='Neg Control', 
  
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
  'No PMA'='No PMA',
  'PMA'='PMA',
  'PPR'='PPR',
  'PHSF Highbay floor'='PHSF Highbay floor',
  'Field Control'='Field Control',
  'PHSF Airlock floor'='PHSF Airlock floor',
  'Cotton'='Cotton',
  'Macrofoam'='Macrofoam',
  'Water'='Water',
  'NTC'='NTC',


  'cotton'='Cotton',
  'macrofoam'='Macrofoam',
  'Polyester Wipe'='Polyester Wipe',
  'SALSA'='SALSA',
  'Mastermix'='Mastermix',
  'Environmental'='Environmental',
  'Device'='Swab',
  'Metal'='Metal',
  'Metal Dry'='Metal Dry',
  'Metal Wet'='Metal Wet',
  'dpcr'='dpcr',
  'qpcr'='qpcr', 
  'both'='Shared',
  'gray1'='gray1',
  'gray2'='gray2',
  'gray3'='gray3',
  'gray4'='gray4',
  'gray5'='gray5',
  'black'='black',
  palette_helper_label)

palette_gray <- (c('#BBBBBB', '#999999', '#777777', '#555555','#333333'))
palette_green <- (c('#D2EDB2', '#89C189', '#60A740', '#107810', '#004B00'))
palette_gold  <- (c('#FFDA65', '#FFCB25', '#DAA520', '#A37B18', '#7D6220'))
palette_purple <- (c('#D2AEFA', '#A267E8', '#7A1FD2', '#5B17A1', '#3C1070'))
palette_red <- (c('#FFCCCC', '#E47C7C', '#C94141', '#993333', '#6A2B2B'))
palette_blue <- (c('#C8D6FF', '#92AEFF', '#5C86FF', '#3A62BE', '#20407F'))

assign_palette(order_vector)

palette_color=c(
  '<i>Escherichia</i> ' = palette_common[1],
  '<i>Cereibacter</i> ' = palette_common[2],
  '<i>Bacillus</i> ' = palette_common[3],
  '<i>Lactobacillus</i> ' = palette_common[4],
  '<i>Clostridium</i> ' = palette_common[5],
  '<i>Streptococcus</i> ' = palette_common[6],
  '<i>Enterococcus</i> ' = palette_common[7],
  '<i>Bifidobacterium</i> ' = palette_common[8],
  '<i>Staphylococcus</i> ' = palette_common[9],
  '<i>Deinococcus</i> ' = palette_common[10],
  'Other' = palette_gray[1],
  'metal_deposition' = palette_gold[4],
  'swab_head_retention' = palette_blue[4],
  'spin_column'=palette_blue[3],
  'mag_beads'=palette_red[4],
  'macrofoam'=palette_green[2],
  'Macrofoam'=palette_green[2],
  'cotton'=palette_purple[1],
  'Cotton'=palette_purple[1],
  'Water'=palette_blue[4],
  'Salsa'=palette_red[4],
  'SALSA'=palette_red[4],
  'Polyester Wipe'=palette_green[4],
  'PMA'=palette_blue[4],
  'No PMA'=palette_red[4],
  'positive_control'=palette_blue[2],
  'negative_control'=palette_red[2],
  'reagent_control'=palette_purple[5],
  na.value='grey40',
  'NTC'='grey40',
  'magnetic_bead'=palette_gold[5],
  'spin_column'=palette_gold[2],

  'alt_group1'=palette_gray[1],
  'alt_group2'=palette_gray[3],
  'Clipper' = palette_common[1],
  'PHSF Highbay floor' = palette_common[2],
  'Field Control' = palette_common[3],
  'PHSF Airlock floor' = palette_common[4],
  'ATCC' = palette_common[5],
  '2.00E+2'=palette_common[22],
  '2.00E+3'=palette_common[23],
  '2.00E+4'=palette_common[24],
  '2.00E+5'=palette_common[25],
  '2.00E+6'=palette_common[26],
  'gray1'=palette_gray[1],
  'gray2'=palette_gray[2],
  'gray3'=palette_gray[3],
  'gray4'=palette_gray[4],
  'gray5'=palette_gray[5],

  'Polyester Wipe'=palette_green[3],
  'Cotton Wet'=palette_purple[5],
  'Macrofoam Wet'=palette_blue[3],
  'Cotton Dry'=palette_purple[5],
  'Macrofoam Dry'=palette_blue[3],
  'Cotton'=palette_purple[5],
  'Macrofoam'=palette_blue[3],
  'Positive Control'=palette_common[9],

  'No PMA' = palette_common[1],
  'PMA' = palette_common[2],
  'PPR' = palette_common[3],
  'PHSF Highbay floor' = palette_common[4],
  'Field Control' = palette_common[5],
  'PHSF Airlock floor' = palette_common[6],
  'Cotton' = palette_common[7],
  'Macrofoam' = palette_common[8],
  'Water' = palette_common[9],
  'NTC' = palette_common[10],
  'mag_beads' = palette_common[11],
  'spin_column' = palette_common[12],
  'cotton' = palette_common[13],
  'macrofoam' = palette_common[14],
  'Polyester Wipe' = palette_common[15],
  'SALSA' = palette_common[16],
  'Mastermix' = palette_common[17],
  'Environmental' = palette_common[18],
  'Device' = palette_common[19],
  'Metal' = palette_common[20],
  'Metal Dry'=palette_common[20],
  'Metal Wet'=palette_common[20],
  'dpcr' = palette_common[21],
  'qpcr' = palette_common[22],
  
  
  
  
'ns'='gray50',
'sig'='red',
"Archaea" = "#FF4040",
"Bacteria" = "#000080",
"Eukaryote"="darkgreen",
"unknown"="gray30",
is.na='gray',
'black'='black',
palette_helper_color)
# base off above
palette_color=c(palette_color,
                'Cotton Device Control'=palette_color['Cotton'],
                'Macrofoam Device Control'=palette_color['Macrofoam'],
                'Cotton Environmental Control'=palette_color['Cotton'],
                'Macrofoam Environmental Control'=palette_color['Macrofoam'],
                'Cotton Dry + Metal'=palette_color['Cotton'],
                'Macrofoam Dry + Metal'=palette_color['Macrofoam'],
                'Cotton Wet + Metal'=palette_color['Cotton'],
                'Macrofoam Wet + Metal'=palette_color['Macrofoam'],
                'Water Only'=palette_color['Water'],
                #'NTC'='NTC', 
                
                'ATP_RLU'='#DC0000FF',
                'Fungal_CFU'='#3C5488FF',
                'Fungal_dPCR'='#4DBBD5FF',
                'Bacterial_CFU'='#00A087FF',
                'Bacterial_dPCR'='#91D1C2FF',
                'ATP_RLU'='#DC0000FF',
                'Fungal_CFU'='#3C5488FF',
                'Fungal_dPCR'='#4DBBD5FF',
                'Bacterial_CFU'='#00A087FF',
                'Bacterial_dPCR'='#91D1C2FF',
                'Bacterial_dPCR'='green2',
                'Fungal_dPCR'='darkturquoise',
                'Bacterial_CFU'='green4',
                'Bacterial_dPCR'='green',
                'Fungal_CFU'='dodgerblue4',
                'Fungal_dPCR'='dodgerblue',
                'Bacterial_dPCR'='yellow',
                'Fungal_dPCR'='darkturquoise',
                'Bacterial_CFU'='green4',
                'Bacterial_dPCR'='green',
                'Fungal_CFU'='dodgerblue4',
                'Fungal_dPCR'='dodgerblue',
                
                'ATP_RLU'='red3',
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
                'No PMA' = palette_common[1],
                'PMA' = palette_common[2],
                'PPR' = palette_common[3],
                'Clipper' = palette_common[4],
                'PHSF Highbay floor' = palette_common[5],
                'Field Control' = palette_common[6],
                'PHSF Airlock floor' = palette_common[7],
                'Cotton' = palette_common[8],
                'Macrofoam' = palette_common[9],
                'Water' = palette_common[10],
                'NTC' = palette_common[11],
                'mag_beads' = palette_common[12],
                'spin_column' = palette_common[13],
                'cotton' = palette_common[14],
                'macrofoam' = palette_common[15],
                'Polyester Wipe' = palette_common[16],
                'SALSA' = palette_common[17],
                'Mastermix' = palette_common[18],
                'Environmental' = palette_common[19],
                'Device' = palette_common[20],
                'Metal' = palette_common[21],
                'dpcr' = palette_common[22],
                'qpcr' = palette_common[23],
                
                
                
                'NTC'='gray75',
                'unique_features'=palette_common[18],
                'chimeras'=palette_common[19],
                'ns'='gray50',
                'sig'='red',
                
                'panel_outline'='black',
                'general_background'='white',
                'panel_background'='white',
                'plot_outline'='white',
                
                'default'='white')


palette_linetype=c(
  'ATP_RLU'='solid',
  'Fungal_CFU'='solid',
  'Fungal_dPCR'='dotted',
  'Bacterial_CFU'='solid',
  'Bacterial_dPCR'='solid',
  'default'='solid'
)

assign_default(order_vector,21)

palette_shape=c(
  'metal_deposition' = 23,
  'swab_head_retention' = 21,
  'Device' = 23,
  'PHSF Highbay floor' = 22,
  'Field Control' = 24,
  'PHSF Airlock floor' = 21,
  'ATCC' = 23,
  '2.00E+2'=21,
  '2.00E+3'=21,
  '2.00E+4'=21,
  '2.00E+5'=21,
  '2.00E+6'=21,
  'NTC'=21, 
  'magnetic_bead'=21,
  'spin_column'=23,
  'Fungal_dPCR'=21,
  'Bacterial_dPCR'=15,
  
  'ATP_RLU'=17,
  'Fungal_CFU'=15,
  'Fungal_dPCR'=19,#
  'Bacterial_CFU'=15,
  'Bacterial_dPCR'=22,#
  'PMA'=23,
  'No PMA'=21,
  'No PMA' = 21,
  'PMA' = 21,
  'PPR' = 21, 
  'Clipper' = 21,
  'PHSF Highbay floor' = 21,
  'Field Control' = 21,
  'PHSF Airlock floor' = 21,
  'Cotton Wet' = 21,
  'Macrofoam Wet' = 21,
  'Cotton Dry' = 21,
  'Macrofoam Dry' = 21,
  'Cotton' = 21,
  'Macrofoam' = 21,
  'Water' = 21,
  'NTC' = 21,
  'mag_beads' = 21,
  'spin_column' = 21,
  'cotton' = 21,
  'macrofoam' = 21,
  'Polyester Wipe' = 21,
  'SALSA' = 21,
  'Mastermix' = 21,
  'Environmental' = 21,
  'Device' = 21,
  'Metal' = 21,
  'Metal Dry' = 21,
  'Metal Wet' = 21,
  'dpcr' = 21,
  'qpcr' = 21,
  'No PMA' = 21,
  'PMA' = 21,
  'PPR' = 21,
  'PHSF Highbay floor' = 21,
  'Field Control' = 21,
  'PHSF Airlock floor' = 21,
  'Cotton' = 21,
  'Macrofoam' = 21,
  'Water' = 21,
  'NTC' = 21,
  'mag_beads' = 21,
  'spin_column' = 21,
  'cotton' = 21,
  'macrofoam' = 21,
  'Polyester Wipe' = 21,
  'SALSA' = 21,
  'Mastermix' = 21,
  'Environmental' = 21,
  'Device' = 23,
  'Metal' = 21,
  'dpcr' = 21,
  'qpcr' = 21,
  'default'=21
)

assign_default(order_vector,3)

palette_size=c(
  'metal_deposition' = 3,
  'swab_head_retention' = 3,
  '0P'=1,
  '20P'=3,
  '40P'=6,
  '60P'=9,
  '80P'=12,
  
  'less_than_0'=1,
  'log0'=1,
  'log4'=3,
  'log6'=6,
  'log8'=9,
  'log10'=12,
  
  'ng100'=1,
  'g100'=2,
  'g1K'=4,
  'g10K'=6,
  'g100K'=8,
  
  
  
  
  'g0'=1,
  'g10'=3,
  'g20'=6,
  'g40'=9,
  'g80'=12,
  
  'gray1'=3,
  'gray2'=3,
  'gray3'=3,
  'gray4'=3,
  'gray5'=3,
  'PMA'=3,
  'No PMA'=3,
  na.value=3,
  'NTC'=3,
  'magnetic_bead'=3,
  'spin_column'=3,
  'SALSA'=3,
  'Polyester Wipe'=3,
  'Cotton'=3,
  'Macrofoam'=3,
  'Positive Control'=3,
  'Water'=3,
  
  'Cotton Device Control'=3,
  'Macrofoam Device Control'=3,
  'Cotton Environmental Control'=3,
  'Macrofoam Environmental Control'=3,
  'Cotton Dry + Metal'=3,
  'Macrofoam Dry + Metal'=3,
  'Cotton Wet + Metal'=3,
  'Macrofoam Wet + Metal'=3,
  'Water Only'=3,
  
  
  
  
  
  'ATP_RLU'=0,
  'ATP_RLU'=3,
  'Fungal_CFU'=3,
  'Fungal_dPCR'=3,#
  'Bacterial_CFU'=3,
  'Bacterial_dPCR'=3,#
  'PMA'=3,
  'No PMA'=3,
  'No PMA' = 3,
  'PMA' = 3,
  'PPR' = 3,
  'Clipper' = 3,
  'PHSF Highbay floor' = 3,
  'Field Control' = 3,
  'PHSF Airlock floor' = 3,
  'Cotton' = 3,
  'Macrofoam' = 3,
  'Water' = 3,
  'NTC' = 3,
  'mag_beads' = 3,
  'spin_column' = 3,
  'cotton' = 3,
  'macrofoam' = 3,
  'Polyester Wipe' = 3,
  'SALSA' = 3,
  'Mastermix' = 3,
  'Environmental' = 3,
  'Device' = 3,
  'Metal' = 3,
  'Metal Dry' = 3,
  'Metal Wet' = 3,
  'dpcr' = 3,
  'qpcr' = 3,
  'No PMA' = 3,
  'PMA' = 3,
  'PPR' = 3,
  'PHSF Highbay floor' = 3,
  'Field Control' = 3,
  'PHSF Airlock floor' = 3,
  'Cotton' = 3,
  'Macrofoam' = 3,
  'Water' = 3,
  'NTC' = 3,
  'mag_beads' = 3,
  'spin_column' = 3,
  'cotton' = 3,
  'macrofoam' = 3,
  'Polyester Wipe' = 3,
  'SALSA' = 3,
  'Mastermix' = 3,
  'Environmental' = 3,
  'Device' = 3,
  'Metal' = 3,
  'dpcr' = 3,
  'qpcr' = 3,
  'default'=3
)


palette_pattern=c(
  'Bacterial_CFU'='circle',
  'Fungal_CFU'='none',
  
  'ATP_RLU'='none',
  #'Fungal_CFU'='circle',
  'Fungal_CFU'='circle',
  'Fungal_dPCR'='none',
  'Bacterial_CFU'='none',
  'Bacterial_dPCR'='none',
  'default'='none'
)

palette_spacing=c(
  'Fungal_CFU'=.05,
  'Bacterial_CFU'=.1,
  'default'=.05
)
palette_density=c(
  'Fungal_CFU'=.2,
  'Bacterial_CFU'=.4,
  'default'=.2
)
################################################################################
# common Theme and gPlot
# plotting function to add palettes, comoon theme, and ?
# refine more in individual scripts
################################################################################
# library(RColorBrewer)
# brewer_colors <- brewer.pal(8, "Dark2")
# scale_color_manual(values = brewer_colors)



# Define color variables
#background_color <- palette_color['general_background']
background_color <- 'white'
panel_background_color <- "white"
grid_major_color <- "grey90"
grid_minor_color <- "grey95"
strip_background_color <- "grey90"
axis_line_color <- "black"
axis_tick_color <- "black"
legend_key_color <- "white"



################################################################################
theme_common <- 
  #theme(axis.text.x = element_text(angle = 90))+
  theme(
    #  plot.margin = margin(5.5, 5.5, 45, 5.5),
    #panel.background = element_rect(fill = 'slategrey'),
    #  legend.background = element_rect(fill = palette_color['general_background']),
    #legend.position = c(0.5, -.09),
    # plot.title = element_blank(),
    #legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 12,face='bold'),
    axis.title.y = element_markdown(size = 12,face='bold'),
    # axis.title.x = element_text(size = 12,face='bold',margin = margin(t = 10)),
    #  axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    # panel.grid.major.x = element_blank(),
    #  panel.grid.minor.x = element_blank(),
    #theme(strip.text = element_text(size = 14)),
    
    #  legend.background = element_rect(fill = palette_color['general_background']),
    panel.background = element_rect(fill = 'white',color='black'),
    strip.background = element_rect(fill = NA,color=NA),
    legend.position='bottom',
    strip.text = element_text(face = "bold",color='black',size=12),
    # plot.title = element_text(hjust = 0.5),
    plot.title = element_markdown(hjust = 0.5,face='bold'),
    axis.text.x = element_markdown(face='bold'),
    axis.text.y = element_markdown(face='bold')
  )


################################################################################
# common Theme and gPlot
################################################################################
# library(RColorBrewer)
# brewer_colors <- brewer.pal(8, "Dark2")
# scale_color_manual(values = brewer_colors)

theme_elements <- theme(axis.text.x = element_text(angle = 90))+
  theme(
    plot.margin = margin(5.5, 5.5, 45, 5.5),
    plot.background = element_rect(fill = palette_color['general_background']),
    legend.background = element_rect(fill = palette_color['general_background']),
    #legend.position='bottom',
    legend.position = c(0.5, -.09),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(size = 12,face='bold'),
    axis.title.x = element_text(size = 12,face='bold',margin = margin(t = 10)),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    #theme(strip.text = element_text(size = 14)),
    #axis.text.y = element_text(face = "bold"),
    # plot.title = element_text(hjust = 0.5),
    #axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5,face='bold',size=x_text_size),
    caption.title = element_text(hjust = 0.5))
theme_elements <- theme(
  # Text elements
  text = element_text(),            # base text
  title = element_text(),           # base title text
  plot.title = element_text(),      # plot title
  plot.subtitle = element_text(),   # plot subtitle
  plot.caption = element_text(),    # plot caption
  plot.tag = element_text(),        # plot tag
  
  # Axis
  axis.title = element_text(),           # both x and y axis titles
  axis.title.x = element_text(),         # x axis title
  axis.title.y = element_text(),         # y axis title
  axis.text = element_text(),            # both x and y axis text
  axis.text.x = element_text(),          # x axis tick labels
  axis.text.y = element_text(),          # y axis tick labels
  axis.ticks = element_line(),           # tick marks
  axis.ticks.x = element_line(),         # x ticks
  axis.ticks.y = element_line(),         # y ticks
  axis.line = element_line(),            # axis line
  axis.line.x = element_line(),          # x axis line
  axis.line.y = element_line(),          # y axis line
  
  # Legend
  legend.background = element_rect(),    # entire legend area
  legend.key = element_rect(),           # background of legend keys
  legend.title = element_text(),         # legend title
  legend.text = element_text(),          # legend labels
  legend.position = "right",             # "left", "right", "top", "bottom", or c(x, y)
  legend.direction = "vertical",         # or "horizontal"
  legend.justification = "center",       # anchor point inside legend box
  legend.box = NULL,                     # "horizontal" or "vertical"
  legend.spacing = unit(1, "lines"),     # spacing between legends
  
  # Panel (plotting area)
  panel.background = element_rect(),     # panel background
  panel.border = element_rect(),         # border around plotting area
  panel.grid = element_line(),           # major + minor grid lines
  panel.grid.major = element_line(),     # major grid lines
  panel.grid.minor = element_line(),     # minor grid lines
  panel.grid.major.x = element_line(),
  panel.grid.major.y = element_line(),
  panel.grid.minor.x = element_line(),
  panel.grid.minor.y = element_line()
  
  # Strip (for facets)
  ,strip.background = element_rect(),    # background of facet strips
  strip.text = element_text(),           # facet label text
  strip.text.x = element_text(),         # facet label text (x)
  strip.text.y = element_text(),         # facet label text (y)
  
  # Plot background
  plot.background = element_rect()       # entire plot background
)


magnify=1
text_size=10
# Define text size variables

margin_size=10
theme_global <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.tag = element_text(size = rel(2), face = "bold"),
      #  panel.background = element_rect(fill = 'white',color='black'),
      strip.background = element_rect(fill = NA,color=NA),
      #  legend.position='bottom',
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
      plot.background = element_rect(fill = palette_color['general_background'], color = palette_color['general_background'],size=1),
      panel.background = element_rect(fill = palette_color['panel_background'],color = palette_color['panel_outline'],size=1),
      theme(axis.ticks = element_line(size = 4)),
      
      legend.background = element_rect(fill = palette_color['general_background']),
      axis.title = element_markdown(size = rel(1.4), face = 'bold'),
      axis.text = element_markdown(size = rel(1.2), face = 'bold',color='black'),
      legend.title = element_markdown(size = rel(1.2), face = 'bold'),
      legend.text = element_text(size =rel(1.2)),
      legend.key = element_rect(color = "black", size = 0.3),
      #  legend.key.size = unit(.5, "lines"),   
      strip.text = element_text(size = rel(1.2), face = 'bold'),
      plot.title = element_markdown(size = rel(1.6), face = 'bold',hjust=0),
      plot.subtitle = element_text(size = rel(1.2)),
      plot.caption = element_text(size = rel(1.2) , face = 'bold')
    )
}

theme_common<-theme_global(base_size = 11)+
  theme(       
    
  )

common_theme=theme_common
theme_plot<-theme(
  plot.background = element_rect(fill = palette_color['general_background'], color = palette_color['plot_outline'], size = 2)
)



gPlot <- function(p) {
  p=p+
    #   scale_fill_manual(values = palette_color) +
    scale_color_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
    scale_fill_manual(values = palette_color,labels=palette_label,na.translate = FALSE)+
    scale_pattern_manual(values = palette_pattern,labels=palette_label)+
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
  scale_x_discrete(labels = palette_label)+
    theme_common
  print(p)
  return(p)
}
################################################################################
# removes row names,adds meta data
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

