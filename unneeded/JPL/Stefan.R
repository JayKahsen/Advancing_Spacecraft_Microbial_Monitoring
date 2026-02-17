source(paste0(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)),'/_globalStuff.R'))
###############################################################################
# Description Start line 80ish
################################################################################
script_title='Stefan'
options(scipen = 999) # don;t use scientific notation
################################################################################
# Running info and switches
################################################################################

# output_data=paste0('output_data/',script_title,'/')
# output_plot=paste0('output_plot/',script_title,'/')

################################################################################
# Extra Script Variables and functions
################################################################################
test_alpha=.5
spacer=.8
################################################################################
# Loading Data
################################################################################
qLoad('data_tables/mapping_data.csv')

qLoad('data_tables/reference_weights.csv')


df_reference_weights=reference_weights%>%
  mutate(reference = gsub("-", ".", reference))%>%
  mutate(reference = gsub(" ", ".", reference))%>%
  ungroup()
  

df_Zymo=df_reference_weights%>%
  filter(reference!='Lambda_NEB')%>%
  mutate(standard=case_when(
    reference %in% c('Saccharomyces_cerevisiae_strain_Y.567_JAQZRZ010000000','Cryptococcus_neoformans_strain_Y.2534_JAQZRY010000001')~2,
    TRUE~12
  ))%>%
  mutate(is_Zymo=TRUE)%>%

  mutate(reference_original=reference)%>%
  mutate(reference=case_when(
    reference_original=='Lambda_NEB'~'Lambda_NEB',
    reference_original=='Bacillus_B354_CP118021_CP118022'~'Bacillus',
    reference_original=='Escherichia_coli_B.1109_CP117971_CP117972'~'Escherichia_coli',
    reference_original=='Pseudomonas_aeruginosa_B.3509_CP117974.CP117975'~'Pseudomonas_aeruginosa',
    reference_original=='Salmonella_enterica_B.4212_CP117976_CP117977_CP117978'~'Salmonella_enterica',
    reference_original=='Staphylococcus_aureus_B.41012_CP117979'~'Staphylococcus_aureus',
    reference_original=='Saccharomyces_cerevisiae_strain_Y.567_JAQZRZ010000000'~'Saccharomyces_cerevisiae',
    reference_original=='Cryptococcus_neoformans_strain_Y.2534_JAQZRY010000001'~'Cryptococcus_neoformans',
    reference_original=='Listeria.monocytogenes.B.33116'~'Listeria.monocytogenes',
    reference_original=='Lactobacillus.fermentum.B.1840'~'Lactobacillus.fermentum',
    reference_original=='Enterococcus.faecalis.B.537'~'Enterococcus.faecalis',
    TRUE~'NA'
  ))%>%
  arrange(reference.bp)%>%
  select(-reference.bp)%>%
  ungroup()

reference_Zymo=df_Zymo$reference


#reference_Zymo_order=c( 'Lactobacillus.fermentum','Staphylococcus_aureus','Enterococcus.faecalis','Listeria.monocytogenes','Bacillus','Salmonella_enterica','Escherichia_coli','Pseudomonas_aeruginosa','Saccharomyces_cerevisiae','Cryptococcus_neoformans' )

qnvPrint(df_reference_weights$reference)

################################################################################
# Data Manipulations
################################################################################
df1=mapping_data%>%
  pivot_longer(
    cols = -c(sample_id, lambda, zymo, score),
    names_to = "reference",
    values_to = "value"
  )%>%
  pivot_wider(
    names_from = score,
    values_from = value
  )%>%
  mutate(LoZ=as.character(lambda/zymo))%>%
  mutate(LoZ_label=case_when(
    LoZ=='Inf'~'0Z:1L',
    LoZ=='0'~'1Z:0L',
    LoZ=='1'~'1Z:1L',
    LoZ=='3'~'1Z:3L',
    LoZ==NaN~'control',
    TRUE~LoZ
    ))%>%
  mutate(reference_original=reference)%>%
  mutate(reference=case_when(
    reference_original=='Lambda_NEB'~'Lambda_NEB',
    reference_original=='Bacillus_B354_CP118021_CP118022'~'Bacillus',
    reference_original=='Escherichia_coli_B.1109_CP117971_CP117972'~'Escherichia_coli',
    reference_original=='Pseudomonas_aeruginosa_B.3509_CP117974.CP117975'~'Pseudomonas_aeruginosa',
    reference_original=='Salmonella_enterica_B.4212_CP117976_CP117977_CP117978'~'Salmonella_enterica',
    reference_original=='Staphylococcus_aureus_B.41012_CP117979'~'Staphylococcus_aureus',
    reference_original=='Saccharomyces_cerevisiae_strain_Y.567_JAQZRZ010000000'~'Saccharomyces_cerevisiae',
    reference_original=='Cryptococcus_neoformans_strain_Y.2534_JAQZRY010000001'~'Cryptococcus_neoformans',
    reference_original=='Listeria.monocytogenes.B.33116'~'Listeria.monocytogenes',
    reference_original=='Lactobacillus.fermentum.B.1840'~'Lactobacillus.fermentum',
    reference_original=='Enterococcus.faecalis.B.537'~'Enterococcus.faecalis',
    TRUE~'NA'
  ))%>%
  mutate(picograms = sub(".*-(.*)-.*", "\\1", sample_id))%>%
  mutate(picogram_label=case_when(
    picograms=='05'~'0.5 pg',
    picograms=='5'~'5.0 pg',
    picograms=='50'~'50.0 pg',
    TRUE~'NTC'
  ))%>%
  
  mutate(short_name=paste('s',LoZ_label,picograms,sep='_'))%>%
  group_by(short_name,reference)%>%
  mutate(replicate=row_number())%>%
  ungroup()%>%
  mutate(sample_name=paste('s',LoZ_label,picograms,replicate,sep='_'))%>%

  
ungroup()

reference_order=c(reference_Zymo,setdiff(unique(df1$reference),reference_Zymo))
reference_order
LoZ_label_order=unique(df1$LoZ_label)
LoZ_label_order
picogram_label_order=unique(df1$picogram_label)
picogram_label_order

df=df1%>% # short_name, #replicate # meaningful sample_name
  left_join(df_reference_weights)%>%
  left_join(df_Zymo)%>%
  mutate(is_Zymo = if_else(is.na(is_Zymo), FALSE, is_Zymo))%>%
  group_by(sample_name)%>%
  mutate(
    weighted_Zymo = sum(weighted[is_Zymo], na.rm = TRUE),  # assuming is_Zymo marks Zymo features
    weighted_detection_percent = ifelse(is_Zymo, weighted / weighted_Zymo * 100, NA),
    weighted_detection_diff = weighted_detection_percent - standard,
    Ideal_Score = sum(abs(weighted_detection_diff), na.rm = TRUE)
  ) %>%
  group_by(sample_name)%>%
  mutate(reads_total=sum(raw))%>%
  mutate(reads_zymo = sum(raw[is_Zymo]))%>%
  mutate(reads_lambda = sum(raw[!is_Zymo]))%>%
  mutate(reads_check=reads_total-reads_zymo-reads_lambda)%>%
  mutate(reads_ZoL=reads_zymo/reads_lambda)%>%
  mutate(reference=factor(reference,levels = reference_order))%>%
  mutate(LoZ_label=factor(LoZ_label,levels = LoZ_label_order))%>%
  mutate(picogram_label=factor(picogram_label,levels = picogram_label_order))%>%
  ungroup()


names(df)

unique(df$reference)

df_reference_weights$reference
################################################################################
plot_title <-'Ideal Scores'
################################################################################
# df2 <- df1 %>%
#   filter(reference %in% reference_Zymo)%>%
#   filter(zymo > 0)%>%
#   #mutate(picogram_label = as.numeric(picogram_label))%>%
#   ungroup()

df_plot=df%>%
  filter(zymo > 0)%>%
  select(sample_name,picogram_label,LoZ_label,Ideal_Score)%>%
  mutate(title=plot_title)%>%
  distinct()


df_summary <- df_plot %>%
  group_by(title,picogram_label, LoZ_label) %>%
  summarise(
    mean_Ideal_Score = mean(Ideal_Score, na.rm = TRUE),
    sd_Ideal_Score = sd(Ideal_Score, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    ymin = mean_Ideal_Score - sd_Ideal_Score,
    ymax = mean_Ideal_Score + sd_Ideal_Score
  )

library(broom)
library(stringr)

stat_results <- df_plot %>%
  group_by(picogram_label) %>%
  group_modify(~ {
    model <- aov(Ideal_Score ~ LoZ_label, data = .x)
    tukey <- TukeyHSD(model)
    broom::tidy(tukey) %>% 
      filter(adj.p.value < 0.05)  # no mutate needed
  }) %>%
  ungroup()

df_y_position=df_plot%>%
  group_by(picogram_label,LoZ_label)%>%
  mutate(y_position=max(Ideal_Score))

LoZ_levels <- levels(df_plot$LoZ_label)



stat_results_formatted <- stat_results %>%
  filter(adj.p.value < 0.05) %>%
  mutate(
    group1 = str_split(contrast, "-", simplify = TRUE)[,1],
    group2 = str_split(contrast, "-", simplify = TRUE)[,2]
  ) %>%
  rowwise() %>%
  mutate(
    base_y = max(df_plot$Ideal_Score[
      df_plot$picogram_label == picogram_label &
        df_plot$LoZ_label %in% c(group1, group2)
    ], na.rm = TRUE) + spacer/2,
    distance = abs(match(group1, LoZ_levels) - match(group2, LoZ_levels)),
    label = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01 ~ "**",
      adj.p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup()


# Stack bars higher if they're wide (distance ≥ 2) and share the same facet
stat_results_stacked <- stat_results_formatted %>%
  group_by(picogram_label) %>%
  arrange(distance, base_y) %>%
  mutate(
    y.position = base_y + spacer * (cumsum(duplicated(base_y) & distance >= 2))
  ) %>%
  ungroup()

stat_results_formatted 
stat_results_stacked
################################################################################



library(ggpubr)


p=ggplot(df_plot, aes(x = LoZ_label,  color = LoZ_label)) +
  geom_errorbar(data = df_summary, aes(x = LoZ_label, ymin = ymin, ymax = ymax), width = 0.1,linewidth=1.2,alpha=test_alpha) +
  geom_point(aes(y = Ideal_Score,),position = position_jitter(width = 0.2),size=3,alpha=test_alpha) +
  facet_grid(title~picogram_label, scales = "free") +
stat_pvalue_manual(
    stat_results_stacked,
    label = "label",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.015,size=5
  )+
  labs(
  x = "Zymo : Lambda Ratio",
   y = "Ideal Score"
  )+
  theme(  axis.text.x = element_text(size = axis_text_size,face = 'bold'))+
  ggtitle(plot_title)

plot_Ideal=gPlot(p)

ggsave(paste(plot_title,'.png'),height = 8,width=12)
################################################################################
plot_title <-'Zymo Detection Percentage'
################################################################################
df_plot=df%>%
  filter(reference %in% reference_Zymo)%>%
  filter(zymo > 0)%>%
  mutate(width=weighted_detection_percent/standard)%>%
  select(sample_name,reference,picogram_label,LoZ_label,weighted_detection_percent,standard)%>%
  mutate(title=plot_title)%>%
ungroup()

df_summary <- df_plot %>%
  mutate(width=weighted_detection_percent/standard)%>%
  group_by(title,picogram_label, LoZ_label,reference) %>%
  summarise(
    mean_detection = mean(weighted_detection_percent, na.rm = TRUE),
    sd_detection = sd(weighted_detection_percent, na.rm = TRUE),
    mean_width = mean(width, na.rm = TRUE),
    sd_width = sd(width, na.rm = TRUE),
    standard=mean(standard),
    .groups = 'drop'
  ) %>%
  mutate(
    ymin = mean_detection - sd_detection,
    ymax = mean_detection + sd_detection,
    xmin = mean_width - sd_width,
    xmax = mean_width + sd_width
  )%>%
  mutate(standard_color='standard')%>%
  ungroup()

################################################################################
dodge_width=.8
p=ggplot(df_summary) +
  facet_grid(title~picogram_label) +
  
  # Main bars
  
  geom_col(
    aes(x = standard, y = reference, color = standard_color, fill = standard_color),
    position = position_dodge(width = dodge_width),width=.9,alpha=test_alpha
  ) +
  
  geom_col(
    aes(x = mean_detection, y = reference, fill = LoZ_label),
    position = position_dodge(width = dodge_width),width=.7,alpha=test_alpha
  ) +

  # Horizontal error bars
  geom_errorbarh(
    aes(xmin = ymin, xmax = ymax, y = reference, color = LoZ_label),
    position = position_dodge(width = dodge_width),  # Align with dodge in bars
    height = 0.4,linewidth=1
  ) +
  
  # Log-scaled x-axis
  #scale_x_continuous(trans = "log2", labels = log_labels2) +
  coord_cartesian(xlim = c(0.125, NA)) +
  labs(
    x = "Percent Detection",
    y = ""
  )+
  theme(  axis.text.x = element_text(size = axis_text_size,face = 'bold'))+
  ggtitle(plot_title)
p
plot_ZymoB=gPlot(p)
ggsave(paste(plot_title,'.png'),height = 8,width=12)
################################################################################
plot_title='Relative Zymo Detection'
df_plot=df_plot%>%
  mutate(title=plot_title)
  
  df_summary=df_summary%>%
    mutate(title=plot_title)
  
dodge_width=.8
p=ggplot(df_summary) +
  facet_grid(title~picogram_label) +
  
  # Main bars
  geom_col(
    aes(x = mean_width, y = reference, fill = LoZ_label),
    position = position_dodge(width = dodge_width),width=.7,alpha=test_alpha
  ) +
  
  # Horizontal error bars
  geom_errorbarh(
    aes(xmin = xmin, xmax = xmax, y = reference, color = LoZ_label),
    position = position_dodge(width = dodge_width),  # Align with dodge in bars
    height = 0.4,linewidth=1
  ) +
  
  # Log-scaled x-axis
  scale_x_continuous(trans = "log2", labels = log2_labels_bold) +
  coord_cartesian(xlim = c(0.125, NA)) +
  
  # Vertical reference line
  geom_vline(aes(xintercept = 1), color = 'magenta', linetype = 'dashed') +
  labs(
    x = "Log2(Percent Detection / Standard)",
    y = ""
  )+
  theme(axis.text.x = element_markdown(size = axis_text_size,face = 'bold'))+
  ggtitle(plot_title)
p
plot_ZymoM=gPlot(p)

ggsave(paste(plot_title,'.png'),height = 8,width=12)
################################################################################
plot_title <-'Read Depth'
################################################################################

df_plot=df%>%
  filter(zymo > 0)%>%
  select(sample_name,reads_total,reads_zymo,reads_lambda,LoZ_label,picogram_label)%>%
  distinct()%>%
  pivot_longer(cols = c(reads_lambda,reads_zymo),names_to = 'read_type',values_to = 'reads')%>%
  mutate(title=plot_title)%>%
  ungroup()

read_type_order=c('reads_zymo','reads_lambda')

df_summary <- df_plot %>%
  group_by(title,picogram_label, LoZ_label,read_type) %>%
  summarise(
    mean_reads = mean(reads, na.rm = TRUE),
    sd_reads = sd(reads, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    ymin = mean_reads - sd_reads,
    ymax = mean_reads + sd_reads
  )%>%
  mutate(read_type=factor(read_type,levels = read_type_order))%>%
  ungroup()
  

dodge_width=.8

p=ggplot(df_summary) +
  facet_grid(title~picogram_label) +
  
  # Main bars
  geom_col(
    aes(x = LoZ_label, y = mean_reads, fill = read_type),
    position = position_dodge(width = dodge_width),width=.7,alpha=test_alpha
  )+
  geom_errorbar(data = df_summary,
    aes(x=LoZ_label,ymin = ymin, ymax = ymax, color = read_type),
    position = position_dodge(width = dodge_width),  # Align with dodge in bars
    width = 0.2,linewidth=1.2
  ) +
  labs(
    x = "Zymo : Lambda Ratio",
    y = "Total Reads"
  )+
  theme(  axis.text.x = element_text(size = axis_text_size,face = 'bold'))+
  ggtitle(plot_title)
p
plot_reads=gPlot(p)

  plot_reads
  ggsave(paste(plot_title,'.png'),height = 8,width=12)
################################################################################
# Statistics
################################################################################
  
  
################################################################################
# Assembly
################################################################################

################################################################################
plot_layout=plot_reads/plot_ZymoB/plot_ZymoM/plot_Ideal
plot_layout

ggsave('Ideal.png',plot_layout,height = 14,width=24)
################################################################################
stop()
plot_title='Read Depth'
plot_title='ANOVA Ideal Score?'
plot_title='by amount Lambda and Zymo' # show for effect of Zymo and Lambda concentrations
plot_title='Lambda/Z ratio?'
  



################################################################################
# Extra Script Variables and functions
################################################################################