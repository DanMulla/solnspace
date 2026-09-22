# dataviz_reserves.R - figures for the reserve actuator sensitivity check.
#
# Produces the boxplot and parallel-coordinate panels (Fig D in S1 File) for each
# reserve actuator variant sampled by feasiblesolns_reserves.m, across the four
# DOF models for the 50 N downward exertion.
#
#   figures/SCAC5/   SC and AC reserves at 5 Nm
#   figures/SCAC20/  SC and AC reserves at 20 Nm
#   figures/GH01/    GH reserves at 0.1 Nm
#
# Identical plotting code to ../../dataviz.R; only the inputs, output location and
# the set of panels differ. The low effort percentile panels are not produced.

# Libraries
library(conflicted)
library(gdata)
library(gridExtra)
library(tidyverse)

# Clear work space
rm(list = ls())

# Set working directory to current folder
if (requireNamespace('rstudioapi', quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# Reserve actuator variants produced by feasiblesolns_reserves.m
variants <- c('SCAC5', 'SCAC20', 'GH01')
task <- 'task03'
effortpercentile <- 0.05          # unused here; p3 is not produced
dofs <- c('3dof-all', '5dof-all', '8dof-all', '8dof-loweffort', '11dof-all', '11dof-loweffort')

# Muscles, reserves, and groupings
musclenames <- c('Superior Trapezius, Scapula', 'Middle Trapezius, Scapula', 'Inferior Trapezius, Scapula', 'Trapezius, Clavicle', 
                 'Levator Scapula', 'Pectoralis Minor', 'Superior Rhombiod', 'Inferior Rhombiod', 'Inferior Serratus Anterior', 
                 'Middle Serratus Anterior', 'Superior Serratus Anterior', 
                 'Posterior Deltoid', 'Middle Deltoid', 'Anterior Deltoid', 'Coracobrachialis', 'Inferior Infraspinatus', 'Superior Infraspinatus', 'Teres Minor', 
                 'Teres Major', 'Posterior Supraspinatus', 'Anterior Supraspinatus', 'Superior Subscapularis', 'Middle Subscapularis', 
                 'Inferior Subscapularis', 
                 'Biceps Long Head', 'Biceps Short Head', 'Triceps Long Head', 
                 'Superior Latissimus Dorsi', 'Middle Latissimus Dorsi', 'Inferior Latissimus Dorsi', 
                 'Inferior Pectoralis Major', 'Middle Pectoralis Major', 'Superior Pectoralis Major', 
                 'Triceps Medial Head', 'Brachialis', 'Brachioradialis', 'Pronator Teres, Humeral', 
                 'Pronator Teres, Ulnar', 'Supinator', 'Pronator Quadratus', 'Triceps Lateral Head', 'Anconeus')
thoracoscapular <- musclenames[1:11]
glenohumeral <- musclenames[12:24]
thoracohumeral <- musclenames[28:33]
elbow <- musclenames[c(25:27, 34:42)]
costs_colnames <- c("effort1", "effort2", "ghjrf_x", "ghjrf_y", "ghjrf_z", "ghjrf", "sy", "sz")
allreserves <- c('SC_y', 'SC_z', 'SC_x', 'AC_y', 'AC_z', 'AC_x', 'GH_y', 'GH_z', 'GH_yy', 'EL_x', 'PS_y',
                 'SC_y-', 'SC_z-', 'SC_x-', 'AC_y-', 'AC_z-', 'AC_x-', 'GH_y-', 'GH_z-', 'GH_yy-', 'EL_x-', 'PS_y-')
dofs_map <- list(
  "3" = c(7:9, 18:20),
  "5" = c(7:11, 18:22),
  "8" = c(4:11, 15:22),
  "11" = 1:22)
tasktitles <- c('Unloaded', 
                'Down: 20 N', 'Down: 50 N', 'Down: 100 N',
                'Up: 20 N', 'Up: 50 N', 'Up: 100 N', 
                'Left: 20 N', 'Left: 50 N', 'Left: 100 N',
                'Right: 20 N', 'Right: 50 N', 'Right: 100 N')

# Iterating through every task, effort percentile, model and effort landscape
for (variant in variants) {

taskfolder <- paste('outputs/', variant, '/', task, '/solutions/', sep="")
figdir <- paste('figures/', variant, '/', sep="")
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

for (dof in dofs) {
  
  # Set seed for reproducible results (particularly for graphs)
  set.seed(2024)
  
  # File name
  filename <- paste(task, '-samples-', dof, '.txt', sep="")
  filename_costs <- paste(task, '-cost-', dof, '.txt', sep="")
  optsoln_filename <- paste(task, '-mineffort-', str_extract(dof, "[^-]+"), '.txt', sep="")

  # Read in file
  if (file.exists(file.path(taskfolder, filename))) {
    sampledata <- read.delim(paste(taskfolder, filename, sep=""), header=FALSE)
    samplecosts <- read.delim(paste(taskfolder, filename_costs, sep=""), header=FALSE)
    optimaldata <- read.delim(paste(taskfolder, optsoln_filename, sep=""), header=FALSE)
  } else {
    print(paste('Skipping file:', filename))
    next
  }
  
  # Column names for reserves
  ndofs <- sub("([0-9]+)dof.*", "\\1", dof)
  dof_range <- dofs_map[[ndofs]]
  reservenames <- allreserves[dof_range]
  
  # Tidy feasible data (add column names and combine solutions with costs)
  colnames(sampledata) <- c(musclenames, reservenames)
  colnames(samplecosts) <- costs_colnames
  sampledata$Sample <- 1:nrow(sampledata)
  samplecosts$Sample <- 1:nrow(samplecosts)
  sampledata <- merge(sampledata, samplecosts, by="Sample")

  # Tidy optimal data
  optimaldata <- as.data.frame(t(optimaldata))
  colnames(optimaldata) <- c(musclenames, reservenames)
  
  # Finding 'low effort' solutions
  if (grepl("loweffort", dof)) {
    effortthreshold <- quantile(sampledata$effort2, effortpercentile)
    sampledata <- sampledata |> dplyr::filter(effort2 <= effortthreshold)
  }
  
  # Tidy data (pivoting, factoring, muscle groups)
  sampledata_long <- sampledata |> 
    select(-all_of(costs_colnames)) |> 
    pivot_longer(cols=-c(Sample), names_to = 'Muscle', values_to = 'Activation')
  
  optimaldata_long <- optimaldata |> pivot_longer(cols=everything(), names_to = 'Muscle', values_to = 'Activation')
  optimaldata_long$Sample <- -999  # Distinguish from feasible samples
  sampledata_long <- rbind(sampledata_long, optimaldata_long)
  
  # Multiplying activation by 100 (scale from 0 to 100 instead of 0 to 1; easier to interpret (IMO))
  sampledata_long$Activation <- sampledata_long$Activation * 100

  # Factoring the levels (not alphabetized and not the same order as musclenames [which is order of muscles in model])
  sampledata_long$Muscle <- factor(sampledata_long$Muscle, levels=c('Superior Trapezius, Scapula', 'Middle Trapezius, Scapula', 'Inferior Trapezius, Scapula', 'Trapezius, Clavicle', 
                                                                    'Levator Scapula', 'Pectoralis Minor', 'Superior Rhombiod', 'Inferior Rhombiod', 
                                                                    'Superior Serratus Anterior', 'Middle Serratus Anterior', 'Inferior Serratus Anterior', 
                                                                    'Superior Latissimus Dorsi', 'Middle Latissimus Dorsi', 'Inferior Latissimus Dorsi', 
                                                                    'Superior Pectoralis Major', 'Middle Pectoralis Major', 'Inferior Pectoralis Major',
                                                                    'Anterior Deltoid', 'Middle Deltoid', 'Posterior Deltoid', 'Coracobrachialis', 
                                                                    'Superior Infraspinatus', 'Inferior Infraspinatus', 'Teres Minor', 
                                                                    'Anterior Supraspinatus', 'Posterior Supraspinatus',
                                                                    'Superior Subscapularis', 'Middle Subscapularis', 'Inferior Subscapularis', 
                                                                    'Teres Major',
                                                                    'Biceps Long Head', 'Biceps Short Head', 'Triceps Long Head', 'Triceps Lateral Head',
                                                                    'Triceps Medial Head', 'Anconeus', 'Brachialis', 'Brachioradialis',
                                                                    'Pronator Quadratus', 'Pronator Teres, Humeral', 'Pronator Teres, Ulnar', 'Supinator',
                                                                    reservenames))
  
  # Adding in muscle groups
  sampledata_long <- sampledata_long |> mutate(musclegroup = case_when(Muscle %in% thoracoscapular ~ 'thoracoscapular', 
                                                                       Muscle %in% glenohumeral ~ 'glenohumeral',
                                                                       Muscle %in% thoracohumeral ~ 'thoracohumeral',
                                                                       Muscle %in% elbow ~ 'elbow',
                                                                       Muscle %in% allreserves ~ 'reserve'))
  
  # Data Visualization: This is not the cleanest code but works :)
  ## Figure 1 (p1): Box and Whisker plot of the activations
  # Joints each model balances, named in the title for the full solution spaces
  dof_joints <- list('3' = 'GH', '5' = 'GH, Elbow',
                     '8' = 'GH, Elbow, AC', '11' = 'GH, Elbow, AC, SC')
  if (grepl("low", dof)) {
    p1title <- paste(ndofs, ' DOF: Low Effort', sep="")
  } else {
    p1title <- paste(ndofs, ' DOF (', dof_joints[[ndofs]], ')', sep="")
  }
  
  p1 <- ggplot() +
    geom_boxplot(data = sampledata_long |> dplyr::filter(!musclegroup %in% c('reserve') & Sample > 0), 
                 aes(x=Activation, y=fct_rev(Muscle), color=musclegroup, fill=musclegroup),
                 alpha=0.3, outlier.alpha=0.05, outlier.size=1, width=0.5, coef=Inf) +
    labs(x='Activation (% Max)', y='') +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    scale_colour_manual(values = c("thoracoscapular" = "#0072B2", "glenohumeral" = "#009E73", "thoracohumeral" = "#D55E00", "elbow" = "#EE3377", "reserve" = "#BBBBBB")) +
    scale_fill_manual(values = c("thoracoscapular" = "#0072B2", "glenohumeral" = "#009E73", "thoracohumeral" = "#D55E00", "elbow" = "#EE3377", "reserve" = "#BBBBBB")) +
    ggtitle(p1title) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y=element_text(size=12, color='black'),
          axis.text.x=element_text(size=14, color='black'),
          axis.title=element_text(size=16, face="bold"),
          plot.title=element_text(size=16, face="bold"),
          legend.position="none")
  
  # Save figure (only for effort percentile is 0.05 -- avoid repetition)
  p1_filename = paste('Activations-Distributions-', dof, '.png', sep="")
  if (effortpercentile==0.05) {ggsave(paste(figdir, p1_filename, sep=""), p1, height=8, width=5, dpi=900)}
  
  
  ## Figure #2 (p2): Parallel coordinates plot of random 100 samples
  selectsamples <- sample(sampledata$Sample, 100, replace=FALSE)
  p2 <- ggplot(data = sampledata_long |> dplyr::filter(!musclegroup %in% c('reserve') & Sample %in% c(selectsamples))) +
    geom_line(aes(x=fct_rev(Muscle), y=Activation, group=Sample), colour='grey50', alpha=0.1) +
    labs(y='Activation (% Max)', x='') +
    coord_flip() +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    ggtitle(p1title) +
    theme_bw() +
    theme(strip.background = element_rect(fill="grey90", linetype = 0), 
          strip.text = element_text(colour = 'black', size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.5, "lines"), 
          axis.text.y=element_text(size=12, color='black'),
          axis.text.x=element_text(size=14, color='black'),
          axis.title=element_text(size=16, face="bold"),
          plot.title=element_text(size=16, face="bold"),
          legend.title=element_blank(),
          legend.position="none")
  
  # Save figure (only for effort percentile is 0.05 -- avoid repetition)
  p2_filename = paste('Activations-Samples-', dof, '.png', sep="")
  if (effortpercentile==0.05) {ggsave(paste(figdir, p2_filename, sep=""), p2, height=8, width=5, dpi=900)}
  
  # Creating alternative versions with 8 dof model with different titles (p1_tasktitle, p2_tasktitle) or with optimal data plotted (p2_optimal_all) - for publication
  if (dof=='8dof-all') {
    tasknum <- as.numeric(regmatches(task, gregexpr("[0-9]+", task)))
    p1_tasktitle <- p1 + ggtitle(tasktitles[tasknum])
    p2_tasktitle <- p2 + ggtitle(tasktitles[tasknum])
    
    p2_optimal_all <- p2 + 
      geom_line(data = sampledata_long |> dplyr::filter(!musclegroup %in% c('reserve') & Sample==-999),
      aes(x=fct_rev(Muscle), y=Activation, group=Sample), color='black', linewidth=1) +
      ggtitle('All solutions')
    
    p1_tasktitle_filename <- paste('Activations-Distributions-', dof, '-titledtask.png', sep="")
    p2_tasktitle_filename <- paste('Activations-Samples-', dof, '-titledtask.png', sep="")
    p2_optimal_all_filename <- paste('Activations-Samples-', dof, '-titledtask_optimal.png', sep="")
    
    if (effortpercentile == 0.05) {
      ggsave(paste(figdir, p1_tasktitle_filename, sep=""), p1_tasktitle, height=8, width=5, dpi=900)
      ggsave(paste(figdir, p2_tasktitle_filename, sep=""), p2_tasktitle, height=8, width=5, dpi=900)
      ggsave(paste(figdir, p2_optimal_all_filename, sep=""), p2_optimal_all, height=8, width=5, dpi=900)
    }
  }

  # Figure #3 (p3): Creating alternative versions with 8 dof with optimal data plotted (p2_optimal_all)
  if (dof=='8dof-loweffort') {
    
    # Adding in plot title, colours, alpha based on percentile
    tasknum <- as.numeric(regmatches(task, gregexpr("[0-9]+", task)))
    if (effortpercentile == 0.01) {
      superscript <- 'st'
      p3colour <- '#3f007d'
      p3alpha <- 0.2
    } else if (effortpercentile == 0.05) {
      superscript <- 'th'
      p3colour <- '#6a51a3'
      p3alpha <- 0.2
    } else if (effortpercentile == 0.10) {
      superscript <- 'th'
      p3colour <- '#9e9ac8'
      p3alpha <- 0.2
    }
    effortpercent <- effortpercentile * 100
    p3_title <- paste(effortpercent, superscript, ' %ile', sep="")
    
    p3 <- ggplot(data = sampledata_long |> dplyr::filter(!musclegroup %in% c('reserve') & Sample %in% c(selectsamples))) +
      geom_line(aes(x=fct_rev(Muscle), y=Activation, group=Sample), colour=p3colour, alpha=p3alpha) +
      geom_line(data = sampledata_long |> dplyr::filter(!musclegroup %in% c('reserve') & Sample==-999),
                aes(x=fct_rev(Muscle), y=Activation, group=Sample), color='black', linewidth=1) +
      labs(y='Activation (% Max)', x='') +
      coord_flip() +
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
      ggtitle(p3_title) +
      theme_bw() +
      theme(strip.background = element_rect(fill="grey90", linetype = 0), 
            strip.text = element_text(colour = 'black', size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(1.5, "lines"), 
            axis.text.y=element_text(size=12, color='black'),
            axis.text.x=element_text(size=14, color='black'),
            axis.title=element_text(size=16, face="bold"),
            plot.title=element_text(size=19, face="bold"),
            legend.title=element_blank(),
            legend.position="none")
    
    # Saving plot  
    p3_filename <- paste('Activations-Samples-', dof, '-optimalvslow', effortpercent, '.png', sep="")
    ggsave(paste(figdir, p3_filename, sep=""), p3, height=8, width=5, dpi=900)
  }
  
  # Clean work space
  gdata::keep(variants, variant, figdir, task, effortpercentile, dofs, taskfolder, musclenames, thoracoscapular, glenohumeral, thoracohumeral, elbow,
              costs_colnames, allreserves, dofs_map, dof, tasktitles, sure = TRUE)
}

}  # reserve variant
