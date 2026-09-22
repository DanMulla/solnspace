# Libraries
library(conflicted)
library(gdata)
library(gridExtra)
library(tidyverse)
library(viridis)

# Clear work space
rm(list = ls())

# Set working directory to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Output folder for the summary statistics written below
dir.create("statistics", showWarnings = FALSE)

# Setting up data frame
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
reservenames <- c('SC_y', 'SC_z', 'SC_x', 'AC_y', 'AC_z', 'AC_x', 'GH_y', 'GH_z', 'GH_yy', 'EL_x', 'PS_y', 
                  'SC_y-', 'SC_z-', 'SC_x-', 'AC_y-', 'AC_z-', 'AC_x-', 'GH_y-', 'GH_z-', 'GH_yy-', 'EL_x-', 'PS_y-')
dofs_map <- list(
  "3" = c(7:9, 18:20),
  "5" = c(7:11, 18:22),
  "8" = c(4:11, 15:22),
  "11" = 1:22)
costs_colnames <- c("effort1", "effort2", "ghjrf_x", "ghjrf_y", "ghjrf_z", "ghjrf", "sy", "sz")
feasibledata_cols <- c(musclenames, reservenames, 'Sample', 'Task', 'DOF', costs_colnames)
feasibledata <- data.frame(matrix(ncol = length(feasibledata_cols), nrow = 0))
loweffortdata <- data.frame(matrix(ncol = length(feasibledata_cols), nrow=0))
optimaldata <- data.frame(matrix(ncol = length(feasibledata_cols), nrow=0))
colnames(feasibledata) <- feasibledata_cols
colnames(loweffortdata) <- feasibledata_cols
colnames(optimaldata) <- feasibledata_cols
dofs <- c('3dof', '5dof', '8dof', '11dof')

# Grab all files
taskfolders <- list.dirs('outputs/.', full.names = FALSE, recursive = FALSE)

for (task in 1:length(taskfolders)) {
  
  taskname <- taskfolders[task]
  taskpath <- file.path('outputs', taskname, 'solutions')
  for (dof in dofs) {
    filename <- paste(taskname, '-samples-', dof, '-all.txt', sep="")
    costfile <- paste(taskname, '-cost-', dof, '-all.txt', sep="")
    
    loweffort_samplesfile <- paste(taskname, '-samples-', dof, '-loweffort.txt', sep="")
    loweffort_costsfile <- paste(taskname, '-cost-', dof, '-loweffort.txt', sep="")
    optimalfile <- paste(taskname, '-mineffort-', dof, '.txt', sep="")
    
    # Column names for reserves
    ndofs <- sub("([0-9]+)dof.*", "\\1", dof)
    dof_range <- dofs_map[[ndofs]]
    reserve_temp <- reservenames[dof_range]
    
    if (file.exists(file.path(taskpath, filename))) {
      
      # Grab data (all feasible samples, samples across effort landscape, and optimal data)
      filedata <- read.delim(file.path(taskpath, filename), header=FALSE)
      costdata <- read.delim(file.path(taskpath, costfile), header=FALSE)
      loweffort_samples <- read.delim(file.path(taskpath, loweffort_samplesfile), header=FALSE)
      loweffort_costs <- read.delim(file.path(taskpath, loweffort_costsfile), header=FALSE)
      optimaldata_temp <- read.delim(file.path(taskpath, optimalfile), header=FALSE)
      
      # Tidy Feasible Data
      # Add coln names, add sample #, task, model dof, fill in not used reserves as NA, and then bind to data from other tasks
      colnames(filedata) <- c(musclenames, reserve_temp)
      filedata$Sample <- 1:nrow(filedata)
      filedata$Task <- taskname
      filedata$DOF <- dof
      colnames(costdata) <- costs_colnames
      costdata$Sample <- 1:nrow(costdata)
      feasibledata_temp <- merge(filedata, costdata, by="Sample")
      feasibledata_temp[setdiff(names(feasibledata), names(feasibledata_temp))] <- NA  # Filling in not used reserves as NA
      feasibledata <- rbind(feasibledata, feasibledata_temp)
      
      # Tidy Low Effort Data (this is feasible data sampled across the effort landscape to be used to subset to low effort data)
      colnames(loweffort_samples) <- c(musclenames, reserve_temp)
      loweffort_samples$Sample <- 1:nrow(loweffort_samples)
      loweffort_samples$Task <- taskname
      loweffort_samples$DOF <- dof
      colnames(loweffort_costs) <- costs_colnames
      loweffort_costs$Sample <- 1:nrow(loweffort_costs)
      loweffortdata_temp <- merge(loweffort_samples, loweffort_costs, by="Sample")
      loweffortdata_temp[setdiff(names(loweffortdata), names(loweffortdata_temp))] <- NA  # Filling in not used reserves as NA
      loweffortdata <- rbind(loweffortdata, loweffortdata_temp)
      
      # Tidy optimal data
      optimaldata_temp <- as.data.frame(t(optimaldata_temp))
      colnames(optimaldata_temp) <- c(musclenames, reserve_temp)
      optimaldata_temp$Sample <- -999
      optimaldata_temp$Task <- taskname
      optimaldata_temp$DOF <- dof
      optimaldata_temp[setdiff(names(optimaldata), names(optimaldata_temp))] <- NA  # Filling in not used reserves as NA
      optimaldata<- rbind(optimaldata, optimaldata_temp)
    }
  }
}

# Calculate effort costs for the optimal solutions
sum_of_squares <- function(x) { sum(x^2, na.rm=TRUE)}
optimaldata$effort1 <- rowSums(optimaldata[,c(musclenames, reservenames)], na.rm=TRUE)
optimaldata$effort2 <- apply(optimaldata[,c(musclenames, reservenames)], 1, sum_of_squares)
optimaldata <- optimaldata |> mutate(effort1 = case_when(DOF == '3dof' ~ effort1 / 48,
                                                         DOF == '5dof' ~ effort1 / 52,
                                                         DOF == '8dof' ~ effort1 / 58,
                                                         DOF == '11dof' ~ effort1 / 64))
optimaldata <- optimaldata |> mutate(effort2 = case_when(DOF == '3dof' ~ effort2 / 48,
                                                         DOF == '5dof' ~ effort2 / 52,
                                                         DOF == '8dof' ~ effort2 / 58,
                                                         DOF == '11dof' ~ effort2 / 64))

# Tidy low effort and optimal data to remove unnecessary columns and pivot longer so muscles aren't separate columns
loweffortdata_long <- loweffortdata |> select(-c(all_of(reservenames), effort1, ghjrf_x, ghjrf_y, ghjrf_z, ghjrf, sy, sz)) |> pivot_longer(cols=-c(Sample, Task, DOF, effort2), names_to = 'Muscle', values_to = 'Activation')
optimaldata_long <- optimaldata |> select(-c(all_of(reservenames), effort1, ghjrf_x, ghjrf_y, ghjrf_z, ghjrf, sy, sz)) |> pivot_longer(cols=-c(Sample, Task, DOF, effort2), names_to = 'Muscle', values_to = 'Activation')

# Multiplying activation by 100 (scale from 0 to 100 instead of 0 to 1; easier to interpret (IMO))
loweffortdata_long$Activation <- loweffortdata_long$Activation * 100
optimaldata_long$Activation <- optimaldata_long$Activation * 100

# Combining the low effort and optimal data
loweffortoptimaldata <- rbind(loweffortdata_long |> dplyr::filter(DOF %in% c('8dof')), optimaldata_long |> dplyr::filter(DOF %in% c('8dof')))

# Calculating differences between sampled solutions and optimal solution (only using the 8 dof models)
# Summary statistic #1: Activation Difference (to look at each individual muscle)
# Summary statistic #2: Activation RMSD (to compile differences from each solution)
# Summary statistic #3: Effort percentile (to determine where solutions are in the effort landscape; optimal = 0 %ile)
loweffortoptimaldata <- loweffortoptimaldata |> group_by(Task, Muscle) |> mutate(Difference = Activation - Activation[Sample == -999],
                                                                                 effort2percentile = percent_rank(effort2) * 100,
                                                                                 effortdif = effort2 - effort2[Sample == -999])
loweffortoptimaldata_RMSD <- loweffortoptimaldata |> group_by(Task, Sample, effort2, effort2percentile) |> summarize(RMSD = sqrt(mean(Difference^2)))

# Running through each task, generate the following plots:
# Fig #1 (p1): Effort Percentile vs. RMSD
# Fig #2 (p2): Histogram of RMSD (different colours for different percentiles)
# Fig #3 (p3): Histogram of Activation Differences (different colours for different percentiles)
for (taskcondition in unique(loweffortoptimaldata_RMSD$Task)) {
  
  p1 <- ggplot(data = loweffortoptimaldata_RMSD |> dplyr::filter(Task %in% taskcondition)) + 
    geom_point(aes(x=effort2percentile, y=RMSD, colour=effort2)) +
    geom_point(data = loweffortoptimaldata_RMSD |> dplyr::filter(Task %in% taskcondition & Sample == -999),
               aes(x=effort2percentile, y=RMSD), shape=1, stroke=1, colour='black') +
    scale_colour_viridis(option='cividis', direction=-1) +
    xlab('Effort Percentile') +
    ylab('RMSD (% Max)') +
    labs(colour = 'Effort') +
    theme_bw() +
    theme(strip.background = element_rect(fill="grey90", linetype = 0), 
          strip.text = element_text(colour = 'black', size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.5, "lines"),
          axis.text.x=element_text(size=14, color='black'),
          axis.text.y=element_text(size=14, color='black'),
          axis.title=element_text(size=20, face="bold"),
          legend.text=element_text(size=14, color='black'),
          legend.title=element_text(size=18, face='bold')) +
    annotate("text", x=2, y=0, label='Minimum Effort Solution', color="black", size=6, fontface="bold", hjust=0)
  
  p2 <- ggplot() +
    geom_histogram(data = loweffortoptimaldata_RMSD |> dplyr::filter(Task %in% taskcondition & Sample != -999 & effort2percentile < 10), 
                   aes(x=RMSD), fill='#9e9ac8', color='white') + 
    geom_histogram(data = loweffortoptimaldata_RMSD |> dplyr::filter(Task %in% taskcondition & Sample != -999 & effort2percentile < 5), 
                   aes(x=RMSD), fill='#6a51a3', color='white') + 
    geom_histogram(data = loweffortoptimaldata_RMSD |> dplyr::filter(Task %in% taskcondition & Sample != -999 & effort2percentile < 1), 
                   aes(x=RMSD), fill='#3f007d', color='white') +
    theme_bw() +
    scale_x_continuous(limits=c(0, 15), breaks=seq(0, 15, 5)) +
    labs(x='RMSD (% Max)', y='Count') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=14, color='black'),
          axis.title = element_text(size=20, face='bold'))
  
  p3 <- ggplot() +
    geom_histogram(data = loweffortoptimaldata |> dplyr::filter(Task %in% taskcondition & Sample != -999 & effort2percentile < 10), 
                   aes(x=abs(Difference)), fill='#9e9ac8', color='white') + 
    geom_histogram(data = loweffortoptimaldata |> dplyr::filter(Task %in% taskcondition & Sample != -999 & effort2percentile < 5), 
                   aes(x=abs(Difference)), fill='#6a51a3', color='white') + 
    geom_histogram(data = loweffortoptimaldata |> dplyr::filter(Task %in% taskcondition & Sample != -999 & effort2percentile < 1), 
                   aes(x=abs(Difference)), fill='#3f007d', color='white') +
    theme_bw() +
    scale_x_continuous(limits=c(0, 30), breaks=seq(0, 30, 10)) +
    labs(x='Absolute Difference (% Max)', y='Count') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=14, color='black'),
          axis.title = element_text(size=20, face='bold'))
  
  # Plot titles
  p1_filename <- paste('outputs/', taskcondition, '/solutions/Fig-RMSDvsEffortPercentiles.png', sep="")
  p2_filename <- paste('outputs/', taskcondition, '/solutions/Fig-RMSDHistogram.png', sep="")
  p3_filename <- paste('outputs/', taskcondition, '/solutions/Fig-DifferenceHistogram.png', sep="")
  
  # Save plots
  ggsave(p1_filename, p1, height=4, width=6, dpi=900)
  ggsave(p2_filename, p2, height=2, width=6, dpi=900)
  ggsave(p3_filename, p3, height=2, width=6, dpi=900)
}

# Summary statistics on RMSD (optimal vs. low effort): 1st, 5th, 10th percentiles (for each task)
RMSD_1 <- loweffortoptimaldata_RMSD |> dplyr::filter(Sample != -999 & effort2percentile < 1) |> group_by(Task) |> summarize(Mean_RMSD=mean(RMSD), SD_RMSD=sd(RMSD), Min_RMSD=min(RMSD), Max_RMSD=max(RMSD))
RMSD_5 <- loweffortoptimaldata_RMSD |> dplyr::filter(Sample != -999 & effort2percentile < 5) |> group_by(Task) |> summarize(Mean_RMSD=mean(RMSD), SD_RMSD=sd(RMSD), Min_RMSD=min(RMSD), Max_RMSD=max(RMSD))
RMSD_10 <- loweffortoptimaldata_RMSD |> dplyr::filter(Sample != -999 & effort2percentile < 10) |> group_by(Task) |> summarize(Mean_RMSD=mean(RMSD), SD_RMSD=sd(RMSD), Min_RMSD=min(RMSD), Max_RMSD=max(RMSD))

# Summary statistics on absolute difference (optimal vs. low effort): 1st, 5th, 10th percentiles (for each task)
AbsDif_1 <- loweffortoptimaldata |> dplyr::filter(Sample != -999 & effort2percentile < 1) |> group_by(Task, Muscle) |> summarize(Mean_Dif=mean(abs(Difference)), SD_Dif=sd(abs(Difference)), Min_Dif=min(abs(Difference)), Max_Dif=max(abs(Difference)))
AbsDif_5 <- loweffortoptimaldata |> dplyr::filter(Sample != -999 & effort2percentile < 5) |> group_by(Task, Muscle) |> summarize(Mean_Dif=mean(abs(Difference)), SD_Dif=sd(abs(Difference)), Min_Dif=min(abs(Difference)), Max_Dif=max(abs(Difference)))
AbsDif_10 <- loweffortoptimaldata |> dplyr::filter(Sample != -999 & effort2percentile < 10) |> group_by(Task, Muscle) |> summarize(Mean_Dif=mean(abs(Difference)), SD_Dif=sd(abs(Difference)), Min_Dif=min(abs(Difference)), Max_Dif=max(abs(Difference)))

# Summary statistics on difference in effort2 (optimal vs. low effort): 1st, 5th, 10th percentiles (for each task)
EffortDif_1 <- loweffortoptimaldata |> dplyr::filter(Sample != -999 & effort2percentile < 1) |> group_by(Task) |> summarize(Mean_Dif=mean(effortdif), SD_Dif=sd(effortdif), Min_Dif=min(effortdif), Max_Dif=max(effortdif))
EffortDif_5 <- loweffortoptimaldata |> dplyr::filter(Sample != -999 & effort2percentile < 5) |> group_by(Task) |> summarize(Mean_Dif=mean(effortdif), SD_Dif=sd(effortdif), Min_Dif=min(effortdif), Max_Dif=max(effortdif))
EffortDif_10 <- loweffortoptimaldata |> dplyr::filter(Sample != -999 & effort2percentile < 10) |> group_by(Task) |> summarize(Mean_Dif=mean(effortdif), SD_Dif=sd(effortdif), Min_Dif=min(effortdif), Max_Dif=max(effortdif))

# Write summary statistics to file
write.csv(RMSD_1, "statistics/loweffort1vsoptimal_rmsd.csv")
write.csv(RMSD_5, "statistics/loweffort5vsoptimal_rmsd.csv")
write.csv(RMSD_10, "statistics/loweffort10vsoptimal_rmsd.csv")
write.csv(AbsDif_1, "statistics/loweffort1vsoptimal_absdif.csv")
write.csv(AbsDif_5, "statistics/loweffort5vsoptimal_absdif.csv")
write.csv(AbsDif_10, "statistics/loweffort10vsoptimal_absdif.csv")
write.csv(EffortDif_1, "statistics/loweffort1vsoptimal_effortdif.csv")
write.csv(EffortDif_5, "statistics/loweffort5vsoptimal_effortdif.csv")
write.csv(EffortDif_10, "statistics/loweffort10vsoptimal_effortdif.csv")

# Clean work space
gdata::keep(feasibledata, reservenames, costs_colnames, thoracoscapular ,glenohumeral, thoracohumeral, elbow, sure = TRUE)

# Tidy feasible data to remove unnecessary columns and pivot longer so muscles aren't separate columns
feasibledata_long <- feasibledata |> select(-all_of(c(reservenames, costs_colnames))) |> pivot_longer(cols=-c(Sample, Task, DOF), names_to = 'Muscle', values_to = 'Activation')

# Multiplying activation by 100 (scale from 0 to 100 instead of 0 to 1; easier to interpret (IMO))
feasibledata_long$Activation <- feasibledata_long$Activation * 100

# Adding in a column for task intensity and direction
feasibledata_long <- feasibledata_long |> mutate(Intensity = case_when(Task %in% c("task02", "task05", "task08", "task11") ~ 20,
                                                                       Task %in% c("task03", "task06", "task09", "task12") ~ 50,
                                                                       Task %in% c("task04", "task07", "task10", "task13") ~ 100,
                                                                       TRUE ~ 0))
feasibledata_long <- feasibledata_long |> mutate(Direction = case_when(Task %in% c("task02", "task03", "task04") ~ 'Down',
                                                                       Task %in% c("task05", "task06", "task07") ~ 'Up',
                                                                       Task %in% c("task08", "task09", "task10") ~ 'Left',
                                                                       Task %in% c("task11", "task12", "task13") ~ 'Right',
                                                                       TRUE ~ 'Unloaded'))

# Summarizing muscle activations for each task x dof combination
feasibledata_activation_amplitude_summary <- feasibledata_long |> group_by(Task, DOF, Muscle) |> summarize(M=mean(Activation), SD=sd(Activation), Median=median(Activation), Min=min(Activation), Max=max(Activation))

# Determining the activation ranges (min and max (LB, UB)) for each muscle (for each task x DOF combination)
feasibledata_activationrange <- feasibledata_long |> group_by(Task, DOF, Muscle, Direction, Intensity) |> summarize(LB = min(Activation), UB = max(Activation))
feasibledata_activationrange$Range <- feasibledata_activationrange$UB - feasibledata_activationrange$LB

# Summarizing counts of activation ranges nearly spanning physiological limits (0-1) and median / mean activation ranges
# Removing Trapezius clavicle since not involved in task
feasibledata_activationrange_summary <- feasibledata_activationrange |> group_by(Task, DOF) |> summarize(n95 = sum(Range > 95, na.rm = TRUE), n=n(), MedianRange=median(Range), MeanRange=mean(Range))

# Determine other metrics of activation range based on percentiles (middle 50, 80, 90% of the data)
feasibledata_activation_quantiles <- feasibledata_long |> group_by(Task, DOF, Muscle) |> reframe(enframe(quantile(Activation, c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1)), "Percentile", "Activation"))
feasibledata_activation_quantiles <- feasibledata_activation_quantiles |> pivot_wider(names_from=Percentile, values_from=Activation)
feasibledata_activation_quantiles$range100 <- feasibledata_activation_quantiles$`100%` - feasibledata_activation_quantiles$`0%`
feasibledata_activation_quantiles$range90 <- feasibledata_activation_quantiles$`95%` - feasibledata_activation_quantiles$`5%`
feasibledata_activation_quantiles$range80 <- feasibledata_activation_quantiles$`90%` - feasibledata_activation_quantiles$`10%`
feasibledata_activation_quantiles$range50 <- feasibledata_activation_quantiles$`75%` - feasibledata_activation_quantiles$`25%`

# Adding in muscle groups
feasibledata_activation_quantiles <- feasibledata_activation_quantiles |> 
  mutate(musclegroup = case_when(Muscle %in% thoracoscapular ~ 'thoracoscapular',
                                 Muscle %in% glenohumeral ~ 'glenohumeral',
                                 Muscle %in% thoracohumeral ~ 'thoracohumeral',
                                 Muscle %in% elbow ~ 'elbow'))

# Summarizing metrics of activation range based on percentiles
feasibledata_quantiles_summarytaskmuscle <- feasibledata_activation_quantiles |> dplyr::filter(DOF %in% c('8dof')) |> group_by(Task, musclegroup) |> summarize(M_range50=mean(range50), M_range80=mean(range80), M_range90=mean(range90), M_range100=mean(range100))
feasibledata_quantiles_summarytask <- feasibledata_activation_quantiles |> dplyr::filter(DOF %in% c('8dof')) |> group_by(Task) |> summarize(M_range50=mean(range50), M_range80=mean(range80), M_range90=mean(range90), M_range100=mean(range100))

# Write summary statistics to file
write.csv(feasibledata_activation_amplitude_summary, "statistics/feasibleactivations_amplitude.csv")
write.csv(feasibledata_activationrange_summary, "statistics/feasibleactivations_range.csv")
write.csv(feasibledata_quantiles_summarytaskmuscle, "statistics/feasibleactivations_quantiles_task_musclegroup.csv")
write.csv(feasibledata_quantiles_summarytask, "statistics/feasibleactivations_quantiles_task.csv")