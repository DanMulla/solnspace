# Principal component analysis of the sampled solution spaces
#
# Two quantities are computed for every task x model combination:
#
#   1. ncomp90/95/99 - the number of principal components required to explain
#      >90/95/99% of the variance in the 42 muscle activations of the sampled
#      solutions. This is an empirical description of how concentrated the
#      solution space is. PCA is run on the correlation matrix (prcomp with
#      scale.=TRUE) so that every muscle contributes equally regardless of the
#      width of its own activation range; the covariance version is reported
#      alongside for comparison.
#
#   2. dim_theoretical - the exact affine dimension of the muscle activation
#      space implied by the mechanical equilibrium constraints alone, i.e.
#      42 - rank(R*F) over the balanced degrees of freedom. Each independent
#      DOF removes exactly one dimension. This is algebra, not statistics, and
#      does not depend on the sampling.
#
# Outputs: statistics/pca_ncomponents.csv
#          statistics/pca_cumulativevariance.csv
#          statistics/Fig-PCA-CumulativeVariance.png

# Libraries
library(conflicted)
library(tidyverse)

# Clear work space
rm(list = ls())

# Set working directory to current folder
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Output folder for the summary statistics written below (write.csv and ggsave
# do not create missing directories)
dir.create("statistics", showWarnings = FALSE)

nmuscles <- 42                      # muscle elements, excluding the conoid ligament
dofs <- c('3dof', '5dof', '8dof', '11dof')
dof_rows <- list('3dof' = 13:15,    # rows of the 19-DOF moment matrix that each
                 '5dof' = 13:17,    # model balances (see feasiblesolns.m)
                 '8dof' = 10:17,
                 '11dof' = 7:17)
ndof_map <- c('3dof' = 3, '5dof' = 5, '8dof' = 8, '11dof' = 11)

ncomp_for <- function(pc, threshold) {
  cumvar <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  which(cumvar > threshold)[1]
}

results <- list()
cumvar_long <- list()
taskfolders <- list.dirs('outputs/.', full.names = FALSE, recursive = FALSE)

for (taskname in taskfolders) {

  # Muscle moment matrix for this task (19 DOF x 43 force elements)
  momentfile <- file.path('outputs', taskname, paste(taskname, '-maxmoments.csv', sep=""))
  if (!file.exists(momentfile)) next
  musclemoments <- as.matrix(read.csv(momentfile, header = FALSE))
  musclemoments[abs(musclemoments) < 1e-5] <- 0      # as in feasiblesolns.m
  musclemoments <- musclemoments[, -1]               # drop the conoid ligament

  for (dof in dofs) {

    samplefile <- file.path('outputs', taskname, 'solutions',
                            paste(taskname, '-samples-', dof, '-all.txt', sep=""))
    if (!file.exists(samplefile)) next

    # Sampled solutions: keep the muscle columns only, drop the reserve actuators
    samples <- as.matrix(read.delim(samplefile, header = FALSE))[, 1:nmuscles]

    # 1. empirical PCA
    pc_cor <- prcomp(samples, center = TRUE, scale. = TRUE)
    pc_cov <- prcomp(samples, center = TRUE, scale. = FALSE)

    # 2. exact dimension implied by the equilibrium constraints alone
    rows <- dof_rows[[dof]]
    rank_moments <- qr(musclemoments[rows, , drop = FALSE])$rank
    dim_theoretical <- nmuscles - rank_moments

    results[[length(results) + 1]] <- tibble(
      Task = taskname, DOF = dof, nDOF = ndof_map[[dof]], nSamples = nrow(samples),
      ncomp90 = ncomp_for(pc_cor, 0.90),
      ncomp95 = ncomp_for(pc_cor, 0.95),
      ncomp99 = ncomp_for(pc_cor, 0.99),
      ncomp90_cov = ncomp_for(pc_cov, 0.90),
      PC1_pct = 100 * pc_cor$sdev[1]^2 / sum(pc_cor$sdev^2),
      rank_moments = rank_moments,
      dim_theoretical = dim_theoretical)

    cumvar_long[[length(cumvar_long) + 1]] <- tibble(
      Task = taskname, DOF = dof, PC = seq_len(nmuscles),
      CumVar = 100 * cumsum(pc_cor$sdev^2) / sum(pc_cor$sdev^2))

    rm(samples, pc_cor, pc_cov)
  }
}

pca_results <- bind_rows(results)
pca_cumvar  <- bind_rows(cumvar_long)
pca_results$DOF <- factor(pca_results$DOF, levels = dofs)
pca_cumvar$DOF  <- factor(pca_cumvar$DOF,  levels = dofs)

# Summary across tasks
pca_summary <- pca_results |> group_by(DOF, nDOF) |>
  summarize(n = n(),
            Median_ncomp90 = median(ncomp90), Min_ncomp90 = min(ncomp90), Max_ncomp90 = max(ncomp90),
            Median_ncomp90_cov = median(ncomp90_cov),
            dim_theoretical = median(dim_theoretical), .groups = 'drop')
print(as.data.frame(pca_summary))

write.csv(pca_results, "statistics/pca_ncomponents.csv", row.names = FALSE)
write.csv(pca_cumvar,  "statistics/pca_cumulativevariance.csv", row.names = FALSE)

# Cumulative variance curves, one line per task, panelled by model
p <- ggplot(pca_cumvar, aes(x = PC, y = CumVar, group = Task)) +
  geom_hline(yintercept = 90, linetype = 'dashed', colour = 'grey40') +
  geom_line(alpha = 0.5, colour = '#0077BB') +
  facet_wrap(~DOF, nrow = 1) +
  labs(x = 'Principal component', y = 'Cumulative variance explained (%)') +
  scale_x_continuous(breaks = seq(0, 42, 10)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", linetype = 0),
        strip.text = element_text(colour = 'black', size = 13),
        axis.text = element_text(size = 11, color = 'black'),
        axis.title = element_text(size = 14, face = 'bold'))
ggsave("statistics/Fig-PCA-CumulativeVariance.png", p, height = 3, width = 10, dpi = 600)
