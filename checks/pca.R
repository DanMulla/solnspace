# Two quantities per task x model:
#
#   1. ncomp90/95/99 - how many principal components are needed to explain
#      >90/95/99% of the variance in the 42 muscle activations of the sampled
#      solutions
#
#   2. dim_theoretical - the exact affine dimension of the muscle activation
#      space implied by the mechanical equilibrium constraints alone, i.e.
#      42 - rank(moment matrix).

rm(list = ls())

# Set working directory to script's folder, then read from ../outputs
if (requireNamespace('rstudioapi', quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
outroot <- file.path('..', 'outputs')
if (!dir.exists(outroot)) stop('Cannot find ', normalizePath(outroot, mustWork = FALSE), ' - run this script from the checks/ folder.')

nmuscles <- 42                      
dofs     <- c('3dof', '5dof', '8dof', '11dof')
dof_rows <- list('3dof' = 13:15, '5dof' = 13:17, '8dof' = 10:17, '11dof' = 7:17)
ndof_map <- c('3dof' = 3, '5dof' = 5, '8dof' = 8, '11dof' = 11)

ncomp_for <- function(pc, threshold) {
  cumvar <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
  which(cumvar > threshold)[1]
}

taskfolders <- sort(list.dirs(outroot, full.names = FALSE, recursive = FALSE))
res <- data.frame()

for (taskname in taskfolders) {
  momentfile <- file.path(outroot, taskname, paste0(taskname, '-maxmoments.csv'))
  if (!file.exists(momentfile)) next
  musclemoments <- as.matrix(read.csv(momentfile, header = FALSE))
  musclemoments[abs(musclemoments) < 1e-5] <- 0
  musclemoments <- musclemoments[, -1]

  for (dof in dofs) {
    samplefile <- file.path(outroot, taskname, 'solutions',
                            paste0(taskname, '-samples-', dof, '-all.txt'))
    if (!file.exists(samplefile)) next

    samples <- as.matrix(read.delim(samplefile, header = FALSE))[, 1:nmuscles]
    pc_cor <- prcomp(samples, center = TRUE, scale. = TRUE)
    pc_cov <- prcomp(samples, center = TRUE, scale. = FALSE)

    rank_moments <- qr(musclemoments[dof_rows[[dof]], , drop = FALSE])$rank

    res <- rbind(res, data.frame(
      Task = taskname, DOF = dof, nDOF = ndof_map[[dof]], nSamples = nrow(samples),
      ncomp90 = ncomp_for(pc_cor, 0.90),
      ncomp95 = ncomp_for(pc_cor, 0.95),
      ncomp99 = ncomp_for(pc_cor, 0.99),
      ncomp90_cov = ncomp_for(pc_cov, 0.90),
      PC1_pct = 100 * pc_cor$sdev[1]^2 / sum(pc_cor$sdev^2),
      rank_moments = rank_moments,
      dim_theoretical = nmuscles - rank_moments,
      stringsAsFactors = FALSE))

    rm(samples, pc_cor, pc_cov)
  }
}

if (nrow(res) == 0) stop('No sampled solutions found under ', normalizePath(outroot))
res$DOF <- factor(res$DOF, levels = dofs)

line <- function(ch = '-', n = 74) cat(strrep(ch, n), '\n', sep = '')

cat('\n'); line('=')
cat('PCA OF THE SAMPLED SOLUTION SPACES\n')
cat(sprintf('%d task x model combinations, %d muscle activations each\n',
            nrow(res), nmuscles))
line('=')

## 1. components needed for >90% variance, per task and model
cat('\nComponents needed to explain >90% of variance (of ', nmuscles, ')\n\n', sep = '')
cat(sprintf('%-9s', 'task')); cat(sprintf('%7s', dofs)); cat('\n')
for (t in unique(res$Task)) {
  cat(sprintf('%-9s', t))
  for (d in dofs) {
    v <- res$ncomp90[res$Task == t & res$DOF == d]
    cat(sprintf('%7s', if (length(v)) as.character(v) else '-'))
  }
  cat('\n')
}
cat('\n("-" = no feasible solution for that task and model)\n')

## 2. summary by model
cat('\n'); line()
cat('SUMMARY BY MODEL\n'); line()
cat(sprintf('\n%-7s%5s%7s%22s%9s%9s%8s%14s\n', 'model', 'nDOF', 'tasks',
            'ncomp90 med[min-max]', 'ncomp95', 'ncomp99', 'PC1 %', 'theoretical'))
for (d in dofs) {
  s <- res[res$DOF == d, ]
  if (!nrow(s)) next
  cat(sprintf('%-7s%5d%7d%22s%9.0f%9.0f%7.1f%%%14d\n',
              d, s$nDOF[1], nrow(s),
              sprintf('%.0f [%d-%d]', median(s$ncomp90), min(s$ncomp90), max(s$ncomp90)),
              median(s$ncomp95), median(s$ncomp99), median(s$PC1_pct),
              median(s$dim_theoretical)))
}

## 3. equilibrium alone vs the full constraint set
cat('\n'); line()
cat('EQUILIBRIUM ALONE vs THE SAMPLED SPACE\n'); line()
cat('\n"theoretical" counts only the dimensions equilibrium removes. The gap below\n')
cat('is the further concentration produced by the activation bounds and the\n')
cat('glenoid stability constraint, plus the fact that variance is not uniform.\n\n')
cat(sprintf('%-7s%14s%14s%10s\n', 'model', 'theoretical', 'ncomp90 (med)', 'gap'))
for (d in dofs) {
  s <- res[res$DOF == d, ]
  if (!nrow(s)) next
  th <- median(s$dim_theoretical); em <- median(s$ncomp90)
  cat(sprintf('%-7s%14d%14.0f%10.0f\n', d, th, em, th - em))
}

## 4. sanity checks
cat('\n'); line()
cat('CHECKS\n'); line('-')
bad_rank <- res[res$rank_moments != res$nDOF, ]
cat(sprintf('  moment matrix full rank for every model (rank == nDOF) : %s\n',
            if (nrow(bad_rank) == 0) 'yes' else
              paste0('NO - ', nrow(bad_rank), ' rows differ')))
cat(sprintf('  ncomp90 <= theoretical dimension everywhere             : %s\n',
            if (all(res$ncomp90 <= res$dim_theoretical)) 'yes' else
              paste0('NO - ', sum(res$ncomp90 > res$dim_theoretical), ' exceed it')))
cat(sprintf('  ncomp90 <= ncomp95 <= ncomp99 everywhere                : %s\n',
            if (all(res$ncomp90 <= res$ncomp95 & res$ncomp95 <= res$ncomp99)) 'yes' else 'NO'))
cat(sprintf('  samples per combination                                 : %d - %d\n',
            min(res$nSamples), max(res$nSamples)))
