library('dplyr')
library('sleuth')

diffexp = function(name, subdir="./") {

  so = sleuth_load(paste(subdir, "/", name, '.rds', sep=''))
  
  # Obtain gene-level differential expression results for LRT
  sleuthTableGene = sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
  write.csv(sleuthTableGene, file=paste0(subdir, '/', name, '-gene.csv'), row.names=FALSE)
  
  # Obtain consistent transcript-level differential expression results
  sleuthTableTx <- sleuth_results(so, 'reduced:full', 'lrt', show_all=FALSE, pval_aggregate=FALSE)

  # Retrieve Wold test results with beta values
  dm = so$fits[["full"]]$design_matrix
  test = colnames(dm)[ncol(dm)]
  sleuthTableWold <- sleuth_results(so, test, 'wt', show_all=FALSE, pval_aggregate=FALSE)
  
  # Merge Beta values to transcript level LRT analysis
  sleuthTableTx = merge(sleuthTableTx, sleuthTableWold[c('target_id', 'b', 'se_b')], by='target_id')
  write.csv(sleuthTableTx, file=paste0(subdir, '/', name, '-Tx.csv'), row.names=FALSE)
  
  # Write filtered transcript (q < 0.05)
  write.csv(sleuthTableTx[sleuthTableTx$qval < 0.05,], file=paste0(subdir, '/', name, '-Tx-q0.05.csv'), row.names=FALSE)
  
  # Write normalised count data for each transcript for post-processing
  write.csv(so$obs_norm, file=gzfile(paste0(subdir, '/', name, '-obsNormCounts.csv.gz')), row.names=FALSE)
  
}

for (name in c('treatment', 'response')) {
  diffexp(so, name, name)
}

metadata = read.table(
  '../config/sleuth-table.tsv',  header=TRUE, colClasses="character")
for (patient in unique(metadata$patient)) {
  subdir = paste("patient/", patient, "/" ,sep="")
  diffexp(so, subdir, patient)
}