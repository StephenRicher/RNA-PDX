library('dplyr')
library('sleuth')
library('edgeR')
library('biomaRt')
library('ggplot2')
library('ggpubr')

runSleuth = function(metadata, name, num_cores=1) {
  # Create subdirectory for named analysis
  subdir = paste(sleuthDir, name, "/" ,sep="")
  print(paste0('Saving results to ', subdir))
  dir.create(subdir, showWarnings=FALSE)
  
  so <- sleuth_prep(metadata, target_mapping=t2g, aggregation_column='ens_gene',
                    read_bootstrap_tpm=TRUE, extra_bootstrap_summary=TRUE)
  

  if (name == 'treatment') {
    so <- sleuth_fit(so, ~patient, 'reduced')
    so <- sleuth_fit(so, ~patient + treatment, 'full')
    testName = 'treatmentCTRL'
  } else {
    so <- sleuth_fit(so, ~treatment, 'reduced')
    so <- sleuth_fit(so, ~treatment + response, 'full')
    testName = 'responseresponsive'
  }
  
  so <- sleuth_lrt(so, 'reduced', 'full')
  so <- sleuth_wt(so, testName)
  sleuth_save(so, paste(subdir, name, '.rds', sep=''))
  
  #### PLOTS ####
  which_model = "full"
  which_var <- so$fits[[which_model]]$which_var
  x_label = ifelse(which_var == "obs_counts",  "counts", "tpms")
  sig_level = 0.01 # Threshold to colour plots
  
  ## PCA by group ##
  pcaData = plot_pca(so, color_by='group')$data
  pcaPlot = ggscatter(
    data=pcaData , x='PC1', y='PC2', color='patient', shape=name, size=5) +
    scale_color_brewer(palette = "Dark2")
  ggsave(paste(subdir, name, '-PCA.png', sep=''), pcaPlot, dpi = 300)
  
  ## Mean-Variance relationship ##
  meanVarData = plot_mean_var(so, which_model=which_model)$data
  meanVarData$sqsq_sigma_sq_pmax = sqrt(sqrt(meanVarData$sigma_sq_pmax))
  meanVarPlot = ggscatter(
    data=meanVarData, x='mean_obs', y="sqsq_sigma_sq_pmax", color='iqr',
    palette="Dark2", alpha=0.1, legend.title="Used in shrinkage estimation",
    xlab="mean( log( counts + 0.5 ) )", ylab="sqrt( sigma )",
    title=paste("Mean-Variance relationship of transcripts -", name))
  ggsave(paste(subdir, name, '-meanVar.png', sep=''), meanVarPlot, dpi=300)
  
  # Hack to fix but in plot_ma and plot_volcano - https://github.com/pachterlab/sleuth/issues/233
  so$pval_aggregate = FALSE

  ## MA for wald test ##
  maData = plot_ma(so, test=testName, test_type="wt", which_model=which_model, sig_level=sig_level)$data
  maDataPlot = ggscatter(
    data=maData , x='mean_obs', y="b", color="significant", palette="Dark2", alpha=0.5,
    legend.title=paste('Significant (Q < ', sig_level, ')', sep=''),
    xlab=paste('mean( log(', x_label, '+ 0.5 ) )'), ylab="Beta", 
    title=paste("MA plot - Treatment vs. Control (Wald Test) -", name))
  ggsave(paste(subdir, name, '-MAplot.png', sep=''), maDataPlot, dpi=300)
  
  ## Volcano plot for wald test (full model) ##
  volcanoData = plot_volcano(so, test=testName, test_type="wt", which_model=which_model, sig_level=sig_level)$data
  volcanoData$minusLog10Q = -log10(volcanoData$qval)
  volcanoPlot = ggscatter(
    data=volcanoData , x="b", y="minusLog10Q", color="significant",
    legend.title=paste('Significant (Q < ', sig_level, ')', sep=''),
    xlab=paste("Beta"), ylab="-log10(Qval)", palette="Dark2", alpha=0.5,
    title=paste("Volcano plot - Treatment vs. Control (Wald Test) -", name)) +
    geom_vline(xintercept=0, colour='black', linetype='longdash')
  ggsave(paste(subdir, name, '-volcanoPlot.png', sep=''), volcanoPlot, dpi=300)
}

mart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL",
                         dataset="hsapiens_gene_ensembl",
                         host="https://dec2016.archive.ensembl.org")
t2g <- biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id",
                                   "external_gene_name"), mart=mart)
t2g <- dplyr::rename(t2g, target_id=ensembl_transcript_id,
                     ens_gene=ensembl_gene_id, ext_gene=external_gene_name)
saveRDS(t2g, 'annotation/ensemblHumanTranscript2Gene.rds')
t2g = readRDS('annotation/ensemblHumanTranscript2Gene.rds')

metadata = read.table('../config/sleuth-table.tsv', 
                      header=TRUE, colClasses="character")
metadata$group = paste(metadata$treatment, metadata$patient, sep='-')
metadata$path = paste0('../analysis/kallisto/', metadata$sample)

# Run model of treatment controlled by 'patient'
runSleuth(metadata, "treatment")

# Run model of 'response' controlled by 'treatment'
runSleuth(metadata, "response")
