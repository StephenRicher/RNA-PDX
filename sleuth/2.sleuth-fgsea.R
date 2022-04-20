library('fgsea')
library('dplyr')
library('tidyverse')
library('data.table')

level = 'gene'
test_stat = "log_qval"

for (name in c('treatment', 'response')) {
  deg = read.csv(paste0(name, '/', name, '-', level, '.csv'))
  deg$log_qval = -log(deg$qval)
  
  # Remove any genes with no valid gene symbol mapping and
  # remove any genes with duplicate symbols for different IDs
  res2 <- deg %>% 
    dplyr::select(ext_gene, !!as.symbol(test_stat)) %>%
    na_if("") %>%
    na.omit() %>% 
    group_by(ext_gene) %>%
    filter(n() == 1) %>%
    arrange(by=!!as.symbol(test_stat))
  ranks <- deframe(res2)
  
  
  pathways <- gmtPathways('annotation/ReactomePathways.gmt')
  fgseaRes <- fgseaSimple(pathways=pathways, stats=ranks, nperm=10000, scoreType='pos')
  out = paste0(name, '/', name, '-',  level, '-ReactomePathways-FGSEA.csv')
  fwrite(fgseaRes, file=out, sep="\t", sep2=c("", " ", ""))
  
  pathways <- gmtPathways('annotation/h.all.v7.4.symbols.gmt')
  fgseaRes <- fgseaSimple(pathways=pathways, stats=ranks, nperm=10000, scoreType='pos')
  out = paste0(name, '/', name, '-', level, '-Hallmark-FGSEA.csv')
  fwrite(fgseaRes, file=out, sep="\t", sep2=c("", " ", ""))
}

