setwd('~/d/sci/src/prnp_penetrance')
source('../exac_papers/exac_constants.R')
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}

# get an_adj from ExAC callset to use in supp table 3
prnp_calls = exac[exac$symbol=='PRNP' & !is.na(exac$symbol),c("chrom","pos","ref","alt","ac_adj","an_adj","pos_id")]
write.table(prnp_calls, 'data_nosync/exac_prnp_calls.tsv', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
