options(stringsAsFactors=FALSE)
require(sqldf)
setwd('~/d/sci/src/prnp_penetrance')

percent = function(number,digits=0) {
  return (gsub(' ','',paste(formatC(number*100,format='f',digits=digits),'%',sep='')))
}

# get real column names to avoid things like "X2.OPRD"
surveillance_colnames = as.character(read.table('data_nosync/surveillance_data_raw.tsv',sep='\t',header=FALSE)[1,])

surveillance = read.table('data_nosync/surveillance_data_raw.tsv',sep='\t',header=TRUE)
colnames(surveillance) = surveillance_colnames
surveillance = rbind(surveillance,rep(NA,dim(surveillance)[2]))
nrow = dim(surveillance)[1]
ncol = dim(surveillance)[2]
surveillance$country[nrow] = 'TOTAL'
surveillance[nrow,6:ncol] = colSums(surveillance[1:(nrow-1),6:ncol],na.rm=TRUE)
surveillance$definite_plus_probable[nrow] = sum(surveillance$definite_plus_probable,na.rm=TRUE)
surveillance$prnp_sequenced[nrow] = sum(surveillance$prnp_sequenced,na.rm=TRUE)
surveillance$proportion_sequenced = percent(surveillance$prnp_sequenced / surveillance$definite_plus_probable)
surveillance$with_rare_variants = rowSums(surveillance[,6:ncol],na.rm=TRUE)
surveillance$proportion_with_rare_variants = percent(surveillance$with_rare_variants / surveillance$definite_plus_probable)
surveillance = surveillance[,c(1:5,ncol+(1:3),6:ncol)]
colnames(surveillance)[10] = '1-OPRI'
write.table(surveillance,'supplement/table_s01_surveillance_allele_counts.tsv',sep='\t',na='',col.names=TRUE,row.names=FALSE,quote=FALSE)

surveillance = t(surveillance)

write.table(surveillance,'data_nosync/table_s01_to_copy_into_document.tsv',sep='\t',na='',col.names=FALSE,row.names=TRUE,quote=FALSE)


# raw screenshot annotations from Google spreadsheet
annos = read.table('data_nosync/screenshots_annotated.tsv',header=TRUE,sep='\t',quote='',comment.char='')
# the above screenshot annotations were done on release v0.1, which had ~3,000 additional
# individuals who later removed upon further QC. Therefore, restrict the list to
# samples included in v0.3 release of 60,706 individuals
list60k = read.table('data_nosync/60.5k.list',header=FALSE,sep='\t')
exac_prnp_calls = read.table('data_nosync/exac_prnp_calls.tsv',header=TRUE,sep='\t')

# spaces to underscores in sample list
list60k$sample = gsub(" ","_",list60k$V1)

write.table(list60k[,"sample"],'data_nosync/list_60_5k_no_space.tsv',col.names=FALSE,row.names=FALSE,quote=FALSE)

annos$sampleid = gsub(",.*","",gsub("Call(sample=","",annos$details,fixed=TRUE))
annos$sampleid = gsub(" ","_",annos$sampleid)
annos$hgvs_pos = as.integer(gsub('[^0-9].*','',substr(annos$HGVS,3,15)))
annos$chrom = gsub(":.*","",annos$allele_id)
annos$pos = as.integer(gsub(".*:","",gsub("_.*","",annos$allele_id)))
annos$ref = gsub(".*_","",gsub(">.*","",annos$allele_id))
annos$alt = gsub(".*>","",annos$allele_id)
annos = annos[annos$sampleid %in% list60k$sample,]
# Monkol's suggestion - join in AN_Adj from ExAC
annos$pos_id = paste(annos$chrom, formatC(annos$pos,width=9,flag='0'), annos$ref, annos$alt, sep='_')
annos$an_adj = exac_prnp_calls$an_adj[match(annos$pos_id,exac_prnp_calls$pos_id)]
annos$call_rate = percent(annos$an_adj / (2*60706))

annos$class = ''
for (i in 1:dim(annos)[1]) {
  first_char = substr(annos$aa_code[i],1,1)
  last_char = substr(annos$aa_code[i],nchar(annos$aa_code[i]),nchar(annos$aa_code[i]))
  if (nchar(annos$aa_code[i]) == 0) {
    annos$class[i] = 'non-coding'
  } else if (first_char =='X') {
    annos$class[i] = 'read-through'
  } else if(first_char == last_char) {
    annos$class[i] = 'synonymous'
  } else if (last_char == 'X') {
    annos$class[i] = 'nonsense'
  } else if (nchar(annos$class[i]) > 5) {
    annos$class[i] = 'indel'
  } else {
    annos$class[i] = 'missense'
  }
}

anno_summary = sqldf("
select   chrom, pos, ref, alt, HGVS hgvs, aa_code variant, class, call_rate, count(*) ac
from     annos a
where    exact_gt_correct = 'y'
group by 1,2,3,4,5,6,7,8
order by 1,2
;")

write.table(anno_summary,'supplement/table_s03_exac_allele_counts.tsv',sep='\t',quote=FALSE,row.names=FALSE)

anno_summary_summary = sqldf("
select   class, sum(ac) total_ac
from     anno_summary
group by 1
order by 1
;")

write.table(anno_summary_summary,'supplement/table_s04_exac_functional_summary.tsv',sep='\t',quote=FALSE,row.names=FALSE)

reportedly_pathogenic = read.table('supplement/table_s02_reportedly_pathogenic_variants.tsv',header=TRUE,sep='\t',quote='',comment.char='')

indivs_with_path_alleles = sqldf("
select   a.sampleid
from     annos a
where    a.aa_code in (select variant from reportedly_pathogenic)
and      a.exact_gt_correct = 'y'
;")

exac_pops = read.table('data_nosync/exac_60706_pops.tsv',header=FALSE,sep='\t',quote='',comment.char='')
colnames(exac_pops) = c('sampleid','population_code')

desc_1kg_pops = read.table('data_nosync/1kg_pops.tsv',header=TRUE,sep='\t',quote='',comment.char='')

exac_pop_summary = sqldf("
select   d.population_code, d.description, d.super_population_code, e.n n_exac
from     desc_1kg_pops d left outer join (
    select   population_code, count(*) n
    from     exac_pops 
    group by 1
    order by 1) e
on       d.population_code = e.population_code
order by 1
;")


path_c129 = read.table('data_nosync/path_allele_codon129.txt',header=FALSE,sep='\t')
colnames(path_c129) = c('sampleid','codon129')

path_pops = sqldf("
select   sampleid, population_code
from     exac_pops
where    sampleid in (select sampleid from indivs_with_path_alleles)
;")

path_c129$aa129 = gsub('G','V',gsub('A','M',path_c129$codon129))

grouped_by_pop = sqldf("
select   a.aa_code, p.population_code, count(*) n
from     path_pops p, annos a
where    p.sampleid = a.sampleid
group by 1,2
;")

summary_by_pop = sqldf("
select   aa_code, group_concat(' '||n||' '||population_code) pops
from     grouped_by_pop
group by 1
;")

summary_by_c129 = sqldf("
select   pos, aa_code, group_concat(' '||n||' '||aa129) codon129
from (
select   a.pos, a.aa_code, c.aa129, count(*) n
from     path_c129 c, annos a
where    c.sampleid = a.sampleid
group by 1,2,3
) subq
group by 1
;")

pop_and_c129_summary = sqldf("
select   p.aa_code variant, p.pops, c.codon129
from     summary_by_pop p, summary_by_c129 c
where    p.aa_code = c.aa_code
order by c.pos
;")

pop_and_c129_summary$use_pop = ''
pop_and_c129_summary$use_pop[pop_and_c129_summary$variant=='M232R'] = 'JPT'
pop_and_c129_summary$use_pop[pop_and_c129_summary$variant=='V180I'] = 'JPT'
pop_and_c129_summary$use_pop[pop_and_c129_summary$variant=='V210I'] = 'TSI'

pop_and_c129_summary$pop_ac[pop_and_c129_summary$variant=='M232R'] = grouped_by_pop$n[grouped_by_pop$population_code=='JPT' & grouped_by_pop$aa_code=='M232R']
pop_and_c129_summary$pop_ac[pop_and_c129_summary$variant=='V180I'] = grouped_by_pop$n[grouped_by_pop$population_code=='JPT' & grouped_by_pop$aa_code=='V180I']
pop_and_c129_summary$pop_ac[pop_and_c129_summary$variant=='V210I'] = grouped_by_pop$n[grouped_by_pop$population_code=='TSI' & grouped_by_pop$aa_code=='V210I']

# note: here i am using the total number of individuals in the population as the "N" or the denominator to which AC
# is compared for penetrance calculations. this is because in JPT and TSI for these three variants of interest, 
# the call rate is exactly 100%. see master.bash where i use summarize_genotypes.py to check this.
pop_and_c129_summary$pop_n[pop_and_c129_summary$variant=='M232R'] = exac_pop_summary$n_exac[exac_pop_summary$population_code=='JPT']
pop_and_c129_summary$pop_n[pop_and_c129_summary$variant=='V180I'] = exac_pop_summary$n_exac[exac_pop_summary$population_code=='JPT']
pop_and_c129_summary$pop_n[pop_and_c129_summary$variant=='V210I'] = exac_pop_summary$n_exac[exac_pop_summary$population_code=='TSI']


write.table(pop_and_c129_summary,"supplement/table_s07_exac_path_ancestry_and_c129.tsv",sep='\t',row.names=FALSE,quote=FALSE)

write.table(exac_pop_summary,"supplement/table_s08_exac_pop_summary.tsv",sep='\t',row.names=FALSE,quote=FALSE)

# QC metrics for screenshot evaluation
sum(annos$exact_gt_correct=='y')/dim(annos)[1]
gt_correct = annos$exact_gt_correct=='y'
hist(annos$this_alt_ab[gt_correct])
hist(annos$this_alt_ad[gt_correct])
hist(annos$this_gt_gq[gt_correct])
sum(annos$this_gt_gq[gt_correct] >= 95,na.rm=TRUE)/sum(gt_correct,na.rm=TRUE)
sum(annos$this_alt_ab[gt_correct] <= .7 & annos$this_alt_ab[gt_correct] >= .3,na.rm=TRUE)/sum(gt_correct,na.rm=TRUE)
sum(annos$this_alt_ad[gt_correct] >= 10,na.rm=TRUE)/sum(gt_correct,na.rm=TRUE)

# Figure S1
meta = read.table('../exac_papers/misc_data/exac_meta_with_age.tsv',header=TRUE,sep='\t')
meta$id = gsub(" ","_",meta$vcf_sampleid)
indivs_with_path_alleles$id = gsub(" ","_",indivs_with_path_alleles$sampleid)
indivs_with_path_alleles$age = meta$age[match(indivs_with_path_alleles$id, meta$id)]
sum(!is.na(indivs_with_path_alleles$age))
cairo_pdf('figures/figure_s1.pdf',width=8,height=6,bg='#FFFFFF00')
par(mar=c(5,5,1,5))
overall_color = '#FF4400'
prnp_color = '#0000CD'
n_age = sum(!is.na(meta$age))
n_age_text = formatC(n_age, big.mark=',')
n_prnp = sum(!is.na(indivs_with_path_alleles$age))
hist(meta$age, breaks=20, prob=TRUE, col=alpha(overall_color,.5), 
     xlim=c(16,87), ylim=c(0,.035),
     border=NA, axes=FALSE, xlab='Age', ylab='', main='')
abline(h=0, col='black')
axis(side=1, at=c(18,30,40,50,60,70,85), labels=c('\u226418', '30', '40', '50', '60', '70', '\u226585'), lwd=0, lwd.ticks=1)
axis(side=2, col.axis = overall_color, col.ticks = overall_color, at=(0:3)/100, labels=percent((0:3)/100), lwd=0, lwd.ticks=1, las=1,)
mtext(side=2, line=3, text=paste('Overall ExAC distribution density, n = ', n_age_text), col=overall_color, cex=.9)
par(new=TRUE)
hist(indivs_with_path_alleles$age, breaks=6, xlim=c(16,87), col=alpha(prnp_color, .5), border=NA, axes=FALSE, xlab='', ylab='', main='')
axis(side=4, at=(0:12), labels=0:12, lwd=0, lwd.ticks=1, las=1, col.axis=prnp_color, col.ticks=prnp_color)
mtext(side=4, line=3, text=paste('Individuals with reportedly pathogenic PRNP variants, n = ',n_prnp,sep=''), col=prnp_color, cex=.9)
dev.off()

# stats to quote in Figure S1 legend
meta$has_prnp_mut = meta$id %in% indivs_with_path_alleles$id
wilcox.test(meta$age[meta$has_prnp_mut], meta$age[!meta$has_prnp_mut])
t.test(meta$age[meta$has_prnp_mut], meta$age[!meta$has_prnp_mut])
m = lm(age ~ consortium + has_prnp_mut, data=meta)
summary(m)
range(meta$age[meta$has_prnp_mut],na.rm=TRUE)
