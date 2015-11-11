options(stringsAsFactors=FALSE)
require(sqldf)
require(binom)

# Table S1. Allele counts of rare PRNP variants in 16,025 definite and probable prion disease cases in 9 countries
ac_case = read.table('supplement/table_s01_surveillance_allele_counts.tsv',header=TRUE,sep='\t',quote='',comment.char='')
# replace NA with 0 in allele counts
ac_case[,9:65][is.na(ac_case[,9:65])] = 0
rownames(ac_case) = ac_case$country

# Table S2. Rare PRNP variants reported in peer-reviewed literature to cause prion disease
reportedly_pathogenic = read.table('supplement/table_s02_reportedly_pathogenic_variants.tsv',header=TRUE,sep='\t',quote='',comment.char='')

# Table S3. Allele counts of rare PRNP variants in 60,706 individuals in ExAC
ac_exac = read.table('supplement/table_s03_exac_allele_counts.tsv',header=TRUE,sep='\t',quote='',comment.char='')

# Table S4. Summary of rare PRNP variants by functional class in ExAC
exac_class_summary = read.table('supplement/table_s04_exac_functional_summary.tsv',header=TRUE,sep='\t',quote='',comment.char='')

# Table S5. Allele counts of 16 reportedly pathogenic PRNP variants in >500,000 23andMe customers
ac_23andme = read.table('supplement/table_s05_23andme_ac.tsv',header=TRUE,sep='\t',quote='',comment.char='')
rownames(ac_23andme) = ac_23andme$variant

# Table S6. Phenotypes investigated in studies in which ExAC individuals with reportedly pathogenic PRNP variants were ascertained
exac_phenos_reportedly_pathogenic = read.table('supplement/table_s06_exac_phenos_reportedly_pathogenic.tsv',header=TRUE,sep='\t',quote='',comment.char='')

# Table S7. Inferred ancestry and codon 129 genotypes of ExAC individuals with reportedly pathogenic variants
exac_pop_c129_summary = read.table('supplement/table_s07_exac_path_ancestry_and_c129.tsv',header=TRUE,sep='\t',quote='',comment.char='')
rownames(exac_pop_c129_summary) = exac_pop_c129_summary$variant

# Table S8. Inferred ancestry of all ExAC individuals
exac_pop_summary = read.table('supplement/table_s08_exac_pop_summary.tsv',header=TRUE,sep='\t',quote='',comment.char='')

# Table S9. Inferred ancestry of 23andMe customers
ancestry_23andme = read.table('supplement/table_s09_23andme_ancestry.tsv',header=TRUE,sep='\t',quote='',comment.char='')

# Table S10. Details of Japanese prion disease cases
japan_details = read.table('supplement/table_s10_japanese_case_details.tsv',fileEncoding='utf8',header=TRUE,sep='\t',quote='',comment.char='')

# Table S11. 
ptv_phenotypes = read.table('supplement/table_s11_truncating_variant_phenotypes.tsv',header=TRUE,sep='\t',quote='',comment.char='')

# Figure 1B.
exac_path_acs = sqldf("
select   variant, ac
from     ac_exac
where    variant in (select variant from reportedly_pathogenic)
order by pos
;")


pdf('figures/figure1b.pdf',width=12,height=3)
k = '#1A6966'
par(mar=c(3,5,2,2))
plot(NA,NA,axes=FALSE,xlab='',ylab='ExAC allele count',ylim=c(0,10),xlim=c(1,dim(exac_path_acs)[1]),cex.lab=1.5)
abline(h=0,col='black',lwd=.5)
points(1:dim(exac_path_acs)[1],exac_path_acs$ac,type='h',lwd=40,lend=1,col=k)
mtext(side=1,at=1:dim(exac_path_acs)[1],text=exac_path_acs$variant,col=k,font=2,cex=1.35)
axis(side=2,at=5*(0:2),labels=5*(0:2),cex.axis=1,las=1,lwd=0,lwd.ticks=1,font=2)
dev.off()

# Figure 2.
case_totals = data.frame(t(ac_case[ac_case$country=="TOTAL",9:dim(ac_case)[2]]))
colnames(case_totals) = 'ac'
case_totals$variant = rownames(case_totals)
case_totals$ac[case_totals$variant == "V180I"] = case_totals$ac[case_totals$variant == "V180I"] + case_totals$ac[case_totals$variant == "V180I_and_M232R_in_trans"]
case_totals$ac[case_totals$variant == "M232R"] = case_totals$ac[case_totals$variant == "M232R"] + case_totals$ac[case_totals$variant == "V180I_and_M232R_in_trans"]
case_control = sqldf("
select   p.variant, 
         case when k.ac is null then 0 else k.ac end ac_case,
         case when e.ac is null then 0 else e.ac end ac_exac
from     reportedly_pathogenic p left outer join case_totals k
on       p.variant = k.variant
left outer join ac_exac e
on       p.variant = e.variant
;")
# specify which mutations to label with text and where to put the labels
label_subset = case_control$variant %in% c("P39L","R148H","V203I","R208C","R208H","M232R","V180I","T188R","V210I","V180I","E200K","A117V","P102L","D178N","E196A")
case_control_text = case_control[label_subset,]
case_control_text$pos = 4
case_control_text$pos[case_control_text$variant=="D178N"] = 2
case_control_text$pos[case_control_text$variant=="E200K"] = 2
case_control_text$pos[case_control_text$variant=="E196A"] = 3
case_control_text$pos[case_control_text$variant=="R148H"] = 1
case_control_text$pos[case_control_text$variant=="T188R"] = 3
case_control_text$pos[case_control_text$variant=="P39L"] = 2
highlight_subset = case_control$variant %in% c("P102L","A117V","D178N","E200K")
pdf('figures/figure2.pdf',width=10,height=4.7)
par(mar=c(5,6,1,1))
plot(NA,NA,axes=FALSE,
     xlim=c(-40,600),ylim=c(0,10.5),
     xlab='cases reported to prion surveillance centers',
     ylab='allele count in ExAC\npopulation controls',
     main='',
     cex.lab=1.5,cex.main=2)
abline(h=0:10,col='#CCCCCC',lwd=.25)
abline(h=0,col='#CCCCCC',lwd=3)
abline(v=(0:6)*100,col='#CCCCCC')
abline(v=0,col='#CCCCCC',lwd=3)
rect(xleft=-100,xright=0,ybottom=-10,ytop=20,col='#FFFFFF',border=NA)
rect(xleft=-100,xright=600,ybottom=-10,ytop=0,col='#FFFFFF',border=NA)
points(case_control$ac_case,case_control$ac_exac,pch=19,cex=1.5)
axis(side=1,at=(0:6)*100,labels=(0:6)*100,lwd=0,lwd.ticks=1)
axis(side=2,at=0:10,labels=0:10,lwd=0,lwd.ticks=1,cex=.5,las=1)
points(case_control$ac_case[highlight_subset],case_control$ac_exac[highlight_subset],cex=2.5,pch=2,col='red')
text(x=case_control_text$ac_case,y=case_control_text$ac_exac,pos=case_control_text$pos,labels=case_control_text$variant,cex=1.2)
rect(xleft=350,xright=700,ybottom=8.1,ytop=11,col='#FFFFFF',border=NA)
legend('topright',c("reportedly pathogenic variant","segregation in multiple multigenerational families;\nspontaneous disease in mouse models"),col=c('black','red'),pt.cex=c(1.5,2),cex=.8,pch=c(19,2),bty='n',bg='#FFFFFF')
dev.off()

# Figure 3.
# formulae for penetrance and 95% confidence intervals as per Kirov et al. 2014
# see http://www.biologicalpsychiatryjournal.com/article/S0006-3223(13)00676-8/pdf
assumed_baseline_risk = 2e-4
penetrance = function(af_case, af_control, baseline_risk=assumed_baseline_risk) {
  calculated_penetrance = af_case * baseline_risk / af_control
  estimated_penetrance = pmin(1,pmax(0,calculated_penetrance)) # trim to [0,1] support
  return (estimated_penetrance)
}
penetrance_confint = function (ac_case, n_case, ac_control, n_control, baseline_risk=assumed_baseline_risk) {
  # for a genotypic model, use 1*n_case; for allelic, use 2*n_case
  # here, results are virtually identical.
  case_confint = binom.confint(x=ac_case,n=2*n_case,method='wilson')
  control_confint = binom.confint(x=ac_control,n=2*n_control,method='wilson')
  lower_bound = penetrance(case_confint$lower,control_confint$upper,baseline_risk)
  best_estimate = penetrance(case_confint$mean,control_confint$mean,baseline_risk)
  upper_bound = penetrance(case_confint$upper,control_confint$lower,baseline_risk)
  return ( c(lower_bound, best_estimate, upper_bound) )
}
wrap = function(text,width) {
  return ( sapply(lapply(strwrap(as.list(text),width=width,simplify=FALSE),paste,collapse="\n"),"[[",1) )
}

mendelians = c("P102L","A117V","D178N","E200K")
forest_data = data.frame(yvals=c(3.75,3.25,2.75,2.25,1.5,0.75,0.25),
                         comparisons=c("M232R_ExAC","M232R_23andMe","V180I_ExAC","V180I_23andMe","V210I_ExAC","Mendelians_ExAC","Mendelians_23andMe"))
forest_data$ac_case = 0
forest_data$ac_case[1:2] = ac_case["Japan","M232R"] + ac_case["Japan","V180I_and_M232R_in_trans"]
forest_data$ac_case[3:4] = ac_case["Japan","V180I"] + ac_case["Japan","V180I_and_M232R_in_trans"]
forest_data$ac_case[5] = ac_case["Italy","V210I"]
forest_data$ac_case[6:7] = sum(ac_case["TOTAL",mendelians])
forest_data$n_case = 0
forest_data$n_case[1:4] = ac_case["Japan","prnp_sequenced"]
forest_data$n_case[5] = ac_case["Italy","prnp_sequenced"]
forest_data$n_case[6:7] = ac_case["TOTAL","prnp_sequenced"]
forest_data$ac_control = 0
forest_data$ac_control[1] = exac_pop_c129_summary["M232R","pop_ac"]
forest_data$ac_control[2] = ac_23andme["M232R","pop_ac_use"]
forest_data$ac_control[3] = exac_pop_c129_summary["V180I","pop_ac"]
forest_data$ac_control[4] = ac_23andme["V180I","pop_ac_use"]
forest_data$ac_control[5] = exac_pop_c129_summary["V210I","pop_ac"]
forest_data$ac_control[6] = sum(ac_exac$ac[ac_exac$variant %in% mendelians]) # 0
forest_data$ac_control[7] = ac_23andme[ac_23andme$display_grouping=="mendelian","group_ac_use"][1]
forest_data$n_control = 0
forest_data$n_control[1] = exac_pop_c129_summary["M232R","pop_n"]
forest_data$n_control[2] = ac_23andme["M232R","pop_n"]
forest_data$n_control[3] = exac_pop_c129_summary["V180I","pop_n"]
forest_data$n_control[4] = ac_23andme["V180I","pop_n"]
forest_data$n_control[5] = exac_pop_summary$n_exac[exac_pop_summary$population_code=="TSI"]
forest_data$n_control[6] = sum(exac_pop_summary$n_exac)
forest_data$n_control[7] = round(mean(ac_23andme$called_genotypes[ac_23andme$variant %in% mendelians]))
# compute penetrance estimates and confidence intervals
forest_data$lower95 = 0
forest_data$best = 0
forest_data$upper95 = 0
for (i in 1:dim(forest_data)[1]) {
  forest_data[i,c("lower95","best","upper95")] = do.call(penetrance_confint,as.list(as.integer(forest_data[i,3:6])))
}
forest_data$case_af = forest_data$ac_case / (2*forest_data$n_case)
forest_data$control_af = forest_data$ac_control / (2*forest_data$n_control)
forest_data$case_af_display = paste(formatC(100*forest_data$case_af,format='fg',digits=2),'%',sep='')
forest_data$control_af_display = paste(formatC(100*forest_data$control_af,format='fg',digits=2),'%',sep='')
forest_data$control_af_display[c(4,7)] = paste("<",forest_data$control_af_display[c(4,7)],"*",sep='')
forest_data$case_display = 'Cases'
forest_data$control_display = ''
forest_data$control_display[c(1,3,5,6)] = 'ExAC'
forest_data$control_display[c(2,4,7)] = '23andMe'
forest_data$comparison_text = paste(forest_data$case_display," (",
                                    forest_data$case_af_display,") vs. ",
                                    forest_data$control_display," (",
                                    forest_data$control_af_display,")",sep='')

cols = c(-5.2,-4.2,-3.7,0,4.5)
names(cols) = c("Variant(s)","Ancestry","Comparison","Forest","Bar")
headers = c("Variant(s)","Ancestry","Comparison (allele frequencies)","Lifetime risk (95%CI)","Positive family history in cases")

# cairo_pdf() is better than pdf() at handling non-ASCII utf8 characters correctly
# see http://stackoverflow.com/a/12775087/3806692
# note that saving this plot as PDF takes ~30 seconds; it doesn't mean the script
# has hung.
cairo_pdf('figures/figure3.pdf',width=15,height=6,pointsize=15,bg='#FFFFFF00') # default pointsize = 12. background transparent.
par(mar=c(5,2,3,2))
xlims = c(-5.7,8.5)
vertical_lines = c(log10(2e-4) - log10(1e-4), 1:4)
plot(NA,NA,xlim=xlims,ylim=c(0,4),axes=FALSE,xlab='',ylab='',xaxs='i',yaxs='i')

# plot text at left
mtext(side=3,at=cols+c(-.5,-.5,0,.3,0)+.05,text=headers,font=2,cex=0.9,adj=0)
text(x=rep(cols["Variant(s)"],4),y=c(3.5,2.5,1.5,0.5),font=2,cex=.9,
     labels=c("M232R","V180I","V210I","P102L\nA117V\nD178N\nE200K"))
text(x=rep(cols["Ancestry"],4),y=c(3.5,2.5,1.5,0.5),font=2,cex=.9,
     labels=c("Japanese","Japanese","Italians","Global"))
text(x=rep(cols["Comparison"],6),y=forest_data$yvals,pos=4,font=2,cex=.9,
     labels=forest_data$comparison_text)
# plot the forest plot
x_offset = -log10(1e-4) # plot starts at .01%, which is now below the baseline assumed risk
points(x=x_offset+log10(cols["Forest"]+forest_data$best),y=forest_data$yval,pch=19,cex=1.8)
segments(x0=x_offset+log10(cols["Forest"]+forest_data$lower95),x1=x_offset+log10(cols["Forest"]+forest_data$upper95), y0=forest_data$yval, y1=forest_data$yval, lwd=4)
# for (i in 1:dim(forest_data)[1]) {
#   points(x=x_offset+log10(cols["Forest"]+forest_data$best[i]),y=forest_data$yval[i],pch=19,cex=2)
#   points(x=x_offset+log10(cols["Forest"]+c(forest_data$lower95[i],forest_data$upper95[i])),y=rep(forest_data$yval[i],2),type='l',lwd=4)
# }
abline(v=vertical_lines,lwd=1)
axis(side=1,at=vertical_lines,labels=c(".02%",".1%","1%","10%","100%"),lwd=0,lwd.ticks=0,font=2)
mtext(side=1,at=vertical_lines[c(1,5)],text=c('population\nbaseline risk','complete\npenetrance'),line=2,font=2,padj=1)

# plot family history barplot
famhx_data = data.frame(type=c('M232R','V180I','V210I','E200K','GSS','FFI'),
                        hx_positive=c(2,5,7,56,23,44),
                        out_of=c(63,218,57,114,33,50),
                        y=c(3.5,2.5,1.5,0.75,0.50,0.25),
                        text=c("3%","2%","12%","49% (E200K)","70% (GSS\u2020)","88% (FFI\u2021)"), # \u2020 is &dagger;. \u2021 is &Dagger;
                        color=rep('#CD2626',6))
confint_results = binom.confint(x=famhx_data$hx_positive, n=famhx_data$out_of, method='wilson')
famhx_data$x = confint_results$mean
famhx_data$lower95 = confint_results$lower
famhx_data$upper95 = confint_results$upper

segments(x0=cols["Bar"], x1=cols["Bar"]+famhx_data$x*3, y0=famhx_data$y, y1=famhx_data$y,lend=1,lwd=25,col=famhx_data$color)
arrows(x0=cols["Bar"]+famhx_data$lower95*3, x1=cols["Bar"]+famhx_data$upper95*3, y0=famhx_data$y, y1=famhx_data$y, code=3, length=.05, angle=90, col='#000000', lwd=1.5)
text(x=cols["Bar"]+famhx_data$upper95*3, y=famhx_data$y, pos=4, label=famhx_data$text, font=2, cex=.9)

abline(h=0:4,lwd=.5)

abline(v=cols["Bar"],lwd=4,col='black')
abline(h=c(0,4),lwd=4,col='black')
abline(v=xlims,lwd=4,col='black')

dev.off()


pdf('figures/figure4.pdf',width=12,height=4)
par(mar=c(4,4,4,4))
lofs = data.frame(codons=c(20,37,75,131,145,160,163,178,186,226,227),
                  names=c("G20Gfs84X","R37X","Q75X","G131X","Y145X","Q160X","Y163X","D178Efs25X","Q186X","Y226X","Q227X"),
                  yvals=c(.8,.8,.8,.8,.7,.65,.85,.8,.7,.8,.8),
                  pos=c(3,3,3,3,3,2,2,3,4,2,4),
                  seen_in=c("healthy control","healthy control","healthy control","unknown control",rep("prion disease patient",7)))
cleaved_color = '#444444'
prp_color = '#EDA23C'
disp_color = '#108F5F'
plot(NA,NA,xlim=c(1,254),ylim=c(-.3,1),axes=FALSE,xlab='',ylab='')
points(lofs$codons,lofs$yvals,type='h',lwd=2,lend=1)
points(lofs$codons,lofs$yvals,pch=20)
text(x=lofs$codons,y=lofs$yvals,pos=lofs$pos,labels=lofs$names,cex=.8,font=2)
points(1:253,rep(0,253),col=cleaved_color,type='l',lwd=25,lend=1)
points(23:230,rep(0,208),col=prp_color,type='l',lwd=25,lend=1)
points(23:94,rep(0,94-23+1),col=disp_color,type='l',lwd=25,lend=1)
axis(side=1,at=c(1,(1:4)*50,253),labels=c(1,(1:4)*50,253),lwd=NA,lwd.ticks=1,cex.axis=.8)
mtext(side=1,text=expression(italic('PRNP')~'codon number'),line=2)
points(c(15,80),c(1,1),type='l',lwd=5,lend=1)
mtext(side=3,at=(23+75)/2,text='seen in\nhealthy controls',cex=1.5)
points(c(140,230),c(1,1),type='l',lwd=5,lend=1)
mtext(side=3,at=(140+230)/2,text='seen in\nprion disease patients',cex=1.5)
points(c(128,133),c(1,1),type='l',lwd=5,lend=1)
mtext(side=3,at=131,adj=1,text='unknown phenotype',cex=.8)
#text(131,.95,labels='unknown phenotype',cex=.8,pos=2)
text((1+23)/2,-.1,pos=1,labels='signal peptide',col=cleaved_color,cex=.8,font=2)
text((23+93)/2,-.1,pos=1,labels='deleted in susceptible mice',col=disp_color,cex=.8,font=2)
text((231+253)/2,-.1,pos=1,labels='GPI signal',col=cleaved_color,cex=.8,font=2)
dev.off()



#### Calcluations for the descriptive statistics quoted in the main text

# "In our surveillance cohorts, 65% of cases underwent PRNP open reading frame sequencing, with 12% of all cases, or 18% of sequenced cases, possessing a rare variant"
ac_case["TOTAL","prnp_sequenced"]/ac_case["TOTAL","definite_plus_probable"]
ac_case["TOTAL","with_rare_variants"]/ac_case["TOTAL","definite_plus_probable"]
ac_case["TOTAL","with_rare_variants"]/ac_case["TOTAL","prnp_sequenced"]

# "0.4% of ExAC individuals harbor a rare (<0.1%) missense variant"
sum(ac_exac$ac[(ac_exac$class=='missense')]) / sum(exac_pop_summary$n_exac)

# "...four missense variants - P102L, A117V, D178N and E200K... account for >50% of genetic prion disease cases..."
sum(ac_case["TOTAL",mendelians])/ac_case["TOTAL","with_rare_variants"]

# "...with confidence intervals including all reported estimates of E200K penetrance based on survival analysis, which range from ~60% to ~90%"
forest_data[forest_data$comparisons=='Mendelians_ExAC',c('lower95','upper95')]

# "Rounding down to 1 instead would raise our estimates of penetrance..." (Materials and Methods)
penetrance_confint(ac_case["Japan","V180I"] + ac_case["Japan","V180I_and_M232R_in_trans"],
                   ac_case["Japan","prnp_sequenced"],
                   1,
                   ac_23andme["V180I","pop_n"])

penetrance_confint(ac_case["Japan","V180I"] + ac_case["Japan","V180I_and_M232R_in_trans"],
                   ac_case["Japan","prnp_sequenced"],
                   5,
                   ac_23andme["V180I","pop_n"])

penetrance_confint(sum(ac_case["TOTAL",mendelians]),
                   ac_case["TOTAL","prnp_sequenced"],
                   1,
                   round(mean(ac_23andme$called_genotypes[ac_23andme$variant %in% mendelians])))

penetrance_confint(sum(ac_case["TOTAL",mendelians]),
                   ac_case["TOTAL","prnp_sequenced"],
                   5,
                   round(mean(ac_23andme$called_genotypes[ac_23andme$variant %in% mendelians])))

# "Other variants" paragraph in the Supplementary Discussion
discussed = c(mendelians,"M232R","V180I","V210I","P39L","R148H","R208C","E196A","V203I","T188R","R208H")
missense_cases_discussed = sum(ac_case["TOTAL",colnames(ac_case) %in% discussed])
missense_cases_not_discussed = ac_case["TOTAL",(1:dim(ac_case)[2] > 19 & !(colnames(ac_case) %in% discussed) & nchar(colnames(ac_case)) <= 5 & !grepl('X',colnames(ac_case)))]
sum(missense_cases_not_discussed)

controls_with = sum(ac_exac$ac[(ac_exac$class=='missense' & !(ac_exac$variant %in% discussed))])
controls_without = sum(exac_pop_summary$n_exac) - sum(ac_exac$ac[(ac_exac$class=='missense' & !(ac_exac$variant %in% discussed))])
cases_with = sum(missense_cases_not_discussed)
cases_without = ac_case["TOTAL","prnp_sequenced"] - sum(missense_cases_not_discussed)
fisher.test(matrix(c(cases_with, cases_without, controls_with, controls_without),nrow=2),alternative='two.sided')

