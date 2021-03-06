### Quantifying prion disease penetrance using large population control cohorts

This repository holds the source code, data, and pre-print for the following manuscript:

[Minikel et al. Quantifying prion disease penetrance using large population control cohorts. Sci. Transl. Med. 8, 322ra9 (2016). DOI: 10.1126/scitranslmed.aad5169](http://stm.sciencemag.org/content/8/322/322ra9.full)

#### What's here

In this repository you can:

+ Read the pre-print version of the paper in Markdown format: [manuscript.md](manuscript.md).
+ Download tab-delimited text (.tsv) versions of all supplementary tables: [supplement](/supplement).
+ Browse IGV screenshots for all rare *PRNP* variants deemed to be genuine: [igv](/supplement/igv).
+ Re-produce the main figures in the manuscript by analyzing data from the supplementary tables in R, per the following instructions.  

#### Instructions for running code

The four figures for the main text can be re-generated solely from the data stored in the supplementary tables. First, make sure you are running R 3.1.2 or later and have the [`sqldf`](http://cran.r-project.org/web/packages/sqldf/index.html) and [`binom`](http://cran.r-project.org/web/packages/binom/index.html) packages installed for R:

```r
install.packages("sqldf")
install.packages("binom")
```

Then, on the command line, clone this repo, delete the figures and re-run the [generate_figures.r](/src/generate_figures.r) script:

```bash
git clone git@github.com:ericminikel/prnp_penetrance.git
cd prnp_penetrance
rm figures/figure1b.pdf
rm figures/figure2.pdf
rm figures/figure3.pdf
rm figures/figure4.pdf
Rscript src/generate_figures.r
```

#### See also

+ [Robert Green's perspective piece](http://stm.sciencemag.org/content/8/322/322fs3.full) about this work, including how it changed one patient's prognosis.
+ [My CureFFI.org blog post](http://www.cureffi.org/2016/01/20/does-this-mean-ill-definitely-get-the-disease/) telling the personal story behind this study.
+ Coverage in [STAT](http://www.statnews.com/2016/01/20/prion-disease-genes/), [AAAS](http://www.aaas.org/news/genetic-risk-prion-disease-reevaluated-thanks-big-data), [Nature Reviews Genetics](http://www.nature.com/nrg/journal/vaop/ncurrent/full/nrg.2016.9.html), [Neurology Today](http://journals.lww.com/neurotodayonline/_layouts/15/oaks.journals.mobile/post.aspx?blogId=1&postId=536), and [El Mundo](http://www.elmundo.es/salud/2016/01/21/569fdbff22601d10308b465f.html).

