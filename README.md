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


