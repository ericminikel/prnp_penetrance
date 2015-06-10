### Quantifying penetrance in a dominant disease gene using large population control cohorts

This repository holds the source code and data for a manuscript in progress. The four figures for the main text can be re-generated solely from the data stored in the supplementary tables. First, make sure you are running R 3.1.2 or later and have the [`sqldf`](http://cran.r-project.org/web/packages/sqldf/index.html) and [`binom`](http://cran.r-project.org/web/packages/binom/index.html) packages installed for R:

```r
install.packages("sqldf")
install.packages("binom")
```

Then, on the command line, clone this repo, delete the figures and re-run the [generate_figures.r](/src/generate_figures.r) script:

```bash
git clone git@github.com:ericminikel/prnp_penetrance.git
rm figures/figure1b.pdf
rm figures/figure2.pdf
rm figures/figure3.pdf
rm figures/figure4.pdf
Rscript src/generate_figures.r
```

Tab-delimited text versions of all supplementary tables are available in [supplement](/supplement), and IGV screenshots for all rare *PRNP* variants deemed to be genuine are in the [igv](/supplement/igv) directory.