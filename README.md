# amplican: R package for the analysis of CRISPR induced mutations using sequencing amplicon data

#### This package is still under development, although this version is stable and can be used already.

#### About

Our R package’s main function performs alignment of the amplicon reads, calculates multiple statistics (eg. cut rates, frameshifts) and presents results in form of HTML pages. Whole analysis starting from fastq files, ending with HTML reports is compacted to one simple function. Each report contains plots and tables in high quality format, that can be immediately further used by users. We provide a wide range of reports that allow breakdown of experiments by different constraints. We support experiment, barcode, group/user defined, guide and amplicon as possible levels of statistics aggregation. 

![Conceptual map of the amplican](vignettes/figures/amplican_conceptual_map.png)


#### Installation

Package is available on Bioconductor. It is possible to install from here, but it requires installation of many dependencies and is not recomended.

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("amplican")
```  

#### More information

After installation run:
```r
library(amplican)

# main functions
?amplican
?amplicanPipeline
?amplicanAlign
?amplicanReport
?plot_cuts
?metaplot_deletions

# read vignette
browseVignettes("amplican")
```

#### Feedback

Please feel free to provide feedback or desired functionality. My contact address is Kornel.Labun at uib.no.

#### Benchmarks

Since my publication [Labun et al. 2019, Accurate analysis of genuine CRISPR editing events with ampliCan](https://www.ncbi.nlm.nih.gov/pubmed/30850374) it seems CRISPREsso authors have updated their tool into [CRISPresso v2](https://www.nature.com/articles/s41587-019-0032-3). I have quickly rerun one of my benchmarks with their latest tool (used default settings, and installation has to work correctly as I used CRISPresso2 docker image). 

![Infuence of size of the mutation on estimating true editing efficiency](vignettes/figures/indel_rate_vs_indel_size.png)  

If you wish to replicate or comment on this small benchmark, contact me through email.
