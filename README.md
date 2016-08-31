# amplican: R package for the analysis of CRISPR induced mutations using sequencing amplicon data

#### This package is still under development although this version is stable and can already be used.

#### About

Our R package’s main function performs alignment of the amplicon reads, calculates multiple statistics (eg. cut rates, frameshifts) and presents results in form of html pages. Whole analysis starting from fastq files, ending with html reports is compacted to one simple function. Each report contains plots and tables in high quality format, that can be immediately further used by users. We provide a wide range of reports that allow breakdown of experiments by different constraints. We support experiment, barcode, group/user defined, guide and amplicon as possible levels of statistics aggregation. 


#### Installation

When final version of the package will be realeasted, it will be available on Bioconductor. For this version dependencies must be installed by the users:

- When on windows, make sure you have latest Rtools. When installing Rtools, it is sufficient to choose the “Package authoring installation” option. Also during the installation, you must tick the “edit system PATH” box.
- Dependencies from CRAN:  
```r
install.packages(c("Rcpp", "R.utils", "doParallel", "foreach", "ggplot2", "stringr", "rmarkdown", "knitr", "devtools"))
```  
- Dependencies from Bioconductor:  
```r
source("https://bioconductor.org/biocLite.R")  
biocLite(c("seqinr", "ShortRead", "IRanges", "GenomicRanges", "S4Vectors", "ggbio"))
```  
- Install amplican using devtools:  
```r
devtools::install_bitbucket("valenlab/amplican", build_vignettes = TRUE)
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
?amplican_plot_cuts

# read vignette
browseVignettes("amplican")
```

#### Feedback

Please feel free to provide feedback or wanted functionality. My contact adress is Kornel.Labun at uib.no.