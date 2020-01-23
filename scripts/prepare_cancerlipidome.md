Load and Prepare the cancerlipidome Data Set
================
Tommi Suvitaival, Steno Diabetes Center Copenhagen,
<TSUV0001@RegionH.DK>
2020-01-23

# Introduction

This document describes, how the lipidomer::cancerlipidome data set was
prepared.

The data come from the project PR000742 in the Metabolomics Workbench.

## Citation of the Data

### Publication

Purwaha, P., et al. Unbiased Lipidomic Profiling of Triple-Negative
Breast Cancer Tissues Reveals the Association of Sphingomyelin Levels
with Patient Disease-Free Survival. Metabolites 2018, 8, 41.
<https://doi.org/10.3390/metabo8030041>

### Repository

This data is available at the NIH Common Fund’s National Metabolomics
Data Repository (NMDR) website, the Metabolomics Workbench,
<https://www.metabolomicsworkbench.org>, where it has been assigned
Project ID PR000742. The data can be accessed directly via it’s Project
DOI: 10.21228/M8RX01 This work is supported by NIH grant, U2C-DK119886.

URL: <http://dx.doi.org/10.21228/M8RX01>

# Download the Data

  - Download the data set in mzTab format from the project PR000742 in
    the Metabolomics Workbench.

## Specify the File Paths

``` r
url.download <-
  c(
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST001111&ANALYSIS_ID=AN001805&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST001111&ANALYSIS_ID=AN001806&MODE=d"
    
  )
```

## Download the Data

``` r
sample.info.loaded <-
  read.table(
    file = url.download[ 1 ],
    check.names = FALSE,
    header = FALSE,
    nrows = 148-29-1,
    sep = "\t",
    skip = 29,
    stringsAsFactors = FALSE
  )

data.loaded.pos <-
  read.table(
    file = url.download[ 1 ],
    check.names = FALSE,
    header = TRUE,
    nrows = 618-171-2,
    sep = "\t",
    skip = 171,
    stringsAsFactors = FALSE
  )

data.loaded.neg <-
  read.table(
    file = url.download[ 2 ],
    check.names = FALSE,
    header = TRUE,
    nrows = 505-171-2,
    sep = "\t",
    skip = 171,
    stringsAsFactors = FALSE
  )
```

# Prepare the Data

  - Reformat the data sets
  - Merge lipidomics data from multiple measurements
  - Reformat the lipid names
  - Combine duplicate measurements from the same lipid into a single
    value
  - Omit classes with only one lipid
  - Merge lipidomics data with sample information
  - Convert the data into long format
  - Save the resulting data set

## Reformat

### Sample Information

``` r
sample.info <- sample.info.loaded

sample.info$"V3" <- as.character( sample.info$"V3" )

colnames( sample.info )[ colnames( sample.info ) == "V3" ] <- "ID"

sample.info$"V4" <-
  stringr::str_replace(
    string = sample.info$"V4",
    pattern = "[A-za-z]+\\:",
    replacement = ""
  )

colnames( sample.info )[ colnames( sample.info ) == "V4" ] <- "Group"

tmp <-
  stringr::str_split_fixed(
    string = sample.info$"V5",
    pattern = "; ",
    n = 3
  )

colnames( tmp ) <-
  stringr::str_extract(
    string = tmp[ 1, ],
    pattern = "[A-Za-z]+\\="
  )

colnames( tmp ) <-
  stringr::str_replace(
    string = colnames( tmp ),
    pattern = "\\=",
    replacement = ""
  )

tmp <-
  apply(
    X = tmp,
    MAR = 2,
    FUN = stringr::str_replace,
    # string = tmp,
    pattern = "[A-Za-z]+\\=",
    replacement = ""
  )

tmp <- data.frame( tmp, check.names = FALSE )

sample.info <-
  dplyr::bind_cols(
    sample.info,
    tmp
  )

sample.info[ sample.info$"Stage" == "NA", "Stage" ] <- NA
sample.info[ grepl( x = sample.info$"Type", pattern = "NA" ), "Type" ] <- NA
sample.info$"Stage" <- droplevels( sample.info$"Stage" )

sample.info <-
  sample.info[ 
    ,
    !grepl( x = colnames( sample.info ), pattern = "^V[0-9]" )
    ]

sample.info$"Race" <-
  stringr::str_replace_all(
    string = sample.info$"Race",
    pattern = "\\.",
    replacement = "-"
  )

sample.info$"Type" <-
  stringr::str_replace_all(
    string = sample.info$"Type",
    pattern = "\\.",
    replacement = "-"
  )
```

### Lipids from the Positive Ion Mode

``` r
data.pos <- data.loaded.pos

rownames( data.pos ) <- data.pos$"Samples"

data.pos <- t( data.pos[ -1, -1 ] )

mode( data.pos ) <- "numeric"
```

### Lipids from the Negative Ion Mode

``` r
data.neg <- data.loaded.neg

rownames( data.neg ) <- data.neg$"Samples"

data.neg <- t( data.neg[ -1, -1 ] )

mode( data.neg ) <- "numeric"
```

## Merge Lipid Data

``` r
data <-
  merge(
    x = data.pos,
    y = data.neg,
    by = "row.names"
  )

rownames( data ) <- data$"Row.names"

data$"Row.names" <- NULL
```

## Reformat Lipid Names

``` r
tmp <-
  stringr::str_detect(
    string = colnames( data ),
    pattern = "N\\-"
  )

colnames( data )[ tmp ] <-
  stringr::str_replace(
    string = colnames( data )[ tmp ],
    pattern = "N\\-\\(",
    replacement = "N-"
  )

colnames( data )[ tmp ] <-
  stringr::str_replace(
    string = colnames( data )[ tmp ],
    pattern = "noyl\\)",
    replacement = "noyl"
  )

colnames( data ) <-
  stringr::str_replace(
    string = colnames( data ),
    pattern = "\\([0-9].+$",
    replacement = ""
  )

tmp <-
  stringr::str_detect(
    string = colnames( data ),
    pattern = "[0-9]$"
  )

colnames( data )[ tmp ] <-
  paste0(
    colnames( data )[ tmp ],
    ")"
  )

tmp <-
  stringr::str_detect(
    string = colnames( data ),
    pattern = " [0-9]"
  )

colnames( data )[ tmp ] <-
  stringr::str_replace(
    string = colnames( data )[ tmp ],
    pattern = " ",
    replacement = "("
  )

tmp <-
  stringr::str_detect(
    string = colnames( data ),
    pattern = "[0-9] "
  )

colnames( data )[ tmp ] <-
  stringr::str_replace(
    string = colnames( data )[ tmp ],
    pattern = " ",
    replacement = ")"
  )

colnames( data ) <-
  stringr::str_replace(
    string = colnames( data ),
    pattern = "\\[.+\\]\\+",
    replacement = ""
  )
```

## Combine Duplicate Lipids

``` r
names.duplicates <- names( which( table( colnames( data ) ) > 1 ) )

sum.duplicates <-
  lapply(
    X = names.duplicates,
    FUN =
      function( x ) {
        
        idx.i <- which( colnames( data ) == x )
        
        tmp <-
          apply(
            X = data[ , idx.i ],
            MAR = 1,
            FUN = prod,
            na.rm = TRUE
          )
        
        return( tmp )
        
      }
  )

sum.duplicates <- simplify2array( x = sum.duplicates )
colnames( sum.duplicates ) <- names.duplicates

sum.duplicates <-
  data.frame(
    sum.duplicates,
    check.names = FALSE
  )

data <-
  merge(
    x = data[ , !( colnames( data ) %in% names.duplicates ) ],
    y = sum.duplicates,
    by = "row.names"
  )

rownames( data ) <- data$"Row.names"

data$"Row.names" <- NULL

data <- data[ , order( colnames( data ), decreasing = FALSE ) ]
```

## Omit Classes with Only One Lipid

``` r
data <- data[ , !grepl( x = colnames( data ), pattern = "SB N" ) ]
```

## Round Observations

``` r
data <- signif( x = data, digits = 3 )
```

## Merge Lipids and Sample Information

``` r
names.lipids <- colnames( data )

data <-
  merge(
    x = sample.info,
    y = data,
    by.x = "ID",
    by.y = "row.names"
  )

data$"Group" <- as.factor( x = data$"Group" )
```

# Convert the Data into a Longer Format

``` r
cancerlipidome <- 
  tidyr::pivot_longer(
    data = data,
    cols = names.lipids,
    names_to = "Lipid_Name",
    values_to = "Lipid_Level"
  )
```

# Save the Resulting Data

``` r
save(
  cancerlipidome,
  file = "../data/cancerlipidome.rda",
  compress = "xz"
)
```

# SessionInfo

``` r
utils::sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 17763)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3       knitr_1.27       magrittr_1.5     tidyselect_0.2.5
    ##  [5] R6_2.4.1         rlang_0.4.2      stringr_1.4.0    dplyr_0.8.3     
    ##  [9] tools_3.6.2      xfun_0.12        htmltools_0.4.0  yaml_2.2.0      
    ## [13] assertthat_0.2.1 digest_0.6.23    tibble_2.1.3     lifecycle_0.1.0 
    ## [17] crayon_1.3.4     purrr_0.3.3      tidyr_1.0.0      vctrs_0.2.1     
    ## [21] zeallot_0.1.0    glue_1.3.1       evaluate_0.14    rmarkdown_2.1   
    ## [25] stringi_1.4.4    compiler_3.6.2   pillar_1.4.3     backports_1.1.5 
    ## [29] pkgconfig_2.0.3
