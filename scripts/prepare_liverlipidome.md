Load and Prepare the liverlipidome Data Set
================
Tommi Suvitaival, Steno Diabetes Center Copenhagen,
<TSUV0001@RegionH.DK>
2020-01-23

# Introduction

This document describes, how the lipidomer::liverlipidome data set was
prepared.

The data come from the project PR000633 in the Metabolomics Workbench.

## Citation of the Data

### Publication

Gorden, D. Lee, et al. Biomarkers of NAFLD Progression: a Lipidomics
Approach to an Epidemic. J Lip Res. 56(3) 722-36 (2015).
<https://dx.doi.org/10.1194/jlr.P056002>

### Repository

This data is available at the NIH Common Fund’s National Metabolomics
Data Repository (NMDR) website, the Metabolomics Workbench,
<https://www.metabolomicsworkbench.org>, where it has been assigned
Project ID PR000633. The data can be accessed directly via it’s Project
DOI: 10.21228/M8V961 This work is supported by NIH grant, U2C-DK119886.

URL: <http://dx.doi.org/10.21228/M8V961>

# Download the Data

  - Download the data set in mzTab format from the study ST000915 of the
    project PR000633 in the Metabolomics Workbench.

## Specify the File Paths

``` r
url.download <-
  c(
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000915&ANALYSIS_ID=AN001485&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000915&ANALYSIS_ID=AN001486&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000915&ANALYSIS_ID=AN001487&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000915&ANALYSIS_ID=AN001488&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000915&ANALYSIS_ID=AN001489&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000915&ANALYSIS_ID=AN001490&MODE=d"
    
  )
```

## Download the Data

``` r
sample.info.loaded <-
  read.table(
    file = url.download[ 1 ],
    check.names = FALSE,
    header = FALSE,
    nrows = 158-69-1,
    # sep = "\t",
    skip = 69,
    stringsAsFactors = FALSE
  )

data.loaded <- list()

data.loaded[[ 1 ]] <-
  read.table(
    file = url.download[ 1 ],
    check.names = FALSE,
    header = TRUE,
    nrows = 373-274-1,
    sep = "\t",
    skip = 274,
    stringsAsFactors = FALSE
  )

data.loaded[[ 2 ]] <-
  read.table(
    file = url.download[ 2 ],
    check.names = FALSE,
    header = TRUE,
    na.strings = "ND",
    nrows = 302-261-1,
    sep = "\t",
    skip = 261,
    stringsAsFactors = FALSE
  )

data.loaded[[ 3 ]] <-
  read.table(
    file = url.download[ 3 ],
    check.names = FALSE,
    header = TRUE,
    nrows = 288-276-1,
    sep = "\t",
    skip = 276,
    stringsAsFactors = FALSE
  )

data.loaded[[ 4 ]] <-
  read.table(
    file = url.download[ 4 ],
    check.names = FALSE,
    header = TRUE,
    nrows = 456-265-1,
    sep = "\t",
    skip = 265,
    stringsAsFactors = FALSE
  )

data.loaded[[ 5 ]] <-
  read.table(
    file = url.download[ 5 ],
    check.names = FALSE,
    header = TRUE,
    nrows = 368-269-1,
    sep = "\t",
    skip = 269,
    stringsAsFactors = FALSE
  )

data.loaded[[ 6 ]] <-
  read.table(
    file = url.download[ 6 ],
    check.names = FALSE,
    header = TRUE,
    nrows = 373-261-1,
    sep = "\t",
    skip = 261,
    stringsAsFactors = FALSE
  )
```

# Prepare the Data

  - Reformat the data sets
  - Merge lipidomics data from multiple measurements
  - Reformat the lipid names
  - Round the lipidomics observations
  - Merge lipidomics data with sample information
  - Omit lipids with less than two observations per group
  - Save the resulting data set

## Reformat

### Sample Information

``` r
sample.info <- sample.info.loaded

for ( i in 4:ncol( sample.info ) ) {
  
  tmp <-
    stringr::str_split_fixed(
      string = sample.info[ , i ],
      pattern = "(\\:)|(\\=)",
      n = 2
    )
  
  sample.info[ , i ] <- tmp[ , 2 ]
  
  colnames( sample.info )[ i ] <- tmp[ 1, 1 ]
  
  sample.info[ , i ] <-
    stringr::str_replace_all(
      string = sample.info[ , i ],
      pattern = "\\;",
      replacement = ""
    )
  
}

sample.info[ sample.info == "-"] <- NA

colnames( sample.info )[ colnames( sample.info ) == "V3" ] <- "ID"

sample.info <- sample.info[ , 3:ncol( sample.info ) ]

sample.info$"Diagnosis" <-
  factor(
    x = sample.info$"Diagnosis",
    levels = c( "Normal", "Steatosis", "NASH", "Cirrhosis" )
  )

names.numeric <- c( "BMI", "AGE", "AST", "ALT", "ALKP", "TBIL", "GLUCOSE" )

for ( i in 1:length( names.numeric ) ) {
  
  sample.info[[ names.numeric[ i ] ]] <-
    as.numeric( sample.info[[ names.numeric[ i ] ]] )

}
```

### Lipids

``` r
data <- data.loaded

data <-
  lapply(
    X = data,
    FUN =
      function( x ) {
        
        rownames( x ) <- x$"Samples"
        
        x <- t( x[ -1, -1 ] )
        
        y <-
          apply(
            X = x,
            MAR = 2,
            FUN = as.numeric
          )
        
        rownames( y ) <- rownames( x )
        
        return( y )
        
      }
  )
```

    ## Warning in apply(X = x, MAR = 2, FUN = as.numeric): NAs introduced by coercion
    
    ## Warning in apply(X = x, MAR = 2, FUN = as.numeric): NAs introduced by coercion
    
    ## Warning in apply(X = x, MAR = 2, FUN = as.numeric): NAs introduced by coercion

## Merge Lipid Data

``` r
data.merged <- data[[ 1 ]]

for ( i in 2:length( data ) ) {

  data.merged <-
    merge(
      x = data.merged,
      y = data[[ i ]],
      by = "row.names"
    )

  rownames( data.merged ) <- data.merged$"Row.names"

  data.merged$"Row.names" <- NULL

}
 
data <- data.merged
```

## Reformat Lipid Names

``` r
is.C.name <-
  stringr::str_detect(
    string = colnames( data ),
    pattern = "C[0-9]+" # "C[0-9]+:?[0-9]+?D?H?"
  )

if ( any( is.C.name ) ) {
  
  colnames( data )[ is.C.name ] <-
    stringr::str_replace(
      string = colnames( data )[ is.C.name ],
      pattern = "DH ",
      replacement = " DH-"
    )
  
  name.split <-
    stringr::str_split_fixed(
      string = colnames( data )[ is.C.name ],
      pattern = " ",
      n = 3
    )
  
  number <-
    stringr::str_replace_all(
      string = name.split[ , 1 ],
      pattern = "[A-Z]+",
      replacement = ""
    )
  
  number.split <-
    stringr::str_split_fixed(
      string = number,
      pattern = ":",
      n = 2
    )
  
  number.split[ which( number.split[ , 2 ] == "" ), 2 ] <- "0"
  
  number.split[ , 1 ] <- as.numeric( number.split[ , 1 ] ) + 18
  
  name.split[ , 2 ] <-
    stringr::str_replace(
      string = name.split[ , 2 ],
      pattern = "Sphingomyelin",
      replacement = "SM"
    )
  
  colnames( data )[ is.C.name ] <-
    paste0(
      name.split[ , 2 ],
      "(",
      number.split[ , 1 ],
      ":",
      number.split[ , 2 ],
      ")"
    )
  
}

is.numbered <-
  grepl(
    x = colnames( data ),
    pattern = "[A-Za-z]\\s?\\("
  )

data <- data[ , is.numbered ]

# Fixes to names

is.N3 <-
  grepl(
    x = colnames( data ),
    pattern = " N3" )

colnames( data )[ is.N3 ] <-
  stringr::str_replace(
    string = colnames( data )[ is.N3 ],
    pattern = " N3",
    replacement = ""
  )

colnames( data )[ is.N3 ] <-
  stringr::str_replace(
    string = colnames( data )[ is.N3 ],
    pattern = "\\(",
    replacement = "-N3("
  )

is.N6 <-
  grepl(
    x = colnames( data ),
    pattern = " N6" )

colnames( data )[ is.N6 ] <-
  stringr::str_replace(
    string = colnames( data )[ is.N6 ],
    pattern = " N6",
    replacement = ""
  )

colnames( data )[ is.N6 ] <-
  stringr::str_replace(
    string = colnames( data )[ is.N6 ],
    pattern = "\\(",
    replacement = "-N6("
  )

is.N9 <-
  grepl(
    x = colnames( data ),
    pattern = " N9" )

colnames( data )[ is.N9 ] <-
  stringr::str_replace(
    string = colnames( data )[ is.N9 ],
    pattern = " N9",
    replacement = ""
  )

colnames( data )[ is.N9 ] <-
  stringr::str_replace(
    string = colnames( data )[ is.N9 ],
    pattern = "\\(",
    replacement = "-N9("
  )

tmp <-
  grep(
    x = colnames( data ),
    pattern = "p\\)"
  )

colnames( data )[ tmp ] <-
  stringr::str_replace(
    string = colnames( data )[ tmp ],
    pattern = "\\(",
    replacement = "-P("
  )

colnames( data )[ tmp ] <-
  stringr::str_replace(
    string = colnames( data )[ tmp ],
    pattern = "p\\)",
    replacement = ")"
  )

colnames( data ) <-
  stringr::str_replace(
    string = colnames( data ),
    pattern = " \\(",
    replacement = "("
  )
```

``` r
colnames( data ) <-
  stringr::str_replace(
    string = colnames( data ),
    pattern = "\\.[0-9]",
    replacement = ""
  )
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
```

## Omit Lipids with Less Than Two Observations Per Group

``` r
N.found <-
  by(
    data = !is.na( data[ , names.lipids ] ),
    INDICES = data$"Diagnosis",
    FUN = colSums
  )

N.found <- simplify2array( N.found )

is.found <- ( N.found > 1 )

is.found <-
  apply(
    X = is.found,
    MAR = 1,
    FUN = all
  ) 

names.missing <- names.lipids

names.missing <- names.missing[ !is.found ]

print( names.missing )
```

    ##  [1] "CE(14:0)"     "CE(15:1)"     "CE(15:0)"     "CE(16:2)"     "CE(16:1)"    
    ##  [6] "CE(17:1)"     "CE(17:0)"     "CE(18:0)"     "CE(20:1)"     "CE(22:1)"    
    ## [11] "CE(22:0)"     "TG(40:2)"     "TG(40:1)"     "TG(40:0)"     "TG(41:5)"    
    ## [16] "TG(41:4)"     "TG(41:3)"     "TG(41:2)"     "TG(41:1)"     "TG(41:0)"    
    ## [21] "TG(42:5)"     "TG(42:4)"     "TG(42:3)"     "TG(42:2)"     "TG(43:2)"    
    ## [26] "TG(43:1)"     "TG(43:0)"     "TG(44:5)"     "TG(44:4)"     "TG(46:5)"    
    ## [31] "TG(47:5)"     "TG(47:4)"     "TG(49:4)"     "TG(56:0)"     "TG(58:0)"    
    ## [36] "TG(59:2)"     "TG(59:1)"     "TG(59:0)"     "TG(60:2)"     "TG(60:1)"    
    ## [41] "TG(60:0)"     "1,2-DG(32:5)" "1,2-DG(37:5)" "1,2-DG(37:4)" "PA(34:0)"

``` r
data <- data[ , !( colnames( data ) %in% names.missing ) ]

names.lipids <- names.lipids[ is.found ]
```

# Convert the Data into a Longer Format

``` r
liverlipidome <- 
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
  liverlipidome,
  file = "../data/liverlipidome.rda",
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
    ## [13] digest_0.6.23    assertthat_0.2.1 tibble_2.1.3     lifecycle_0.1.0 
    ## [17] crayon_1.3.4     purrr_0.3.3      tidyr_1.0.0      vctrs_0.2.1     
    ## [21] zeallot_0.1.0    glue_1.3.1       evaluate_0.14    rmarkdown_2.1   
    ## [25] stringi_1.4.4    compiler_3.6.2   pillar_1.4.3     backports_1.1.5 
    ## [29] pkgconfig_2.0.3
