Load and Prepare the humanlipidome Data Set
================
Tommi Suvitaival, Steno Diabetes Center Copenhagen,
<TSUV0001@RegionH.DK>
2019-12-04

# Introduction

This document describes, how the lipidomer::humanlipidome data set was
prepared.

The data come from the project PR00004 in the Metabolomics Workbench.

## Citation of the Data

### Publication

Quehenberger, O. et al. Lipidomics reveals a remarkable diversity of
lipids in human plasma. J Lipid Res. 51, 3299-3305 (2010).
<https://doi.org/10.1194/jlr.M009449>

### Repository

This data is available at the NIH Common Fund’s National Metabolomics
Data Repository (NMDR) website, the Metabolomics Workbench,
<https://www.metabolomicsworkbench.org>, where it has been assigned
Project ID PR000004. The data can be accessed directly via it’s Project
DOI: 10.21228/M8MW26 This work is supported by NIH grant, U2C- DK119886.

# Download the Data

  - Download the data set in mzTab format from the project PR00004 in
    the Metabolomics Workbench.

## Specify the File Paths

``` r
url.download <-
  c(
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000006&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000011&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000010&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000007&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000005&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000004&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000009&MODE=d",
    "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?STUDY_ID=ST000004&ANALYSIS_ID=AN000008&MODE=d"
    
  )
```

## Download the Data

``` r
data.loaded <-
  lapply(
    X = url.download,
    FUN = readLines
  )
```

# Prepare the Data

  - Reshape the data into value-unit pairs
  - Convert all values to the unit umol/mL
  - Bind data from multiple measurements
  - Clean up the lipid names
  - Subset the lipids that can be enumerated by size and saturation
  - Save the resulting data set

## Reshape the Data into Value-Unit Pairs

``` r
data <-
  lapply(
    X = data.loaded,
    FUN =
      function( x ) {
        
        unit <- 
          x[
            grep(
              x = x,
              pattern = "MS_METABOLITE_DATA:UNITS"
            )
            ]
        
        unit <- stringr::str_replace(
          string = unit,
          pattern = "MS_METABOLITE_DATA:UNITS",
          replacement = ""
        )
        
        unit <- stringr::str_trim( string = unit )
        
        line.start <- which( x == "MS_METABOLITE_DATA_START" ) + 3
        line.end <- which( x == "MS_METABOLITE_DATA_END" ) - 1
        
        x.data <- x[ line.start : line.end ]
        
        x.data <-
          stringr::str_replace(
            string = x.data,
            pattern = "\\( ",
            replacement = "("
          )
        
        location.value <-
          stringr::str_locate(
            string = x.data,
            pattern = "\\s[0-9]+\\.?([0-9]+)?"
          )
        
        y <-
          data.frame(
            Name =
              stringr::str_sub(
                string = x.data,
                start = 1,
                end = location.value[ , "start" ] - 1
              ),
            Value =
              stringr::str_sub(
                string = x.data,
                start = location.value[ , "start" ] + 1,
                end = location.value[ , "end" ]
              ),
            stringsAsFactors = FALSE
          )
        
        y$"Value" <- as.numeric( y$"Value" )
        
        y <- list( data = y, unit = unit )
        
        return( y )
        
      }
  )
```

## Convert All Values to the Same Unit of umol/mL

``` r
data <-
  lapply(
    X = data,
    FUN =
      function( x ) {
        
        y <- x$"data"
        
        if ( x$"unit" == "pmol/ml" ) {
          
          y$"Value" <- y$"Value" * 0.000001
          
        } else if ( x$"unit" == "nmol/ml" ) {
          
          y$"Value" <- y$"Value" * 0.001
          
        }
        
        return( y )
        
      }
  )
```

## Bind the Data from Multiple Measurements

``` r
data <- dplyr::bind_rows( data )

data <- data[ !is.na( data$"Name" ), ]
```

## Clean Up the Lipid Names

``` r
is.C.name <-
  stringr::str_detect(
    string = data$"Name",
    pattern = "C[0-9]+"
  )

if ( any( is.C.name ) ) {
  
  name.split <-
    stringr::str_split_fixed(
      string = data[ is.C.name, "Name" ],
      pattern = " ",
      n = 2
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
  
  name.split[ , 2 ] <-
    stringr::str_replace(
      string = name.split[ , 2 ],
      pattern = "Sphingomyelin",
      replacement = "SM"
    )
  
  data[ is.C.name, "Name" ] <-
    paste0(
      name.split[ , 2 ],
      "(",
      number.split[ , 1 ],
      ":",
      number.split[ , 2 ],
      ")"
    )
  
}
```

## Subset the Lipids That Can Be Enumerated to Size and Saturation

``` r
humanlipidome <-
  data[ grepl( x = data$"Name", pattern = "[A-Za-z]+\\(.+\\:.+\\)" ), ]

colnames( humanlipidome )[ colnames( humanlipidome ) == "Value" ] <- 
  "Concentration"
```

# Save the Resulting Data

``` r
save(
  humanlipidome,
  file = "../data/humanlipidome.rda",
  compress = "xz"
)
```

# SessionInfo

``` r
utils::sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18362)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3       crayon_1.3.4     digest_0.6.23    dplyr_0.8.3     
    ##  [5] assertthat_0.2.1 R6_2.4.1         magrittr_1.5     evaluate_0.14   
    ##  [9] pillar_1.4.2     rlang_0.4.2      stringi_1.4.3    rmarkdown_1.18  
    ## [13] tools_3.6.1      stringr_1.4.0    glue_1.3.1       purrr_0.3.3     
    ## [17] xfun_0.11        yaml_2.2.0       compiler_3.6.1   pkgconfig_2.0.3 
    ## [21] htmltools_0.4.0  tidyselect_0.2.5 knitr_1.26       tibble_2.1.3
