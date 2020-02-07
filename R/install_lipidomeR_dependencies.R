#' Install the Bioconductor dependencies for the 'lipidomeR'
#'
#' Run this function, if you encounter an error with missing packages,
#'    such as, \href{https://doi.org/doi:10.18129/B9.bioc.limma}{limma}.
#'
#' @param pkgs.bioconductor Character vector of Bioconductor packages to
#'    install.
#'
#' @export
#'
install_Bioconductor_dependencies <-
  function(
    pkgs.bioconductor = c( "limma" )
  ) {

    is.legacy.R <-
      ( as.numeric( R.version$"major" ) < 3 )

    is.legacy.R <-
      is.legacy.R |
      (
        ( as.numeric( R.version$"major" ) == 3 ) &
          ( as.numeric( R.version$"minor" ) < 5.0 )
      )

    if ( is.legacy.R ) {

      tmp <-
        paste(
          "Legacy R with no BiocManager package available.",
          "Install Bioconductor dependencies",
          "by manually running the following lines:\n"
        )

      tmp <-
        paste(
          tmp,
          "source( 'https://bioconductor.org/biocLite.R'"
        )

      tmp <-
        paste(
          tmp,
          "BiocInstaller::biocLite( c(",
          paste(
            pkgs.bioconductor,
            collapse = "," ),
          ")"
        )

      stop( tmp )

    } else {

      BiocManager::install( pkgs = pkgs.bioconductor )

    }

    return( NULL )

}
