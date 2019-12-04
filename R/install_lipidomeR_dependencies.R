#' Install the Bioconductor dependencies for the lipidomeR
#'
#' Run this function, if you encounter an error with missing packages,
#'    such as, \code{limma}.
#'
#' @param pkgs.bioconductor Character vector of Bioconductor packages to
#'    install.
#' @params ... Other arguments to the \code{BiocManager::install()} function.
#'
#' @inheritParams BiocManager::install
#'
#' @export
#'
install_Bioconductor_dependencies <-
  function(
    pkgs.bioconductor = c( "limma" ),
    ...
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

      print( "WARNING: Legacy R with no BiocManager package available." )
      print( "Install Bioconductor dependencies" )
      print( "by manually running the following lines:" )

      print( "source( 'https://bioconductor.org/biocLite.R' )" )

      print(
        paste(
          "BiocInstaller::biocLite( c(",
          paste(
            pkgs.bioconductor,
            collapse = "," ),
          ")"
        )
      )

    } else {

      BiocManager::install(
        pkgs = pkgs.bioconductor,
        ...
      )

    }

    return( NULL )

}
