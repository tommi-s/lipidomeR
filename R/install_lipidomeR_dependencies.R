#' Install dependencies for the lipidomeR.
#' @param pkgs.bioconductor Character vector of Bioconductor packages to install.
#' @param ... Other arguments to the \code{BiocManager::install()} function.
install_lipidomeR_dependencies <-
  function(
    pkgs.bioconductor = c( "limma" ),
    ...
  ) {

    BiocManager::install(
      pkgs = pkgs.bioconductor,
      ...
    )

    return( NULL )

}
