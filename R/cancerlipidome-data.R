#' Levels of lipids in benign and malignant breast tumors in humans.
#'
#' This data set contains levels of 409 named lipids in 118 human breast tumor
#'    tissue samples.
#'
#' @docType data
#' @format A long-format data frame with 48262 rows and 7 variables:
#' \describe{
#'   \item{ID}{Participant number}
#'   \item{Group}{Diagnosis of the type tumor: benign, cancer, or metastasis}
#'   \item{Race}{Ethnic background of the participant}
#'   \item{Stage}{Diagnosis of the stage of the tumor}
#'   \item{Type}{Sub-type of the breast tumor. IDC: Invasive Ductal Carcinoma}
#'   \item{Lipid_Name}{Name of the lipid. The names are in
#'         the format 'XY(C:D)', where 'XY' is the abbreviation of the lipid
#'         class, 'C' is the total number of carbon atoms in the fatty-acid
#'         chains, and 'D' is the total number of double-bonds in the fatty
#'         acid chains.}
#'   \item{Lipid_Level}{Measured level of the lipid.}
#' }
#' @keywords data datasets human lipidome lipids lipidomics breast cancer
#'    tissue tumor molecule invasive ductal carcinoma diagnosis
#' @references Purwaha, P., et al.
#'    Unbiased lipidomic profiling of triple-negative breast cancer tissues
#'    reveals the association of sphingomyelin levels with patient disease-free
#'    survival.
#'    Metabolites 8, 41 (2018)
#'    (\href{https://doi.org/10.3390/metabo8030041}{doi: 10.3390/metabo8030041})
#' @source This data is available at the NIH Common Fund's National
#'    Metabolomics Data Repository (NMDR) website, the Metabolomics Workbench,
#'    \url{https://www.metabolomicsworkbench.org},
#'    where it has been assigned Project ID PR000742.
#'    The data can be accessed directly via its Project DOI:
#'    \href{http://dx.doi.org/10.21228/M8RX01}{10.21228/M8RX01}.
#'    This work was supported by NIH grant, U2C- DK119886.
#' @usage
#' data( cancerlipidome )
#' @examples
#' # Import the data set.
#' data( cancerlipidome )
#' \dontshow{
#' # Reduce the size of the data set for automated checking.
#' tmp <-
#'     tapply(
#'         X = cancerlipidome$"Lipid_Level",
#'         INDEX = cancerlipidome$"Lipid_Name",
#'         FUN = function( x ){ sd( x ) / mean( x ) }
#'     )
#' tmp2 <- order( tmp, decreasing = TRUE )[ 1:20 ]
#' tmp <- names( tmp )[ tmp2 ]
#' cancerlipidome <-
#'     cancerlipidome[ cancerlipidome$"Lipid_Name" %in% tmp, ]
#' }
#' # Convert the data into wide format, where each lipid is one column and
#' # each sample is one row.
#' cancerlipidome.wide <-
#'    tidyr::pivot_wider(
#'        data = cancerlipidome,
#'        names_from = Lipid_Name,
#'        values_from = Lipid_Level
#'    )
#' # Inspect the data frame.
#' # View( cancerlipidome.wide )
#' # Create a mapping of the lipid names.
#' names.mapping <-
#'    map_lipid_names( x = unique( cancerlipidome$"Lipid_Name" ) )
#' # Compute the regression models.
#' result.limma <-
#'    compute_models_with_limma(
#'        x = cancerlipidome.wide,
#'        dependent.variables = names.mapping$"Name",
#'        independent.variables = c( "Group" )
#'    )
#' \donttest{
#' # Create a figure of all lipids and factors.
#' figure.output <-
#'   heatmap_lipidome_from_limma(
#'     x = result.limma$"model",
#'     names.mapping = names.mapping,
#'     axis.x.carbons = FALSE,
#'     class.facet = "row",
#'     plot.all = TRUE,
#'     plot.individual = FALSE,
#'     print.figure = TRUE,
#'     scales = "free",
#'     space = "free"
#'   )
#' }
#' # Create individual figures for each factor.
#' figure.output <-
#'    heatmap_lipidome_from_limma(
#'        x = result.limma$"model",
#'        names.mapping = names.mapping,
#'        axis.x.carbons = FALSE,
#'        class.facet = "wrap",
#'        omit.class = "PA",
#'        plot.all = FALSE,
#'        plot.individual = TRUE,
#'        print.figure = FALSE,
#'        scales = "free",
#'        space = "free"
#'    )
#' # Print the figure of differences between cancer and benign tumors.
#' print( figure.output[[ "GroupCancer" ]] )
#'
"cancerlipidome"
