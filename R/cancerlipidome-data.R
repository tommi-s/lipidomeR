#' Levels of lipids in benign and malignant breast tumors in humans.
#'
#' This data set contains levels of 409 named lipids in human breast tissue.
#'
#' @docType data
#' @format A data frame with 118 rows and 414 variables:
#' \describe{
#'   \item{ID}{Participant number}
#'   \item{Group}{Diagnosis of the type tumor: benign, cancer, or metastasis}
#'   \item{Race}{Ethnic background of the participant}
#'   \item{Stage}{Diagnosis of the stage of the tumor}
#'   \item{Type}{Subtype of the breast tumor. IDC: Invasive Ductal Carcinoma}
#'   \item{CE(16:0) ... TG(62:7)}{Levels of 409 lipids. The names are in
#'         the format XY(C:D), where XY is the abbreviation of the lipid class,
#'         C is the total number of carbon atoms in the fatty-acid chains, and
#'         D is the total number of double-bonds in the fatty acid chains.}
#' }
#' @keywords data datasets human lipidome lipids lipidomics breast cancer
#'    tissue tumor molecule invasive ductal carcinoma diagnosis
#' @references Purwaha, P., et al.
#'    Unbiased lipidomic profiling of triple-negative breast cancer tissues
#'    reveals the association of sphingomyelin levels with patient diseasefFree
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
#' @usage data( cancerlipidome )
#' @examples
#' # Load the data set.
#' data( cancerlipidome )
#' names.mapping <-
#'    map_lipid_names( x = colnames( cancerlipidome )[ -( 1:5 ) ] )
#' result.limma <-
#'    compute_models_with_limma(
#'        x = cancerlipidome,
#'        independent.variables = c( "Group" ),
#'        dependent.variables = names.mapping$"Name"
#'    )
#' figure.output <-
#'    heatmap_lipidome_from_limma(
#'        x = result.limma$"model",
#'        names.mapping = names.mapping,
#'        axis.x.carbons = FALSE,
#'        omit.class = "PA",
#'        print.figure = FALSE,
#'        class.facet = "wrap",
#'        plot.all = FALSE,
#'        plot.individual = TRUE,
#'        scales = "free",
#'        space = "free"
#'    )
#' print( figure.output$"figures"[[ 2 ]] )
#'
"cancerlipidome"
