#' Concentrations of lipids in a human plasma reference sample
#'
#' This data set contains concentrations of 403 named lipids in the 'National
#'    Institute of Standards and Technology Human Plasma Standard Reference
#'    Material' (NIST SRM 1950).
#'
#' @docType data
#' @format A data frame with 403 rows and 2 variables:
#' \describe{
#'   \item{Name}{Name of the lipid. The names are in
#'         the format 'XY(C:D)', where 'XY' is the abbreviation of the lipid
#'         class, 'C' is the total number of carbon atoms in the fatty-acid
#'         chains, and 'D' is the total number of double-bonds in the fatty
#'         acid chains.}
#'   \item{Concentration}{Concentration of the lipid (umol/mL)}
#' }
#' @keywords data datasets human lipidome lipids lipidomics NIST reference
#'    plasma blood molecule concentration
#' @references Quehenberger, O. et al.
#'    Lipidomics reveals a remarkable diversity of lipids in human plasma.
#'    J Lipid Res. 51, 3299-3305 (2010)
#'    (\href{https://doi.org/10.1194/jlr.M009449}{doi: 10.1194/jlr.M009449})
#' @source This data is available at the NIH Common Fund's National
#'    Metabolomics Data Repository (NMDR) website, the 'Metabolomics Workbench',
#'    \url{https://www.metabolomicsworkbench.org},
#'    where it has been assigned Project ID PR000004.
#'    The data can be accessed directly via its Project DOI:
#'    \href{http://dx.doi.org/10.21228/M8MW26}{10.21228/M8MW26}.
#'    This work was supported by NIH grant, U2C- DK119886.
#' @usage data( humanlipidome )
#' @examples
#' # Load the data set.
#' data( humanlipidome )
#' # Transform the concentrations into log-10 scale.
#' humanlipidome$"Concentration_log10_umol_per_mL" <-
#'    log10( humanlipidome$"Concentration" )
#' # Enumerate the lipid names into values.
#' names.mapping <- map_lipid_names( x = humanlipidome$"Name" )
#' # Create the lipidomeR heatmap of lipid concentrations.
#' heatmap_lipidome(
#'    x = humanlipidome[ , c( "Name", "Concentration_log10_umol_per_mL" ) ],
#'    names.mapping = names.mapping,
#'    class.facet = "wrap",
#'    x.names = "Name",
#'    fill.limits =
#'        range(
#'            x = humanlipidome$"Concentration_log10_umol_per_mL",
#'            na.rm = TRUE
#'        ),
#'    fill.midpoint =
#'        sum(
#'            range(
#'                x = humanlipidome$"Concentration_log10_umol_per_mL",
#'                na.rm = TRUE
#'            )
#'        ) / 2,
#'    melt.value.name = "Concentration_umol_per_mL_log10",
#'    scales = "free"
#' )
#'
"humanlipidome"
