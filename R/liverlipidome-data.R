#' Levels of lipids in the human liver with or without non-alcoholic liver
#'    disease (NAFLD).
#'
#' This data set contains levels of 383 named lipids in 88 liver tissue samples.
#'
#' @docType data
#' @format A long-format data frame with 33704 rows and 13 variables:
#' \describe{
#'   \item{ID}{Participant number}
#'   \item{Diagnosis}{Diagnosis of the liver: normal, steatosis, non-alcoholic
#'         steatohepatitis (NASH), or cirrhosis}
#'   \item{Gender}{Gender of the participant}
#'   \item{BMI}{Body-mass-index (BMI) of the participant}
#'   \item{Ethnicity}{Ethnicity of the participant}
#'   \item{Age}{Age of the participant}
#'   \item{AST}{Aspartate aminotransferase blood test (U/l)}
#'   \item{ALT}{Alanine aminotransferase blood test (U/l)}
#'   \item{ALKP}{Alkaline phosphatase blood test (U/l)}
#'   \item{TBIL}{Total bilirubin blood test (mg/dl)}
#'   \item{Glucose}{Glucose blood test (mg/dl)}
#'   \item{Type}{Sub-type of the breast tumor. IDC: Invasive Ductal Carcinoma}
#'   \item{Lipid_Name}{Name of the lipid. The names are in
#'         the format 'XY(C:D)', where 'XY' is the abbreviation of the lipid
#'         class, 'C' is the total number of carbon atoms in the fatty-acid
#'         chains, and 'D' is the total number of double-bonds in the fatty
#'         acid chains.}
#'   \item{Lipid_Level}{Measured level of the lipid.}
#' }
#' @keywords data datasets human lipidome lipids lipidomics non-alcoholic liver
#'    disease steatohepatitis NASH tissue molecule diagnosis
#' @references Gorden, D. Lee, et al.
#'    Biomarkers of NAFLD Progression: a Lipidomics Approach to an Epidemic.
#'    J Lip Res. 56(3) 722-36 (2015)
#'    (\href{https://dx.doi.org/10.1194/jlr.P056002}{doi: 10.1194/jlr.P056002}
#' @source This data is available at the NIH Common Fund's National
#'    Metabolomics Data Repository (NMDR) website, the 'Metabolomics Workbench',
#'    \url{https://www.metabolomicsworkbench.org},
#'    where it has been assigned Project ID PR000633.
#'    The data can be accessed directly via its Project DOI:
#'    \href{http://dx.doi.org/10.21228/M8MW26}{10.21228/M8MW26}.
#'    This work was supported by NIH grant, U2C- DK119886.
#' @usage data( liverlipidome )
#' @examples
#' # Load the data set.
#' data( liverlipidome )
#' # Convert the data into wide format, where each lipid is one column and
#' # each sample is one row.
#' liverlipidome.wide <-
#'    tidyr::pivot_wider(
#'        data = liverlipidome,
#'        names_from = Lipid_Name,
#'        values_from = Lipid_Level
#'    )
#' # Create a mapping of the lipid names.
#' names.mapping <-
#'    map_lipid_names( x = unique( liverlipidome$"Lipid_Name" ) )
#' # Compute the regression models.
#' result.limma <-
#'    compute_models_with_limma(
#'        x = liverlipidome.wide,
#'        dependent.variables = names.mapping$"Name",
#'        independent.variables = c( "Diagnosis" ),
#'        F.test = TRUE # Compute an F-test for a factor variable.
#'    )
#' # Compute the F-test.
#' result.limma <- compute_F_test_with_limma( x = result.limma )
#' # Print a figure of the F-test.
#' \donttest{
#' figure.output <-
#'   heatmap_lipidome_from_limma(
#'       x = result.limma,
#'       names.mapping = names.mapping,
#'       F.test = TRUE,
#'       axis.x.carbons = FALSE,
#'       class.facet = "wrap",
#'       plot.all = FALSE,
#'       plot.individual = TRUE,
#'       scales = "free",
#'       space = "free"
#'   )
#' }
#' # Compute pairwise post-hoc comparisons between the factor levels for
#' # the dependent variables (i.e., lipids) with a significant F-test result.
#' result.limma <-
#'    compute_post_hoc_test_with_limma(
#'        x = result.limma,
#'        remap.level.names = TRUE
#'    )
#' # Print a figure of all post-hoc comparisons.
#' \donttest{
#' figure.output <-
#'     heatmap_lipidome_from_limma(
#'     x = result.limma$"result.post.hoc.test",
#'     names.mapping = names.mapping,
#'     axis.x.carbons = FALSE,
#'     plot.all = TRUE,
#'     plot.individual = FALSE,
#'     scales = "free",
#'     space = "free"
#' )
#' }
#' # Specify the contrasts of the post-hoc comparison that will be included
#' # in the figure.
#' contrasts.included <-
#'    c( "DiagnosisSteatosis", "DiagnosisNASH", "DiagnosisCirrhosis" )
#' # Get the omitted contrasts based on the above definition.
#' contrasts.omitted <-
#'    colnames( result.limma$"result.post.hoc.test"$"p.value" )[
#'        !(
#'            colnames( result.limma$"result.post.hoc.test"$"p.value" ) %in%
#'            contrasts.included
#'        )
#'    ]
#' # Find dependent variables (i.e., lipids) that have any significant
#' # difference.
#' has.any.significant <-
#'    apply(
#'        X =
#'            result.limma$"result.post.hoc.test"$"p.value"[
#'                ,
#'                contrasts.included
#'            ],
#'        MAR = 2,
#'        FUN = p.adjust,
#'        method = "BH"
#'    )
#' has.any.significant <-
#'    rownames(
#'        has.any.significant[
#'            apply(
#'                X = has.any.significant < 0.05,
#'                MAR = 1,
#'                FUN = any
#'            ),
#'        ]
#'    )
#' # Include in the figure only lipid classes that have at least four
#' # significant differences.
#' classes.included <-
#'    names(
#'        which(
#'            table(
#'                names.mapping[
#'                    make.names( has.any.significant ), "Class"
#'                ]
#'            ) > 4
#'        )
#'    )
#' classes.omitted <- unique( names.mapping$"Class" )
#' classes.omitted <-
#'    classes.omitted[ !( classes.omitted ) %in% classes.included ]
#' # Print a figure of the selected post-hoc-comparisons.
#' figure.output <-
#'    heatmap_lipidome_from_limma(
#'        x = result.limma$"result.post.hoc.test",
#'        names.mapping = names.mapping,
#'        axis.x.carbons = FALSE,
#'        omit.class = classes.omitted,
#'        omit.factor = contrasts.omitted,
#'        plot.all = TRUE,
#'        plot.individual = FALSE,
#'        scales = "free",
#'        space = "free"
#'    )
#'
"liverlipidome"







