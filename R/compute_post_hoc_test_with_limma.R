#' Compute pairwise post-hoc comparisons for a multi-level factor
#'
#' Use this function to achieve the post-hoc comparisons between the
#'    multiple levels of an independent variable. These comparisons follow
#'    analysis of variance (ANOVA) or analysis of covariance (ANCOVA).
#'    The pairwise comparisons are based on a result of an F-test produced by
#'    the function \code{\link{compute_F_test_with_limma}}.
#'    To use this function, first call the functions
#'    \code{\link{compute_models_with_limma}} and
#'    \code{\link{compute_F_test_with_limma}} consecutively.
#'
#' @param x (Required) list of models for which the test will be done.
#'    The pairwise comparisons will be computed using the levels of the first
#'    independent variable that was specified in
#'    the \code{independent.variables} argument to the
#'    function \code{\link{compute_models_with_limma}}. The comparisons will be
#'    computed for the lipids that had an F-test result of statistical
#'    significance from the function \code{\link{compute_F_test_with_limma}}.
#' @param p.adj.threshold (Optional) numeric value specifying the threshold of
#'    statistical significance in the pairwise comparisons after a correction
#'    for multiple testing. We recommend to leave this argument unfilled,
#'    leading to the same threshold to be used as in the preceding F-test.
#' @param remap.level.names (Optional) \code{TRUE} or \code{FALSE}: Should
#'    the levels of the factor independent variable be re-mapped? This feature
#'    can be used to solve a problem with the factor levels, which may occur
#'    in some versions of the \code{limma} package. We recommend to keep this
#'    argument unchanged from the default value.
#'
#' @return List \code{x} supplemented by the results of the pairwise
#'    post-hoc comparisons.
#'
#' @seealso \code{\link{compute_models_with_limma}} for the model computation
#'    step that is required prior to calling this function.
#' @seealso \code{\link{compute_F_test_with_limma}} for the F-test step that is
#'    required prior to calling this function.
#'
#' @export
#'
compute_post_hoc_test_with_limma <-
  function(
    x,
    p.adj.threshold = NULL,
    remap.level.names = FALSE
  ) {

    options.original <- options()
    on.exit( expr = options( options.original ) )

    # Set appropriate contrasts.

    options( contrasts = c( "contr.treatment", "contr.poly" ) )

    if ( is.null( p.adj.threshold ) ) {

      p.adj.threshold <- x$"parameters"$"p.adj.threshold"

    }

    # Re-fit linear regression models

    names.post.hoc <-
      rownames(
        x$"result.F.test"
      )[ which( x$"result.F.test"$"adj.P.Val" < p.adj.threshold ) ]

    if ( length( names.post.hoc ) == 0 ) {

      message(
        paste(
          "No Significant F-test results at p.adj <",
          p.adj.threshold
        )
      )

      return( NULL )

    } else {

      # Re-compute design matrix with contrasts for post-hoc tests.

      x$"data"$"design.post.hoc" <-
        stats::model.matrix(
          object = x$"parameters"$"formula",
          data = x$"data"$"predictors"
        )

      if ( is.null( x$"random.effect" ) ) {

        rFit <-
          limma::lmFit(
            object = x$"data"$"object"[ names.post.hoc, , drop = FALSE ],
            design = x$"data"$"design.post.hoc"
          )

      } else {

        rFit <-
          limma::lmFit(
            object = x$"data"$"object"[ names.post.hoc, , drop = FALSE ],
            design = x$"data"$"design.post.hoc",
            block = y$"random.effect"$"block",
            correlation = cFit$"consensus.correlation"
          )

      }

      rEbFit <- limma::eBayes( fit = rFit )

      # Generate contrasts.

      names.target <-
        colnames( x$"result.F.test" )

      names.target <-
        names.target[
          !(
            names.target %in%
              c(
                "AveExpr",
                "F",
                "P.Value",
                "adj.P.Val"
              )
          )
          ]

      contrasts.posthoc <- names.target

      for ( i in 2:length( names.target ) ) {

        for ( j in 1:( i - 1 ) ) {

          contrasts.posthoc <-
            c(
              contrasts.posthoc,
              paste(
                names.target[ i ],
                "-",
                names.target[ j ]
              )
            )

        }

      }

      contrasts.posthoc <-
        limma::makeContrasts(
          contrasts = contrasts.posthoc,
          levels = x$"data"$"design"
        )

      cFit <-
        limma::contrasts.fit(
          fit = rEbFit,
          contrasts = contrasts.posthoc
        )

      cEbFit <- limma::eBayes( fit = cFit )

      y <- x

      y$"result.post.hoc.test" <- cEbFit

      if ( remap.level.names ) {

        tmp <-
          levels(
            x =
              x$"data"$"predictors"[
                ,
                x$"parameters"$"independent.variables"[ 1 ]
                ]
          )

        mapping.levels <-
          data.frame(
            label = tmp,
            number = 0:( length( tmp ) - 1 )
          )

        mapping.levels$"combination.number" <-
          paste0(
            x$"parameters"$"independent.variables"[ 1 ],
            mapping.levels$"number"
          )

        mapping.levels$"combination.label" <-
          paste0(
            x$"parameters"$"independent.variables"[ 1 ],
            mapping.levels$"label"
          )

        for ( i in 1:nrow( mapping.levels ) ) {

          colnames( y$"result.post.hoc.test"$"coefficients" ) <-
            colnames( y$"result.post.hoc.test"$"p.value" ) <-
            stringr::str_replace(
              string = colnames( y$"result.post.hoc.test"$"coefficients" ),
              pattern = mapping.levels[ i, "combination.number" ],
              replacement = mapping.levels[ i, "combination.label" ]
            )

        }

      }

      y$"result.post.hoc.test"$"formula" <- x$"parameters"$"formula"
      y$"result.post.hoc.test"$"formula.text" <- x$"parameters"$"formula.text"

      return( y )

    }

  }
