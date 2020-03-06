#' Compute an F-test for a multi-level factor independent variable
#'
#' Use this to achieve analysis of variance (ANOVA) or analysis of covariance
#'    (ANCOVA). The F-test is based on a model produced by the function
#'    \code{\link{compute_models_with_limma}}.
#'    To use this function, first call the function
#'    \code{\link{compute_models_with_limma}}.
#'
#' @param x (Required) list of models for which the test will be done.
#'    The F-test will be computed for the first independent variable that was
#'    specified in the \code{independent.variables} argument to the
#'    function \code{\link{compute_models_with_limma}}.
#' @param p.adj.threshold (Optional) numeric value specifying the threshold
#'    for statistical significance of difference after correction for multiple
#'    testing.
#' @param print.table (Optional) \code{TRUE} or \code{FALSE}: Print the results
#'    of the test?
#'
#' @return List \code{x} supplemented by the results of the F-test.
#'
#' @seealso \code{\link{compute_models_with_limma}} for the model computation
#'    step that is required prior to calling this function.
#' @seealso \code{\link{compute_post_hoc_test_with_limma}} for the pairwise
#'    post-hoc comparisons that may follow the F-test done with this function.
#'
#' @export
#'
compute_F_test_with_limma <-
  function( x,
            p.adj.threshold = 0.05,
            print.table = FALSE ) {

    if ( any( x$"parameters"$"contrasts" != "contr.sum" ) ) {

      stop( "Call compute_models_with_limma() first with F.test = TRUE." )

      return( NULL )

    } else {

      idx.target <-
        grep(
          x = colnames( x$"model"$"coefficients" ),
          pattern = paste0("^", x$"parameters"$"independent.variables"[ 1 ] )
        )

      names.target <- colnames( x$"model"$"coefficients" )[ idx.target ]

      y <-
        limma::topTable(
          fit = x$"model",
          coef = names.target,
          confint = TRUE,
          number = Inf,
          p.value = 1,
          adjust.method = "BH"
        )

      if ( print.table ) {

        tmp <- which( y$"adj.P.Val" < p.adj.threshold )

        if ( length( tmp ) > 0 ) {

          y.printed <- y[ tmp, , drop = FALSE ]

          y.printed <-
            signif(
              x = y.printed,
              digits = 3
            )

          y.printed$"Name" <- rownames( y.printed )

          y.printed$"Name" <-
            stringr::str_sub(
              string = y.printed$"Name",
              start = 1,
              end = 25
            )

          y.printed <-
            y.printed[ , c( "Name", "AveExpr", "F", "P.Value", "adj.P.Val" ) ]

          y.printed <-
            knitr::kable(
              x = y.printed,
              caption =
                paste(
                  "F-test for",
                  paste(
                    names.target,
                    collapse = ", "
                  )
                ),
              row.names = FALSE )

          print( y.printed )

        } else {

          message(
            paste(
              "No significant F-tests at p.adj <",
              p.adj.threshold
            )
          )

        }

      }

      y.out <- x
      y.out$"result.F.test" <- y

      y.out$"parameters"$"p.adj.threshold" <- p.adj.threshold

      return( y.out )

    }

  }
