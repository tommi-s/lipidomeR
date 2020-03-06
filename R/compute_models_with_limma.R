#' Compute lipid-specific regression models
#'
#' Use this function to computing multiple regression models that can be
#'    directly supplied to the visualization functions of the 'lipidomeR'.
#'
#' @param x (Required) data matrix.
#' @param dependent.variables
#'    (Required) vector of names of dependent variables.
#'    These should be the names of the lipids.
#' @param independent.variables
#'    (Required) vector of names of the independent variables.
#'    These should be the names of the variables defining the experiment design.
#' @param random.effect (Optional) name of a single variable specifying
#'    the random effect for a random-effects model.
#'    For instance, \code{ID} specifies a random effect as in
#'    \code{limma::duplicateCorrelation( ..., block = x$'ID' )}.
#' @param formula
#'    (Optional) character string of model formula in the format accepted by
#'    the function \code{\link[stats]{model.matrix}} and starting with \code{~}.
#'    Variables mentioned in the formula should be included in
#'    the \code{independent.variables} argument.
#'    For instance, \code{Group * Treatment}.
#' @param F.test (Optional) \code{TRUE} or \code{FALSE}:
#'    Should an F-test for analysis of variance (ANOVA) or
#'    analysis of covariance (ANCOVA) be computed?
#' @param print.table1 (Optional) \code{TRUE} or \code{FALSE}:
#'    Should a summary table of the independent variables be printed?
#' @param scale.dependent.variables (Optional) \code{TRUE} or \code{FALSE}:
#'    Should dependent variables be scaled to zero-mean and unit-variance
#'    prior to model fitting?
#' @param scale.independent.variables (Optional) \code{TRUE} or \code{FALSE}:
#'    Should independent variables be scaled to zero-mean and unit-variance
#'    prior to model fitting?
#' @param verbose (Optional) \code{TRUE} or \code{FALSE}: Print messages from
#'    the model fitting?
#'
#' @return List of regression results the that can be directly supplied as
#'    an argument to the function \code{\link{heatmap_lipidome_from_limma}} and
#'    other visualization functions of the lipidomeR.
#'
#' @seealso \code{\link{heatmap_lipidome_from_limma}} for visualizing the
#'    output of this function.
#'
#' @export
#'
compute_models_with_limma <-
  function(
    x,
    dependent.variables,
    independent.variables,
    random.effect = NULL,
    formula = NULL,
    F.test = FALSE,
    print.table1 = FALSE,
    scale.dependent.variables = TRUE,
    scale.independent.variables = FALSE,
    verbose = TRUE
  ) {

    options.original <- options()
    on.exit( expr = options( options.original ) )

    # Extract independent variables.

    design.test <-
      data.frame(
        x[ , c( independent.variables, random.effect ), drop = FALSE ]
      )

    # Subset samples with complete data of independent variables.

    tmp <-
      apply(
        X = !is.na( design.test ),
        MARGIN = 1,
        FUN = all
      )

    if ( !any ( tmp ) ) {

      stop( "No samples with non-missing independent variables." )

    }

    design.test <- design.test[ tmp, , drop = FALSE ]

    design.test <- droplevels( design.test )

    # Set appropriate contrasts.

    if ( F.test ) {

      options( contrasts = c( 'contr.sum', 'contr.sum' ) )

      design.test[ , independent.variables[ 1 ] ] <-
        factor( x = design.test[ , independent.variables[ 1 ] ] )

    } else {

      options( contrasts = c( 'contr.treatment', 'contr.poly' ) )

    }

    y <- list()

    y$"parameters"$"contrasts" <- options( "contrasts" )$"contrasts"

    if ( print.table1 ) {

      if ( length( unique( independent.variables[ 1 ] ) < 10 ) &
           length( independent.variables > 1 ) ) {

        table1 <-
          tableone::CreateTableOne( vars = independent.variables[ -1 ],
                                    data = design.test,
                                    strata = independent.variables[ 1 ] )

      } else {

        table1 <-
          tableone::CreateTableOne( vars = independent.variables,
                                    data = design.test )

      }

      print( "Table of independent variables:" )
      print( table1 )
      print( "" )

    }

    if ( scale.independent.variables ) {

      tmp <-
        sapply( X = design.test[ 1, ],
                FUN = is.numeric )

      if ( any( tmp ) ) {

        design.test[ , tmp ] <- scale( x = design.test[ , tmp ] )

      }

    }

    # Extract dependent variables

    data.test <- x[ tmp, dependent.variables ]

    if ( scale.dependent.variables ) {

      data.test <- scale( x = data.test )

    }

    data.test <- t( data.test )

    # Create the formula as an additive model of independent variables.

    if ( is.null( formula ) ) {

      formula <-
        as.formula(
          paste(
            "~",
            paste(
              independent.variables,
              collapse = " + "
            )
          )
        )

    } else {

      formula <- as.formula( formula )

    }

    if ( verbose ) {

      tmp <- paste( "Fitting models:" )
      tmp <-
        paste(
          tmp,
          paste( as.character( formula ), collapse = " " )
        )

      if ( !is.null( random.effect ) ) {

        tmp <-
          paste(
            tmp,
            paste( "with random effect:", random.effect )
          )

      }

      message( tmp )

    }

    # Create model matrix.

    model.matrix <-
      stats::model.matrix(
        object = formula,
        data = design.test
      )

    y$"data" <- list()

    y$"data"$"design" <- model.matrix

    y$"data"$"object" <- data.test

    y$"data"$"predictors" <- design.test

    if ( is.null( random.effect ) ) {

      # Fit the linear regression models.

      mFit <-
        limma::lmFit(
          object = data.test,
          design = model.matrix
        )

    } else {

      # Fit the linear mixed effects model.

      y$"random.effect"$"block" <-
        factor( x = design.test[ , random.effect ] )

      cFit <-
        limma::duplicateCorrelation(
          object = data.test,
          design = model.matrix,
          block = y$"random.effect"$"block"
        )

      mFit <-
        limma::lmFit(
          object = data.test,
          design = model.matrix,
          block = y$"random.effect"$"block",
          correlation = cFit$"consensus.correlation"
        )

      y$"random.effect"$"correlation" <- cFit$"consensus.correlation"

    }

    # Fit the Bayesian noise model.

    y$"model" <- limma::eBayes( fit = mFit )

    y$"parameters"$"formula" <-
      y$"model"$"formula" <-
      formula

    y$"parameters"$"formula.text" <-
      y$"model"$"formula.text" <-
      paste(
        as.character( y$"parameters"$"formula" ),
        collapse = " "
      )

    if ( !is.null( random.effect ) ) {

      y$"parameters"$"formula.text" <-
        y$"model"$"formula.text" <-
        paste(
          y$"parameters"$"formula.text",
          " + (1|",
          random.effect,
          ")",
          sep = ""
        )

    }

    y$"parameters"$"independent.variables" <-
      independent.variables

    y$"parameters"$"levels.target" <-
      levels( x = design.test[ , independent.variables[ 1 ] ] )

    # Return the result.

    return( y )

  }
