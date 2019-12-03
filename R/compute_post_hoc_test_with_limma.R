compute_post_hoc_test_with_limma <-
  function(
    x,
    p.adj.threshold = NULL,
    remap.level.names = FALSE
  ) {

    # Set appropriate contrasts.

    options( contrasts = c( "contr.treatment", "contr.poly" ) )

    if ( is.null( p.adj.threshold ) ) {

      p.adj.threshold <- x$"parameters"$"p.adj.threshold"

    }

    # Re-fit linear regression models

    names.post.hoc <-
      rownames(
        x$"result.F.test"
      )[ x$"result.F.test"$"adj.P.Val" < p.adj.threshold ]

    if ( length( names.post.hoc ) == 0 ) {

      print(
        paste(
          "WARNING: No Significant F-test results at p.adj <",
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
