compute_F_test_with_limma <-
  function( x,
            target.variable,
            p.adj.threshold = 0.05,
            print.table = TRUE ) {

    tmp <- options( "contrasts" )

    if ( any( x$"parameters"$"contrasts" != "contr.sum" ) ) {

      print(
        "ERROR: Call compute_models_with_limma() first with F.test = TRUE."
      )

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

          print(
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
