print_tables_from_limma <-
  function(
    x,
    all = FALSE,
    coefs = NULL,
    digits = 3,
    formula.max.nchar = 60,
    p.adj.method = "BH",
    p.adj.threshold = 0.05,
    print = TRUE,
    name.max.nchar = 30,
    omit.first.coef = TRUE
  ) {

    if ( is.null( coefs ) ) {

      if ( omit.first.coef ) {

        coefs <- colnames( x$"coefficients" )[ -1 ]

      } else {

        coefs <- colnames( x$"coefficients" )

      }

    }

    if ( all ) {

      tables.coefs <-
        vector(
          mode = "list",
          length = length( coefs ) + 1
        )

      names( tables.coefs ) <- c( coefs, "All" )

      table.all <- NULL

    } else {

      tables.coefs <-
        vector(
          mode = "list",
          length = length( coefs )
        )

      names( tables.coefs ) <- coefs

    }

    tables.coefs.printed <- tables.coefs

    printout.model <-
      paste(
        as.character( x$"formula" ),
        collapse = " "
      )

    printout.model <-
      stringr::str_wrap(
        string = printout.model,
        width = formula.max.nchar,
        indent = 3,
        exdent = 5
      )

    printout.model <-
      stringr::str_split(
        string = printout.model,
        pattern = "\\n"
      )[[ 1 ]]

    printout.model[ length( printout.model ) ] <-
      paste(
        printout.model[ length( printout.model ) ],
        ")",
        sep = ""
      )

    for ( i in 1:length( coefs ) ) {

      coef.i <- coefs[ i ]

      tables.coefs[[ coef.i ]] <-
        tables.coefs.printed[[ coef.i ]] <-
        limma::topTable(
          fit = x,
          coef = coef.i,
          number = Inf,
          p.value = p.adj.threshold
        )

      if ( all & nrow( tables.coefs[[ coef.i ]] ) > 0 ) {

        if ( is.null( table.all ) ) {

          table.all <- tables.coefs[[ coef.i ]]
          table.all$"Coefficient" <- coef.i
          table.all$"Name" <- rownames( table.all )

        } else {

          tmp <- tables.coefs[[ coef.i ]]
          tmp$"Coefficient" <- coef.i
          tmp$"Name" <- rownames( tmp )

          table.all <- rbind( table.all, tmp )

        }

      }

      if ( print ) {

        print( "" )
        print( paste( "Table:", coef.i ) )
        print( "  (from model: " )

        for ( i in 1:length( printout.model ) ) {

          print( printout.model[ i ] )

        }

        print( "" )

        if ( nrow( tables.coefs[[ coef.i ]] ) > 0 ) {

          if ( nrow( tables.coefs[[ coef.i ]] ) > 1 ) {

            tables.coefs.printed[[ coef.i ]] <-
              apply(
                X = tables.coefs.printed[[ coef.i ]],
                MARGIN = 2,
                FUN = signif,
                digits = 3
              )

          } else if ( nrow( tables.coefs[[ coef.i ]] ) == 1 ) {

            tables.coefs.printed[[ coef.i ]] <-
              signif(
                x = tables.coefs.printed[[ coef.i ]],
                digits = 3
              )

          }

          tables.coefs.printed[[ coef.i ]] <-
            data.frame(
              Name = rownames( tables.coefs.printed[[ coef.i ]] ),
              tables.coefs.printed[[ coef.i ]][
                ,
                c(
                  "logFC",
                  "P.Value",
                  "adj.P.Val"
                )
                ],
              check.names = FALSE,
              stringsAsFactors = FALSE
            )

          colnames(
            tables.coefs.printed[[ coef.i ]]
          )[ colnames( tables.coefs.printed[[ coef.i ]] ) == "logFC" ] <-
            "Coefficient"

          tables.coefs.printed[[ coef.i ]]$"Name" <-
            stringr::str_sub(
              string = tables.coefs.printed[[ coef.i ]]$"Name",
              start = 1,
              end = name.max.nchar
            )

          rownames( tables.coefs.printed[[ coef.i ]] ) <- NULL

          print( tables.coefs.printed[[ coef.i ]] )

        } else {

          print(
            paste(
              "No significant associations at p.adj <",
              p.adj.threshold
            )
          )

        }

      }

    }

    if ( all ) {

      table.all$"adj.P.Val" <-
        p.adjust( p = table.all$"P.Value", method = "BH" )

      tables.coefs$"All" <- table.all

    }

    return( tables.coefs )

  }
