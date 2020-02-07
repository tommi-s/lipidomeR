#' Create 'lipidomeR' heatmaps of model statistics
#'
#' Use this function to creating heatmaps of model statistics from
#'    the output of the \code{\link{compute_models_with_limma}} function or
#'    related functions (\code{\link{compute_F_test_with_limma}} and
#'    \code{\link{compute_post_hoc_test_with_limma}}).
#' @param x (Required) list of output from
#'    the \code{\link{compute_models_with_limma}} function or related functions.
#' @param names.mapping (Optional) mapping of lipid names from
#'    the \code{\link{map_lipid_names}} function.
#' @param axis.x.carbons (Optional) \code{TRUE} or \code{FALSE}: Should
#'    the lipid size (i.e., number of carbon atoms in the fatty acid chain) be
#'    presented on the x-axis or y-axis?
#' @param class.facet (Optional) character string with possible values
#'    \code{'col'}, \code{'row'} or \code{'wrap'}:
#'    Present lipid classes as panels organized into columns, rows or into a
#'    wrapped layout spanning over multiple rows and columns. The alternative
#'    \code{'wrap'} is only available with \code{plot.infividual = TRUE}.
#' @param class.subset (Optional) character vector specifying a subset of
#'    the lipid classes (e.g., \code{c( 'PC', 'SM', 'TG' )}). Same outcome
#'    can be achieved with the argument \code{omit.class}.
#' @param F.test (Optional) \code{TRUE} or \code{FALSE}: Should the result of
#'    an F-test be visualized instead of individual model coefficients?
#'    Result from the function \code{\link{compute_F_test_with_limma}} needs
#'    to be provided as argument \code{x} for the
#'    \code{\link{heatmap_lipidome_from_limma}} function along with
#'    the argument \code{F.test = TRUE}.
#' @param omit.class (Optional) character vector of lipid classes omitted from
#'    the visualization (e.g., \code{c( 'PC', 'SM', 'TG' )}).
#' @param omit.factor (Optional) character vector of lipid classes omitted from
#'    the visualization (e.g., \code{c( 'Time' )}).
#' @param order.factor (Optional) \code{TRUE} or \code{FALSE}: Sort the factor
#'    panels alphabetically?
#' @param p.val.thresholds (Optional) numeric vector with increasing values of
#'    highlighting thresholds for multiple-testing-corrected p-values.
#'    For instance, \code{c( 0.01, 0.05, 0.1 )} gives distinct highlighting of
#'    associations in the ranges \eqn{p < 0.01}, \eqn{0.01 \geq p < 0.05}, and
#'    \eqn{0.05 \geq p < 0.1}. If the number of these categories is changed from
#'    the default value, also the argument \code{p.val.labels} needs to be
#'    modified to match it.
#' @param p.val.labels (Optional) numeric vector of point characters for the
#'    p-value highlighting categories specified in the argument
#'    \code{p.val.thresholds}. The values need to be compatible with
#'    the function \code{\link[ggplot2]{scale_shape_manual}}.
#'    If the number of these labels is changed from the default value, also
#'    the argument \code{p.val.thresholds} needs to be modified to match it.
#' @param p.val.label.bg.size (Optional) numeric value to scale the size of
#'    the colored background of the symbols that indicate values with
#'    statistical significance. The background is used to ensure that
#'    the white symbols are visible also when the color of a heatmap
#'    rectangle is of a bleak color.
#' @param p.val.label.size (Optional) numeric value to scale the size of the
#'    symbols that indicate values with statistical significance.
#' @param p.adj.method (Optional) name of the method to correct p-values for
#'    multiple testing. Accepted values are as in the function
#'    \code{stats::p.adjust()}.
#' @param plot.individual (Optional) \code{TRUE} or \code{FALSE}: Create
#'    specific figures for each coefficient of the model? For models with many
#'    independent variables, this may take some time to complete.
#' @param plot.all (Optional) \code{TRUE} or \code{FALSE}: Create
#'    a combined figure of all the coefficients of the model? For models with
#'    many independent variables, this may lead to too busy a figure.
#' @param print.figure (Optional) \code{TRUE} or \code{FALSE}: Print created
#'    figure(s) to the standard output? Note that when
#'    \code{plot.individual = TRUE}, multiple figures will be printed
#'    consecutively.
#' @param print.formula (Optional) \code{TRUE} or \code{FALSE}: Add
#'    the model formula as a title to the figure? We recommend to use
#'    \code{print.formula = TRUE} to make the result more understandable and
#'    reproducible.
#' @param formula.width (Optional) numeric value to specify the width of a line
#'    in the model formula. Relevant only if \code{print.formula = TRUE}.
#' @param legend.key.size.multiplier (Optional) numeric value to scale the size
#'    of the figure (key) legends.
#' @param range.min.N.carbons (Optional) numeric value to specify the minimum
#'    range of the axis showing the lipid size (number of carbon atoms in the
#'    fatty acid chains). This value can be increased from the default value to
#'    improve readability in situations, where there are lipid classes with
#'    little or no variation in the lipid size.
#' @param range.min.N.double.bonds (Optional) numeric value to specify
#'    the minimum range of the axis showing the lipid saturation (number of
#'    double bonds in the fatty acid chains). This value can be increased from
#'    the default value to improve readability in situtions, where there are
#'    lipid classes with little or no variation in the lipid saturation.
#' @param scales (Optional) character string with possible values
#'    \code{'fixed'}, \code{'free'}, \code{'free_x'} or \code{'free_y'}. This
#'    argument specifies, whether the axes in multiple sub-heatmaps will be in
#'    the same scale (\code{'fixed'}) or in a scale specific to each sub-figure.
#'    See the function \code{\link[ggplot2]{facet_grid}} for details.
#' @param shadowtext (Optional) \code{TRUE} or \code{FALSE}: Should an
#'    alternative method of statistical significance highlighting be used in
#'    the figure? If \code{TRUE}, highlighting will be based on the function
#'    \code{\link[shadowtext]{geom_shadowtext}}.
#' @param space (Optional) character string with possible values
#'    \code{'fixed'}, \code{'free'}, \code{'free_x'} or \code{'free_y'}.
#'    This argument specifies, whether the sub-heatmaps will be of identical
#'    size (\code{'fixed'}) or not.
#'    See the function \code{\link[ggplot2]{facet_grid}} for details.
#' @param verbose (Optional) \code{TRUE} or \code{FALSE}: Print information
#'    about the progress of the function call?
#' @param wrap.contrast.name (Optional) \code{TRUE} or \code{FALSE}: Wrap
#'    the name of a contrast to multiple lines to ensure readability?
#' @param baseline.adjusted \code{TRUE} or \code{FALSE}: Is the model object
#'    that is specified as argument \code{x} a baseline-adjusted regression
#'    model? (Note: this model type is not yet provided in the package.)
#' @param survival (Optional) \code{TRUE} or \code{FALSE}: Is the model object
#'    that is specified as argument \code{x} a survival model?
#'    (Note: this model type is not yet provided in the package.)
#' @return List of 'lipidomeR' heatmap figure(s).
#' @seealso \code{\link{compute_models_with_limma}} for computing the argument
#'    \code{x} for this function.
#'
#' @export
#'
heatmap_lipidome_from_limma <-
  function(
    x,
    names.mapping = NULL,
    axis.x.carbons = TRUE,
    baseline.adjusted = FALSE,
    class.facet = "row",
    class.subset = NULL,
    F.test = FALSE,
    omit.class = NULL,
    omit.factor = NULL,
    order.factor = FALSE,
    p.val.thresholds = c( 0.01, 0.05, 0.1 ),
    p.val.labels = c( 8, 4, 3 ),
    p.val.label.bg.size = 2,
    p.val.label.size = 1,
    p.adj.method = "BH",
    plot.individual = FALSE,
    plot.all = TRUE,
    print.figure = TRUE,
    print.formula = TRUE,
    formula.width = 110,
    legend.key.size.multiplier = 2,
    range.min.N.carbons = 5,
    range.min.N.double.bonds = 5,
    scales = "fixed",
    shadowtext = FALSE,
    space = "free",
    survival = FALSE,
    verbose = FALSE,
    wrap.contrast.name = TRUE
  ) {

    if ( verbose ) {

      tmp <-
        paste(
          "heatmap_lipidome_from_limma() was created by Tommi Suvitaival,",
          "tommi.raimo.leo.suvitaival@regionh.dk"
        )

      message( tmp )

    }

    Class <- Factor <- NULL # Avoid note with 'ggplot2::Vars()'.

    x.in <- x

    p.val.cat <-
      data.frame(
        Label = p.val.labels,
        Relation = rep( x = "p.adj <", times = 3 ),
        p.adj = p.val.thresholds
      )

    if ( is.null( names.mapping ) ) {

      if ( F.test ) {

        names.mapping <- map_lipid_names( x = rownames( x$"result.F.test" ) )

      } else if ( survival ) {

        names.mapping <- map_lipid_names( x = rownames( x$"result.survival" ) )

      } else {

        names.mapping <- map_lipid_names( x = rownames( x$"coefs" ) )

      }

    }

    if ( !is.null( class.subset ) ) {

      names.mapping <-
        names.mapping[ which( names.mapping$"Class" %in% class.subset ), ]

    }

    ### Create heatmap of F-test
    ### or heatmaps of individual coefficients.

    if ( F.test ) { # Result of F-test.

      if ( is.null( x$"result.F.test" ) ) {

        stop( "Call compute_F_test_with_limma() first." )

        return( NULL )

      }

      y <- x$"result.F.test"

      title <-
        colnames( y )[
          !( colnames( y ) %in%
               c(
                 "AveExpr",
                 "F",
                 "P.Value",
                 "adj.P.Val"
               )
          )
          ]

      title <-
        paste(
          title,
          collapse = ", "
        )

      title <- paste( "F:", title )

      y$"Name" <- rownames( y )

      y$"Name.checked" <- make.names( y$"Name" )

      y$"P_Value" <-
        cut(
          x = y$"adj.P.Val",
          breaks = c( 0 - .Machine$"double.eps", p.val.thresholds ),
          labels = paste( "p <", p.val.thresholds ),
          include.lowest = TRUE
        )

      y$"p.adj.range" <-
        cut(
          x = y$"adj.P.Val",
          breaks = c( 0, p.val.thresholds ),
          labels = p.val.labels,
          include.lowest = TRUE
        )

      y$"Significance" <- rep( x = NA, times = nrow( y ) )

      y[ !is.na( y$"P_Value" ), "Significance" ] <- "Significant"

      y <-
        merge(
          x = y,
          y = names.mapping,
          by = "Name.checked",
          all.x = FALSE, # TRUE,
          all.y = FALSE
        )

      if ( print.formula ) {

        printout.formula <-
          paste(
            as.character( x$"parameters"$"formula.text" ),
            collapse = " "
          )

        printout.formula <-
          paste(
            "Model:",
            printout.formula
          )

        printout.formula <-
          split_formula(
            x = printout.formula,
            width = formula.width
          )

      } else {

        printout.formula <- NULL

      }

      if ( axis.x.carbons ) {

        plot.i <-
          ggplot2::ggplot(
            data = y,
            mapping =
              ggplot2::aes_string(
                x = "N.carbons",
                y = "N.double.bonds",
                fill = "F"
              )
          ) +
          ggplot2::ylab( label = "Number of fatty-acid double bonds" ) +
          ggplot2::xlab( label = "Number of fatty-acid carbon atoms" )

      } else {

        plot.i <-
          ggplot2::ggplot(
            data = y,
            mapping =
              ggplot2::aes_string(
                x = "N.double.bonds",
                y = "N.carbons",
                fill = "F"
              )
          ) +
          ggplot2::xlab( label = "Number of fatty-acid double bonds" ) +
          ggplot2::ylab( label = "Number of fatty-acid carbon atoms" )

      }

      plot.i <- plot.i +
        ggplot2::theme_dark() +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
          low = "white",
          high = "red"
        ) +
        ggplot2::ggtitle(
          label = title,
          subtitle = printout.formula
        ) +
        ggplot2::theme(
          legend.position = "top",
          legend.text = ggplot2::element_text( angle = 45 )
        ) +
        ggplot2::theme(
          axis.text.x =
            ggplot2::element_text(
              angle = 45,
              vjust = 0.5
            )
        ) +
        ggplot2::theme( strip.text.y = ggplot2::element_text( angle = 0 ) ) +
        ggplot2::scale_x_continuous(
          breaks = get_integer_breaks,
          expand = ggplot2::expand_scale( mult = 0, add = 0.5 )
        ) +
        ggplot2::scale_y_continuous(
          breaks = get_integer_breaks,
          expand = ggplot2::expand_scale( mult = 0, add = 0.5 )
        )

      if ( shadowtext ) {

        plot.i <-
          plot.i +
          shadowtext::geom_shadowtext(
            mapping = ggplot2::aes_string( label = "p.adj.range" ),
            color = "red",
            bg.colour = "white"
          )

      } else {

        plot.i <-
          plot.i +
          ggplot2::geom_point(
            mapping =
              ggplot2::aes_string(
                color = "Significance"
              ),
            na.rm = TRUE,
            shape = 19,
            size = p.val.label.bg.size
          ) +
          ggplot2::scale_color_manual(
            values = c( "red" ),
            drop = FALSE,
            na.translate = FALSE
          ) +
          ggplot2::geom_point(
            mapping =
              ggplot2::aes_string(
                shape = "P_Value"
              ),
            color = "white",
            size = p.val.label.size
          ) +
          ggplot2::scale_shape_manual(
            values = p.val.labels,
            drop = FALSE,
            na.translate = FALSE
          ) +
          ggplot2::guides(
            fill = ggplot2::guide_colorbar( order = 1 ),
            color =
              ggplot2::guide_legend(
                order = 2,
                override.aes =
                  list(
                    size = legend.key.size.multiplier * p.val.label.bg.size
                  )
              ),
            shape =
              ggplot2::guide_legend(
                order = 3,
                override.aes =
                  list( size = legend.key.size.multiplier * p.val.label.size )
              )
          )

      }

      if ( class.facet == "row" ) {

        plot.i <-
          plot.i +
          ggplot2::facet_grid(
            rows = ggplot2::vars( Class ), # Class ~ .,
            scales = scales
            ,
            space = space
          )

      } else if ( class.facet == "col" ) {

        plot.i <-
          plot.i +
          ggplot2::facet_grid(
            cols = ggplot2::vars( Class ), # ~ Class,
            scales = scales
            ,
            space = space
          )

      } else {

        plot.i <-
          plot.i +
          ggplot2::facet_wrap(
            facets = ggplot2::vars( Class ),
            scales = scales
          )

      }

      if ( class.facet == "wrap" ) {

        if (
          !is.null( range.min.N.carbons ) |
          !is.null( range.min.N.double.bonds )
        ) {

          data.blank <-
            create_blank_data_for_range(
              x = y,
              range.min.N.carbons = range.min.N.carbons,
              range.min.N.double.bonds = range.min.N.double.bonds
            )

          plot.i <-
            plot.i +
            ggplot2::geom_blank( data = data.blank )

        }

      }

      figures <- list()

      figures$"F" <- plot.i

      if ( print.figure ) {

        print( plot.i )

      }


    } else { # Results of individual coefficients.

      if ( survival | baseline.adjusted ) {

        if ( baseline.adjusted ) {

          y <-
            dplyr::bind_rows(
              x$"result.baseline.adjusted",
              .id = "Factor"
            )

          colnames( y )[ colnames( y ) == "Estimate" ] <-
            "Coefficient"

        } else if ( survival ) {

          y <- x$"result.survival"

          y$"Factor" <-
            rep(
              x = x$"parameters"$"event",
              times = nrow( y )
            )

          y$"Coefficient" <- y$"exp(coef)"

        }

        y$"P_Value" <-
          cut( x = y$"p.adj",
               breaks = c( 0 - .Machine$"double.eps", p.val.thresholds ),
               labels = paste( "p <", p.val.thresholds ),
               include.lowest = TRUE )

        y$"p.adj.range" <-
          cut( x = y$"p.adj",
               breaks = c( 0, p.val.thresholds ),
               labels = p.val.labels,
               include.lowest=TRUE )

        y$"p.adj.range.pos" <-
          y$"p.adj.range.neg" <-
          y$"p.adj.range"

        y$"Significant_Coefficient_Sign" <-
          rep( x = NA, times = nrow( y ) )

        if ( survival ) {

          y[ which( y$"Coefficient" <= 1 ), "p.adj.range.pos" ] <- NA
          y[ which( y$"Coefficient" > 1 ), "p.adj.range.neg" ] <- NA

          y[ !is.na( y$"p.adj.range.pos" ), "Significant_Coefficient_Sign" ] <-
            "Increased\nRisk"

          y[ !is.na( y$"p.adj.range.neg" ), "Significant_Coefficient_Sign" ] <-
            "Reduced\nRisk"

          y$"Significant_Coefficient_Sign" <-
            factor(
              x = y$"Significant_Coefficient_Sign",
              levels = c( "Reduced\nRisk", "Increased\nRisk" )
            )

        } else {

          y[ which( y$"Coefficient" <= 0 ), "p.adj.range.pos" ] <- NA
          y[ which( y$"Coefficient" > 0 ), "p.adj.range.neg" ] <- NA

          y[ !is.na( y$"p.adj.range.pos" ), "Significant_Coefficient_Sign" ] <-
            "Positive"

          y[ !is.na( y$"p.adj.range.neg" ), "Significant_Coefficient_Sign" ] <-
            "Negative"

          y$"Significant_Coefficient_Sign" <-
            factor(
              x = y$"Significant_Coefficient_Sign",
              levels = c( "Negative", "Positive" )
            )

        }

      } else {

        ### Gathering

        coefs <- x$"coefficients"

        coefs <-
          data.frame( Name = rownames( coefs ),
                      coefs,
                      stringsAsFactors = FALSE,
                      check.names = FALSE )

        coefs <-
          tidyr::gather( data = coefs,
                         key = "Factor",
                         value = "Coefficient",
                         -1 )


        pvals <-
          apply( X = x$"p.value",
                 MARGIN = 2,
                 FUN = p.adjust,
                 method = p.adj.method )

        pvals <-
          data.frame( Name = rownames( pvals ),
                      pvals,
                      stringsAsFactors = FALSE,
                      check.names = FALSE )

        pvals <-
          tidyr::gather( data = pvals,
                         key = "Factor",
                         value = "p.adj",
                         -1 )

        pvals$"P_Value" <-
          cut( x = pvals$"p.adj",
               breaks = c( 0 - .Machine$"double.eps", p.val.thresholds ),
               labels = paste( "p <", p.val.thresholds ),
               include.lowest = TRUE )

        pvals$"p.adj.range" <-
          cut( x = pvals$"p.adj",
               breaks = c( 0, p.val.thresholds ),
               labels = p.val.labels,
               include.lowest=TRUE )

        pvals$"p.adj.range.pos" <-
          pvals$"p.adj.range.neg" <-
          pvals$"p.adj.range"

        pvals[ which( coefs$"Coefficient" <= 0 ), "p.adj.range.pos" ] <- NA
        pvals[ which( coefs$"Coefficient" > 0 ), "p.adj.range.neg" ] <- NA

        pvals$"Significant_Coefficient_Sign" <-
          rep( x = NA, times = nrow( pvals ) )

        pvals[
          !is.na( pvals$"p.adj.range.pos" ), "Significant_Coefficient_Sign"
          ] <- "Positive"

        pvals[
          !is.na( pvals$"p.adj.range.neg" ), "Significant_Coefficient_Sign"
          ] <- "Negative"

        pvals$"Significant_Coefficient_Sign" <-
          factor(
            x = pvals$"Significant_Coefficient_Sign",
            levels = c( "Negative", "Positive" )
          )

        y <-
          data.frame( coefs,
                      pvals,
                      stringsAsFactors = FALSE,
                      check.names = FALSE )

        if ( !order.factor ) {

          y$"Factor" <-
            factor(
              x = y$"Factor",
              levels = colnames( x$"coefficients" )
            )

        }

      }

      tmp <- table( colnames( y ) )

      tmp <- names( tmp[ which( tmp > 1 ) ] )

      if ( length( tmp ) > 0 ) {

        idx.omitted <- NULL

        for ( i in 1:length( tmp ) ) {

          idx.omitted <-
            c( idx.omitted,
               which( colnames( y ) == tmp[ i ] )[ -1 ] )

        }

        y <- y[ , -idx.omitted ]

      }

      y$"Name.checked" <- make.names( y$"Name" )

      y <-
        merge(
          x = y,
          y = names.mapping,
          by = "Name.checked",
          all.x = FALSE, # TRUE,
          all.y = FALSE
          )

      if ( !is.null( omit.class ) ) {

        tmp <- which( y$"Class" %in% omit.class )

        if ( length( tmp ) > 0 ) {

          y <- y[ -tmp, ]

        }

      }

      if ( !is.null( omit.factor ) ) {

        tmp <- which( y$"Factor" %in% omit.factor )

        if ( length( tmp ) > 0 ) {

          y <- y[ -tmp, ]

        }

      }

      ### Heatmap

      if ( is.factor( y$"Factor" ) ) {

        factors.unique <- levels( y$"Factor" )

      } else {

        factors.unique <- unique( as.character( y$"Factor" ) )

      }

      if ( plot.individual ) {

        if ( order.factor ) {

          factors.unique <- sort( x = factors.unique, decreasing = FALSE )

        }

        if ( plot.all ) {

          figures <-
            vector(
              mode = "list",
              length = length( factors.unique ) + 1
            )

          names( figures ) <- c( factors.unique, "All" )

        } else {

          figures <-
            vector(
              mode = "list",
              length = length( factors.unique )
            )

          names( figures ) <- factors.unique

        }

      } else {

        figures <-
          vector(
            mode = "list",
            length = 1
          )

        names( figures ) <- "All"

      }

      if ( print.formula ) {

        if ( baseline.adjusted | survival ) {

          if ( !is.null( x$"parameters"$"formula.text" ) ) {

            printout.formula <-
              paste(
                as.character( x$"parameters"$"formula.text" ),
                collapse = " "
              )

            printout.formula <-
              paste(
                "Model:",
                printout.formula
              )

            printout.formula <-
              split_formula(
                x = printout.formula,
                width = formula.width
              )

          }

        } else if ( !is.null( x$"formula.text" ) ) {

          printout.formula <-
            paste(
              as.character( x$"formula.text" ),
              collapse = " "
            )

          printout.formula <-
            paste(
              "Model:",
              printout.formula
            )

          printout.formula <-
            split_formula(
              x = printout.formula,
              width = formula.width
            )

        }

      } else {

        printout.formula <- NULL

      }

      if ( plot.individual ) {

        # coefs <- sort( unique( x = y$"Factor", decreasing = FALSE ) )

        coefs <- factors.unique

        for ( i in 1:length( coefs ) ) {

          coef.i <- coefs[ i ]

          data.plot <- y[ y$"Factor" == coef.i, , drop = FALSE ]

          is.contrast <- grepl( x = coef.i, pattern = "\\s\\-\\s" )

          if ( is.contrast ) {

            coef.i <-
              stringr::str_replace(
                string = coef.i,
                pattern = "\\s\\-\\s",
                replacement = " vs. "
              )

          }

          if ( axis.x.carbons ) {

            plot.i <-
              ggplot2::ggplot(
                data = data.plot,
                mapping =
                  ggplot2::aes_string(
                    x = "N.carbons",
                    y = "N.double.bonds",
                    fill = "Coefficient"
                  )
              ) +
              ggplot2::ylab( label = "Number of fatty-acid double bonds" ) +
              ggplot2::xlab( label = "Number of fatty-acid carbon atoms" )

          } else {

            plot.i <-
              ggplot2::ggplot(
                data = data.plot,
                mapping =
                  ggplot2::aes_string(
                    x = "N.double.bonds",
                    y = "N.carbons",
                    fill = "Coefficient"
                  )
              ) +
              ggplot2::xlab( label = "Number of fatty-acid double bonds" ) +
              ggplot2::ylab( label = "Number of fatty-acid carbon atoms" )

          }

          if ( is.contrast ) {

            title.i <- paste( "Contrast:", coef.i )

          } else if ( survival ) {

            title.i <- paste( "Hazard ratio to", coef.i )

          } else {

            title.i <- paste( "Coefficient:", coef.i )

          }

          plot.i <- plot.i +
            ggplot2::theme_dark() +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient2(
              low = "blue",
              high = "red",
              midpoint = ifelse( test = survival, yes = 1, no = 0 )
            ) +
            ggplot2::ggtitle(
              label = title.i,
              subtitle = printout.formula
            ) +
            ggplot2::theme(
              legend.position = "top",
              legend.text = ggplot2::element_text( angle = 45 )
            ) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(
                angle = 45,
                vjust = 0.5
              )
            ) +
            ggplot2::theme(
              strip.text.y = ggplot2::element_text( angle = 0 )
            ) +
            ggplot2::scale_x_continuous(
              breaks = get_integer_breaks,
              expand = ggplot2::expand_scale( mult = 0, add = 0.5 )
            ) +
            ggplot2::scale_y_continuous(
              breaks = get_integer_breaks,
              expand = ggplot2::expand_scale( mult = 0, add = 0.5 )
            )

          if ( class.facet == "row" ) {

            plot.i <-
              plot.i +
              ggplot2::facet_grid(
                rows = ggplot2::vars( Class ),
                scales = scales,
                space = space
              )

          } else if ( class.facet == "col" ) {

            plot.i <-
              plot.i +
              ggplot2::facet_grid(
                cols = ggplot2::vars( Class ),
                scales = scales,
                space = space
              )

          } else {

            plot.i <-
              plot.i +
              ggplot2::facet_wrap(
                facets = ggplot2::vars( Class ),
                scales = scales
              )

          }

          if ( shadowtext ) {

            plot.i <-
              plot.i +
              shadowtext::geom_shadowtext(
                mapping = ggplot2::aes_string( label = "p.adj.range.neg" ),
                color = "blue",
                bg.colour = "white"
              ) +
              shadowtext::geom_shadowtext(
                mapping = ggplot2::aes_string( label = "p.adj.range.pos" ),
                color = "red",
                bg.colour = "white"
              )

          } else {

            plot.i <-
              plot.i +
              ggplot2::geom_point(
                mapping =
                  ggplot2::aes_string(
                    color = "Significant_Coefficient_Sign"
                  ),
                na.rm = TRUE,
                shape = 19
              ) +
              ggplot2::scale_color_manual(
                values = c( "blue", "red" ),
                drop = FALSE,
                na.translate = FALSE
              ) +
              ggplot2::geom_point(
                mapping =
                  ggplot2::aes_string(
                    shape = "P_Value"
                  ),
                color = "white",
                size = p.val.label.size
              ) +
              ggplot2::scale_shape_manual(
                values = p.val.labels,
                drop = FALSE,
                na.translate = FALSE
              ) +
              ggplot2::guides(
                color =
                  ggplot2::guide_legend(
                    order = 2,
                    override.aes =
                      list(
                        size = legend.key.size.multiplier * p.val.label.bg.size
                      ),
                    title = "Sign of\nSignificant\nCoefficient"
                  ),
                shape =
                  ggplot2::guide_legend(
                    order = 3,
                    override.aes =
                      list(
                        size = legend.key.size.multiplier * p.val.label.size
                      ),
                    title = "P-Value"
                  )
              )

            if ( survival ) {

              plot.i <-
                plot.i +
                ggplot2::guides(
                  fill =
                    ggplot2::guide_colorbar(
                      order = 1,
                      title = "Hazard\nRatio"
                    )
                )

            } else {

              plot.i <-
                plot.i +
                ggplot2::guides( fill = ggplot2::guide_colorbar( order = 1 ) )

            }

          }

          if ( class.facet == "wrap" ) {

            if (
              !is.null( range.min.N.carbons ) |
              !is.null( range.min.N.double.bonds )
            ) {

              data.blank <-
                create_blank_data_for_range(
                  x = y,
                  range.min.N.carbons = range.min.N.carbons,
                  range.min.N.double.bonds = range.min.N.double.bonds
                )

              plot.i <-
                plot.i +
                ggplot2::geom_blank( data = data.blank )

            }

          }

          figures[[ coef.i ]] <- plot.i

          if ( print.figure ) {

            print( plot.i )

          }

        }

      }

      if ( plot.all ) {

        data.plot <- y[ y$"Factor" != "(Intercept)", , drop = FALSE ]

        if ( wrap.contrast.name ) {

          levels( data.plot$"Factor" ) <-
            stringr::str_replace(
              string = levels( data.plot$"Factor" ),
              pattern = "\\s\\-\\s",
              replacement = "\nvs.\n"
            )

        } else {

          levels( data.plot$"Factor" ) <-
            stringr::str_replace(
              string = levels( data.plot$"Factor" ),
              pattern = "\\s\\-\\s",
              replacement = " vs. "
            )

        }

        if ( axis.x.carbons ) {

          plot.i <-
            ggplot2::ggplot(
              data = data.plot,
              mapping =
                ggplot2::aes_string(
                  x = "N.carbons",
                  y = "N.double.bonds",
                  fill = "Coefficient"
                )
            ) +
            ggplot2::ylab( label = "Number of fatty-acid double bonds" ) +
            ggplot2::xlab( label = "Number of fatty-acid carbon atoms" )

        } else {

          plot.i <-
            ggplot2::ggplot(
              data = data.plot,
              mapping =
                ggplot2::aes_string(
                  x = "N.double.bonds",
                  y = "N.carbons",
                  fill = "Coefficient"
                )
            ) +
            ggplot2::xlab( label = "Number of fatty-acid double bonds" ) +
            ggplot2::ylab( label = "Number of fatty-acid carbon atoms" )

        }

        plot.i <-
          plot.i +
          ggplot2::theme_dark() +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(
            low = "blue",
            high = "red"
          ) +
          ggplot2::theme(
            legend.position = "top",
            legend.text = ggplot2::element_text( angle = 45 )
          ) +
          ggplot2::theme(
            axis.text.x =
              ggplot2::element_text(
                angle = 45,
                vjust = 0.5
              )
          ) +
          ggplot2::theme( strip.text.y = ggplot2::element_text( angle = 0 ) ) +
          ggplot2::scale_x_continuous(
            breaks = get_integer_breaks,
            expand = ggplot2::expand_scale( mult = 0, add = 0.5 )
          ) +
          ggplot2::scale_y_continuous(
            breaks = get_integer_breaks,
            expand = ggplot2::expand_scale( mult = 0, add = 0.5 )
          )

        if ( shadowtext ) {

          plot.i <-
            plot.i +
            shadowtext::geom_shadowtext(
              mapping = ggplot2::aes_string( label = "p.adj.range.neg" ),
              color = "blue",
              bg.colour = "white",
              show.legend = TRUE
            ) +
            shadowtext::geom_shadowtext(
              mapping = ggplot2::aes_string( label = "p.adj.range.pos" ),
              color = "red",
              bg.colour = "white",
              show.legend = TRUE
            )

        } else {

          plot.i <-
            plot.i +
            ggplot2::geom_point(
              mapping =
                ggplot2::aes_string(
                  color = "Significant_Coefficient_Sign"
                ),
              na.rm = TRUE,
              shape = 19
            ) +
            ggplot2::scale_color_manual(
              values = c( "blue", "red" ),
              drop = FALSE,
              na.translate = FALSE
            ) +
            ggplot2::geom_point(
              mapping =
                ggplot2::aes_string(
                  shape = "P_Value"
                ),
              color = "white",
              size = p.val.label.size
            ) +
            ggplot2::scale_shape_manual(
              values = p.val.labels,
              drop = FALSE,
              na.translate = FALSE
            ) +
            ggplot2::guides(
              fill = ggplot2::guide_colorbar( order = 1 ),
              color =
                ggplot2::guide_legend(
                  order = 2,
                  override.aes =
                    list(
                      size = legend.key.size.multiplier * p.val.label.bg.size
                    ),
                  title = "Sign of\nSignificant\nCoefficient"
                ),
              shape =
                ggplot2::guide_legend(
                  order = 3,
                  override.aes =
                    list(
                      size = legend.key.size.multiplier * p.val.label.size
                    ),
                  title = "P-Value"
                )
            )

        }

        if ( print.formula ) {

          plot.i <-
            plot.i +
            ggplot2::ggtitle( label = printout.formula )

        }

        if ( class.facet == "col" ) {

          plot.i <-
            plot.i +
            ggplot2::facet_grid(
              rows = ggplot2::vars( Factor ),
              cols = ggplot2::vars( Class ),
              scales = scales,
              space = space
            )

        } else {

          plot.i <-
            plot.i +
            ggplot2::facet_grid(
              cols = ggplot2::vars( Factor ),
              rows = ggplot2::vars( Class ),
              scales = scales,
              space = space
            )

        }

        figures[[ "All" ]] <- plot.i

        if ( print.figure ) {

          print( plot.i )

        }

      }

    }

    return( figures )

  }
