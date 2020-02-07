#' Create 'lipidomeR' heatmaps of arbitrary lipid-specific values.
#'
#' Use this function to create a heatmap of any lipid-specific values.
#'    Note: Use the function \code{\link{heatmap_lipidome_from_limma}} to create
#'    heatmaps of model statistics.
#'
#' @param x (Required) named vector of numeric values to create a figure of.
#'    Names need to match to the argument names.mapping through
#'    the function \code{\link{map_lipid_names}}. Alternatively, a data frame
#'    can be supplied. In that case, set \code{melt.x = TRUE}.
#' @param names.mapping (Required) mapping of lipid names from
#'    the \code{\link{map_lipid_names}} function.
#' @param axis.x.carbons (Optional) \code{TRUE} or \code{FALSE}: Should
#'    the lipid size (i.e., number of carbon atoms in the fatty acid chain) be
#'    presented on the x-axis or y-axis?
#' @param class.facet (Optional) character string with possible values
#'    \code{'col'}, \code{'row'} or \code{'wrap'}:
#'    Present lipid classes as panels organized into columns, rows or into a
#'    wrapped layout spanning over multiple rows and columns. The alternative
#'    \code{'wrap'} is only available with \code{plot.infividual = TRUE}.
#' @param fill.direction (Optional) \code{TRUE} or \code{FALSE}: Should color
#'    fill be in an increasing direction?
#' @param fill.limits (Optional) numeric vector of length two, indicating
#'    the limits of the fill scale.
#' @param fill.midpoint (Optional) numeric value specifying the midpoint of
#'    the fill scale.
#' @param melt.value.name (Optional) character string, specifying the name of
#'    the variable that will be shown as fill in the heatmap.
#' @param melt.variable.name (Optional) character string, specifying the name of
#'    of the variable that will be used to creating faceted sub-heatmaps.
#' @param melt.x (Optional) \code{TRUE} or \code{FALSE}: Should the argument
#'    \code{x} be molten by the function \code{\link[reshape2]{melt}}
#'    prior to plotting? Set \code{melt.x = TRUE}, if you want to plot a data
#'    frame instead of a vector of values. In that case, each column of the
#'    data frame will be plotted as an individual facet.
#' @param range.min.N.carbons (Optional) numeric value to specify the minimum
#'    range of the axis showing the lipid size (number of carbon atoms in the
#'    fatty acid chains). This value can be increased from the default value to
#'    improve readability in situtions, where there are lipid classes with
#'    little or no variation in the lipid size.
#' @param range.min.N.double.bonds (Optional) numeric value to specify
#'    the minimum range of the axis showing the lipid saturation (number of
#'    double bonds in the fatty acid chains). This value can be increased from
#'    the default value to improve readability in situtions, where there are
#'    lipid classes with little or no variation in the lipid saturation.
#' @param scale.fill.log (Optional) numeric value specifying the base of
#'    the logarithm, which will be used to creating a logarithmic scale for
#'    the fill scale of the plot.
#' @param scales (Optional) character string with possible values
#'    \code{'fixed'}, \code{'free'}, \code{'free_x'} or \code{'free_y'}. This
#'    argument specifies, whether the axes in multiple sub-heatmaps will be in
#'    the same scale (\code{'fixed'}) or in a scale specific to each sub-figure.
#'    See the function \code{\link[ggplot2]{facet_grid}} for details.
#' @param space (Optional) character string with possible values
#'    \code{'fixed'}, \code{'free'}, \code{'free_x'} or \code{'free_y'}.
#'    This argument specifies, whether the sub-heatmaps will be of identical
#'    size (\code{'fixed'}) or not.
#' @param x.names (Optional) character string specifying the name of the
#'    variable in the argument \code{x}, which will be used to matching
#'    the values to the argument \code{names.mapping}.
#'    Use this argument only together with \code{melt.x = TRUE}.
#' @param x.variables (Optional) character vector specifying the names of
#'    the variables, which will be included as individual facets in
#'    the figure.
#'    Use this argument only together with \code{melt.x = TRUE}.
#' @export
heatmap_lipidome <-
  function(
    x,
    names.mapping,
    axis.x.carbons = TRUE,
    class.facet = "row",
    fill.direction = "increasing",
    fill.limits = c( 0, 40 ),
    fill.midpoint = 20,
    melt.value.name = "CV",
    melt.variable.name = NULL,
    melt.x = TRUE,
    range.min.N.carbons = 5,
    range.min.N.double.bonds = 5,
    scale.fill.log = NULL,
    scales = "free_y",
    space = "free",
    x.names = "row.names",
    x.variables = NULL
  ) {

    Class <- NULL # Avoid note with 'ggplot2::Vars()'.

    x.in <- x

    if ( is.null( melt.value.name ) ) {

      melt.value.name <- "value"

    }

    if ( !is.null( scale.fill.log ) ) {

      melt.value.name.log <-
        paste(
          melt.value.name,
          "log",
          scale.fill.log,
          sep = "_"
        )

    }

    if ( melt.x ) {

      if ( is.null( x.variables ) ) {

        x <-
          reshape2::melt(
            data = x,
            id.vars = x.names,
            value.name = melt.value.name
          )

      } else {

        x <-
          reshape2::melt(
            data = x,
            id.vars = x.names,
            measure.vars = x.variables,
            value.name = melt.value.name
          )

      }

    } else {

      if ( is.null( melt.variable.name ) ) {

        melt.variable.name <- "variable"
        x$"variable" <- rep( x = "", times = nrow( x ) )

      } else {

        colnames( x )[ colnames( x ) == melt.variable.name ] <- "variable"

      }
    }

    y <-
      merge(
        x = x,
        y = names.mapping,
        by.x = x.names,
        by.y = "Name"
      )

    if ( !is.null( scale.fill.log ) ) {

      y[ , melt.value.name ] <-
        log(
          x = y[ , melt.value.name ],
          base = scale.fill.log
        )

      colnames( y )[ colnames( y ) == melt.value.name ] <- melt.value.name.log

      melt.value.name <- melt.value.name.log

    }

    y$"truncated.above" <-
      y$"truncated.below" <-
      rep( x = NA, times = nrow( y ) )

    if ( !is.null( fill.limits ) ) {

      tmp <- which( y[ , melt.value.name ] < fill.limits[ 1 ] )
      y[ tmp, "truncated.below" ] <- "-"
      y[ tmp, melt.value.name ] <- fill.limits[ 1 ]

      tmp <- which( y[ , melt.value.name ] > fill.limits[ 2 ] )
      y[ tmp, "truncated.above" ] <- "+"
      y[ tmp, melt.value.name ] <- fill.limits[ 2 ]

    }

    if ( axis.x.carbons ) {

      plot <-
        ggplot2::ggplot(
          data = y,
          mapping =
            ggplot2::aes_string(
              x = "N.carbons",
              y = "N.double.bonds",
              fill = melt.value.name
            )
        ) +
        ggplot2::ylab( label = "Number of fatty-acid double bonds" ) +
        ggplot2::xlab( label = "Number of fatty-acid carbon atoms" )

    } else {

      plot <-
        ggplot2::ggplot(
          data = y,
          mapping =
            ggplot2::aes_string(
              x = "N.double.bonds",
              y = "N.carbons",
              fill = melt.value.name
            )
        ) +
        ggplot2::xlab( label = "Number of fatty-acid double bonds" ) +
        ggplot2::ylab( label = "Number of fatty-acid carbon atoms" )

    }

    plot <-
      plot +
      ggplot2::geom_tile() +
      ggplot2::theme_dark() +
      ggplot2::theme(
        legend.position = "top",
        legend.text = ggplot2::element_text( angle = 45 ) ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text( angle = 45, vjust = 0.5 )
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

    if ( class.facet == "row" ) {

      plot <-
        plot +
        ggplot2::facet_grid(
          rows = ggplot2::vars( Class ),
          scales = scales
          ,
          space = space
        )

    } else if ( class.facet == "col" ) {

      plot <-
        plot +
        ggplot2::facet_grid(
          cols = ggplot2::vars( Class ),
          scales = scales
          ,
          space = space
        )

    } else {

      plot <-
        plot +
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

        plot <-
          plot +
          ggplot2::geom_blank( data = data.blank )

      }

    }

    if ( is.numeric( y[ , melt.value.name ] ) ) {

      if ( fill.direction == "increasing" ) {

        plot <-
          plot +
          ggplot2::scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red",
            midpoint = fill.midpoint,
            limits = fill.limits,
            na.value = "black"
          ) +
          ggplot2::geom_text(
            mapping = ggplot2::aes_string( label = "truncated.below" ),
            color = "red" ) +
          ggplot2::geom_text(
            mapping = ggplot2::aes_string( label = "truncated.above" ),
            color = "blue" )

      } else {

        plot <-
          plot +
          ggplot2::scale_fill_gradient2(
            low = "red",
            mid = "white",
            high = "blue",
            midpoint = fill.midpoint,
            limits = fill.limits,
            na.value = "black"
          ) +
          ggplot2::geom_text(
            mapping = ggplot2::aes_string( label = "truncated.below" ),
            color = "blue" ) +
          ggplot2::geom_text(
            mapping = ggplot2::aes_string( label = "truncated.above" ),
            color = "red" )

      }

    } else {

      plot <-
        plot +
        ggplot2::scale_fill_discrete( na.value = "black" )

    }

    print( plot )

  }
