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
          reshape2::melt.data.frame(
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
            mapping = ggplot2::aes( label = truncated.below ),
            color = "red" ) +
          ggplot2::geom_text(
            mapping = ggplot2::aes( label = truncated.above ),
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
            mapping = ggplot2::aes( label = truncated.below ),
            color = "blue" ) +
          ggplot2::geom_text(
            mapping = ggplot2::aes( label = truncated.above ),
            color = "red" )

      }

    } else {

      plot <-
        plot +
        ggplot2::scale_fill_discrete( na.value = "black" )

    }

    print( plot )

  }
