#' @importFrom BiocManager install
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 expand_scale
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_blank
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_discrete
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_dark
#' @importFrom ggplot2 vars
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom knitr kable
#' @importFrom limma contrasts.fit
#' @importFrom limma duplicateCorrelation
#' @importFrom limma eBayes
#' @importFrom limma lmFit
#' @importFrom limma topTable
#' @importFrom reshape2 melt
#' @importFrom shadowtext geom_shadowtext
#' @importFrom stats as.formula
#' @importFrom stats median
#' @importFrom stats model.matrix
#' @importFrom stats p.adjust
#' @importFrom stringr str_locate
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
#' @importFrom stringr str_sub
#' @importFrom stringr str_wrap
#' @importFrom tableone CreateTableOne
#' @importFrom tidyr gather



get_integer_breaks <-
  function(
    x,
    include.zero = FALSE
  ) {

    max <- ceiling( max( x = x, na.rm = TRUE ) )

    min <- floor( min( x = x, na.rm = TRUE ) )

    if ( include.zero ) {

      min <- min( 0, min )

    }

    x.range <- max - min

    if ( x.range > 15 ) {

      y <-
        base::pretty(
          x = c( min, max ),
          n = 4,
          min.n = 4
        )

    } else if ( x.range > 4 ) {

      y <-
        base::pretty(
          x = c( min, max ),
          n = 3,
          min.n = 3
        )

    } else {

      y <- round( x = median( x = c( min, max ), na.rm = TRUE ) )

    }

    names( y ) <- attr( y, "labels" )

    return( y )

  }



create_blank_data_for_range <-
  function(
    x,
    range.min.N.carbons,
    range.min.N.double.bonds
  ) {

    tmp <-
      by(
        data = x[ , c( "N.carbons", "N.double.bonds" ) ],
        INDICES = x$"Class",
        FUN =
          function( x ) {
            apply(
              X = x,
              MARGIN = 2,
              FUN =
                function( x ) {

                  sum( range( x = x, na.rm = TRUE ) ) / 2

                }
            )
          }
      )

    tmp <- simplify2array( x = tmp, higher = FALSE )

    y <- x[ NULL, ]

    if ( !is.null( range.min.N.carbons ) ) {

      y <-
        dplyr::bind_rows(
          y,
          data.frame(
            Class = colnames( tmp ),
            t(
              tmp[ "N.carbons", , drop = FALSE ] -
                range.min.N.carbons / 2
            ),
            stringsAsFactors = FALSE
          )

        )

      y <-
        dplyr::bind_rows(
          y,
          data.frame(
            Class = colnames( tmp ),
            t(
              tmp[ "N.carbons", , drop = FALSE ] +
                range.min.N.carbons / 2
            ),
            stringsAsFactors = FALSE
          )

        )

    }

    if ( !is.null( range.min.N.double.bonds ) ) {

      y <-
        dplyr::bind_rows(
          y,
          data.frame(
            Class = colnames( tmp ),
            t(
              tmp[ "N.double.bonds", , drop = FALSE ] -
                range.min.N.double.bonds / 2
            ),
            stringsAsFactors = FALSE
          )

        )

      y <-
        dplyr::bind_rows(
          y,
          data.frame(
            Class = colnames( tmp ),
            t(
              tmp[ "N.double.bonds", , drop = FALSE ] +
                range.min.N.double.bonds / 2
            ),
            stringsAsFactors = FALSE
          )

        )

    }

    return( y )

  }



split_formula <-
  function(
    x,
    width = 110,
    width.indent = 16
  ) {

    y <- x

    if ( nchar( x ) > width ) {

      y <-
        stringr::str_replace_all(
          string = y,
          pattern = " \\* ",
          replacement = "*"
        )

      lines.split <-
        stringr::str_wrap(
          string = y,
          width = width
        )

      lines.split <-
        stringr::str_split(
          string = lines.split,
          pattern = "\\\n"
        )

      lines.split <- lines.split[[ 1 ]]

      y.first <- lines.split[ 1 ]

      if ( length( lines.split ) > 1  ) {

        whitespace <-
          paste(
            rep( x = " ", times = width.indent ),
            collapse = ""
          )

        lines.remaining <- lines.split[ -1 ]

        lines.remaining <-
          paste(
            lines.remaining,
            collapse = " "
          )

        last.char <-
          substr(
            x = y.first,
            start = nchar( y.first ),
            stop = nchar( y.first )
          )

        if ( last.char %in% c( "+", "*", ":" ) ) {

          lines.remaining <-
            paste(
              last.char,
              lines.remaining,
              sep = " "
            )

        } else {

          y.first <-
            paste(
              y.first,
              substr(
                x = lines.remaining[ 1 ],
                start = 1,
                stop = 1
              ),
              sep = " "
            )
        }

        y.remaining <-
          stringr::str_wrap(
            string = lines.remaining,
            width = width - width.indent - 10
          )

        y.remaining <-
          stringr::str_replace(
            string = y.remaining,
            pattern = "^\\+\\\n",
            replacement =
              paste0(
                whitespace,
                "... + "
              )
          )

        y.remaining <-
          stringr::str_replace(
            string = y.remaining,
            pattern = "^\\+",
            replacement =
              paste0(
                whitespace,
                "... +"
              )
          )

        y.remaining <-
          stringr::str_replace_all(
            string = y.remaining,
            pattern = "\\\n\\+\\\n",
            replacement = " +\n"
          )

        y.remaining <-
          stringr::str_replace_all(
            string = y.remaining,
            pattern = "\\+\\\n",
            replacement =
              paste0(
                "+ ...\n",
                whitespace,
                "... + "
              )
          )

        y.remaining <-
          stringr::str_replace_all(
            string = y.remaining,
            pattern = "\\\n\\+",
            replacement =
              paste0(
                " + ...\n",
                whitespace,
                "... +"
              )
          )

        y <-
          paste(
            y.first,
            y.remaining,
            sep = " ...\n"
          )

      } else {

        y <- y.first

      }

      y <-
        stringr::str_replace_all(
          string = y,
          pattern = "\\*",
          replacement = " * "
        )

    }

    return( y )

  }
