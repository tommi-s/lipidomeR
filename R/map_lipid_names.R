#' Enumerate lipid names into values about lipid class, size and saturation
#'
#' Use this function to create a mapping of the lipids to values, which can be
#'    used to creating 'lipidomeR' heatmaps.
#'
#' @param x (Required) character vector of lipid names. The names are expected
#'    in the format 'XY(C:D)', where 'XY' is the abbreviation of the lipid
#'    class, 'C' is the total number of carbon atoms in the fatty-acid
#'    chains, and 'D' is the total number of double-bonds in the fatty
#'    acid chains.
#'
#' @return Data frame with lipid names in various formats for R and
#'    the enumerated values of lipid class (\code{Class}), lipid size
#'    (\code{N.carbons}) and lipid unsaturation (\code{N.double.bonds}).
#'
#' @seealso \code{\link{heatmap_lipidome_from_limma}} for creating 'lipidomeR'
#'    heatmaps of regression model results based on the output of this function.
#' @seealso \code{\link{heatmap_lipidome}} for creating 'lipidomeR' heatmaps
#'    of any lipid-specific values based on the output of this function.
#'
#' @export
#'
map_lipid_names <-
  function( x ) {

  colnames.y <-
    c(
      "Name",
      "Name.checked",
      "Name.clean",
      "Name.clean.checked",
      "Name.simple",
      "Name.simple.checked",
      "Class",
      "N.carbons",
      "N.double.bonds"
    )

  y <-
    array(
      dim =
        c(
          length( x ),
          length( colnames.y )
        )
    )

  y <- as.data.frame( y )

  colnames( y ) <- colnames.y

  y[ , "N.double.bonds" ] <- 0

  y$"Name" <-
    y$"Name.clean" <-
    x

  y$"Name.checked" <- make.names( y$"Name" )

  rownames( y ) <- y$"Name.checked"

  # For instance, "SM(d42:0)" becomes "SM(42:0)".

  y$"Name.clean" <-
    stringr::str_replace_all(
      string = y$"Name.clean",
      pattern = "\\([det]{0,1}[\\-]{0,1}",
      replacement = "("
    ) # TODO: e as O...

  y$"Name.clean" <-
    stringr::str_replace(
      string = y$"Name.clean",
      pattern = "\\(9Z\\)",
      replacement = ""
    )

  y$"Name.clean" <-
    stringr::str_replace(
      string = y$"Name.clean",
      pattern = "PC\\(O-",
      replacement = "PC-O/P("
    )

  y$"Name.clean" <-
    stringr::str_replace(
      string = y$"Name.clean",
      pattern = "PC\\(P-",
      replacement = "PC-O/P("
    )

  y$"Name.clean" <-
    stringr::str_replace(
      string = y$"Name.clean",
      pattern = "PE\\(O-",
      replacement = "PC-O/P("
    )

  y$"Name.clean" <-
    stringr::str_replace(
      string = y$"Name.clean",
      pattern = "PE\\(P-",
      replacement = "PC-O/P("
    )

  # Remove the odd character(s) preceding the carbon number.

  # For instance, "PC(16:0/O-1:0)" becomes "PC(16:0/1:0)".
  y$"Name.clean" <-
    stringr::str_replace_all(
      string = y$"Name.clean",
      pattern = "/[A-Za-z]{0,1}[\\-]{0,1}",
      replacement = "/"
    ) # TODO: FIX

    # For instance, "PC(38:2p)" becomes "PC-O/P(38:2)".

    tmp <-
      stringr::str_locate(
        string = y$"Name.clean",
        pattern = "\\([0-9]+\\:[0-9]+(e|p)"
      )

    tmp2 <- which( !is.na( tmp[ , "start" ] ) )

    if ( length( tmp2 ) > 0 ) {

      y[ tmp2, "Name.clean" ] <-
        paste0(
          stringr::str_sub(
            string = y[ tmp2, "Name.clean" ],
            start = 1,
            end = tmp[ tmp2, "start" ] - 1
          ),
          "-O/P",
          stringr::str_sub(
            string = y[ tmp2, "Name.clean" ],
            start = tmp[ tmp2, "start" ],
            end = tmp[ tmp2, "end" ] - 1
          ),
          stringr::str_sub(
            string = y[ tmp2, "Name.clean" ],
            start = tmp[ tmp2, "end" ] + 1
          )
        )

    }

  # properties <- list()
  # properties$N.carbons <- rep( x=NA, times=length( x ) )
  # properties$N.double.bonds <- rep( x=NA, times=length( x ) )

    ## Algorithm:
    #
    # 1) Find the first parentheses:
    #     start <- "("
    #     end <- ")"
    # 2) Find all ":"s between the parentheses.
    # 3) For each ":":
    #     Convert the sequence between "(" or "/" and ":" to a number.
    # 4) Sum up the numbers.

    position <- list()

    # Find the first opening parenthesis of the number sequence.

    position$"start" <-
      stringr::str_locate(
        string = y$"Name.clean",
        pattern = "\\([0-9]+:"
      )

    # Find all the middle dividers separating the number pairs from other pairs.

    position$"middle" <-
      stringr::str_locate_all(
        string = y$"Name.clean",
        pattern = "/[0-9]+:[0-9]+"
      )

    # Find the first closing parenthesis of the number sequence.

    position$"end" <-
      stringr::str_locate(
        string = y$"Name.clean",
        pattern = "[0-9]\\)"
      ) + 1

    # Find all the middle dividers separating the pair of numbers from each other.

    position$"colon" <-
      stringr::str_locate_all(
        string = y$"Name.clean",
        pattern = ":[0-9]"
      )

    # Initialize the vectors for the result.
    # properties.TG <- list()
    # properties.TG$N.carbons <- rep( x=NA, times=length( idx.TGs ) )
    # properties.TG$N.double.bonds <- rep( x=0, times=length( idx.TGs ) )

    for ( i in 1:nrow( y ) ) { # Go through all items.

      colons.before.end.i <-
        which( position$"colon"[[ i ]][ , 1 ] < position$"end"[ i, 1 ] )

      for ( j in 1:length( colons.before.end.i ) ) {

        # The first carbon number after the start of the description.
        if ( j == 1 ) {

          y[ i, "N.carbons" ] <-
            as.numeric(
              base::substr(
                x = y[ i, "Name.clean" ],
                start = position$"start"[ i, 1 ] + 1,
                stop = position$"colon"[[ i ]][ 1, 1 ] - 1
              )
            )

        } else { # Carbon numbers following the first one.

          y[ i, "N.carbons" ] <-
            y[ i, "N.carbons" ] +
            as.numeric(
              base::substr(
                x = y[ i, "Name.clean" ],
                start = position$"middle"[[ i ]][ j - 1, 1 ] + 1,
                stop = position$"colon"[[ i ]][ j, 1 ]-1
              )
            )

        }

        # Double bond number before the last number.
        if ( j < length( colons.before.end.i ) ) {

          y[ i, "N.double.bonds" ] <-
            y[ i, "N.double.bonds" ] +
            as.numeric(
              base::substr(
                x = y[ i, "Name.clean" ],
                start = position$"colon"[[ i ]][ j, 1 ] + 1,
                stop = position$"middle"[[ i ]][ j, 1 ] - 1
              )
            )

        } else { # The last double bond number.

          y[ i, "N.double.bonds" ] <-
            y[ i, "N.double.bonds" ] +
            as.numeric(
              base::substr(
                x = y[ i, "Name.clean" ],
                start = position$"colon"[[ i ]][ j, 1 ] + 1,
                stop = position$"end"[ i, 1 ] - 1
              )
            )

        }
      }
    }

    y$"Name.clean.checked" <-
      make.names( y$"Name.clean" )

    y$"Class" <-
      stringr::str_sub(
        string = y$"Name.clean",
        start = 1,
        end = position$"start"[ , 1 ] - 1
      )

    y$"Name.simple" <-
      paste(
        y$"Class",
        "(",
        y$"N.carbons",
        ":",
        y$"N.double.bonds",
        stringr::str_sub(
          string = y$"Name.clean",
          start = position$"end"[ , 1 ],
          end = nchar( y$"Name.clean" )
        ),
        sep = ""
      )

    y$"Name.simple.checked" <- make.names( y$"Name.simple" )

    return( y )

}
