#' Default ABC Plot Options
#'
#' Internal helper returning the default list of ABC plot customization
#' options. This is the single source of truth for option names and
#' defaults, used both by \code{\link{cABC_plot_style}} and by
#' \code{\link{cABC_resolve_plot_args}}.
#'
#' @return A named list of default plot options.
#' @keywords internal
cABC_default_plot_args <- function() {
  list(
    LineWidth         = 1.25,
    ShowUniform       = TRUE,
    ShowBoundary      = TRUE,
    Plot_title        = 'ABC plot',
    UniformColor      = 'green',
    IdentityColor     = grDevices::colors()[452],
    EquilibriumColor  = grDevices::colors()[175],
    CurveColor        = 'blue',
    PointColor        = 'blue',
    ABoundaryColor    = 'red',
    BBoundaryColor    = 'green',
    CBoundaryColor    = 'blue',
    BoundaryLineColor = 'red',
    LabelColor        = 'red',
    LegendTextSize    = 4,
    Theme = ggplot2::theme_light() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = 'grey', fill = NA, linewidth = 1),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
  )
}

#' Resolve ABC Plot Arguments
#'
#' Internal helper that merges a user-supplied plot options list with the
#' package defaults, filling in anything the user did not specify and
#' warning about any names it does not recognize.
#'
#' @param plot_args \code{NULL}, or a named list of plot options (as
#'   returned by \code{\link{cABC_plot_style}}, or a raw named list with
#'   the same element names).
#'
#' @return A complete named list of plot options, with every default
#'   element present.
#' @keywords internal
cABC_resolve_plot_args <- function(plot_args = NULL) {
  defaults <- cABC_default_plot_args()
  
  if (is.null(plot_args)) return(defaults)
  
  if (!is.list(plot_args)) {
    stop("'plot_args' must be a list or NULL, see ?cABC_plot_style")
  }
  
  supplied_names <- names(plot_args)
  unknown <- setdiff(supplied_names, names(defaults))
  if (length(unknown) > 0) {
    warning(sprintf(
      "Ignoring unrecognized plot_args element(s): %s. See ?cABC_plot_style for valid options.",
      paste(unknown, collapse = ", ")
    ))
    plot_args <- plot_args[setdiff(supplied_names, unknown)]
  }
  
  # Merge: user-supplied values override defaults, everything else keeps default
  merged <- utils::modifyList(defaults, plot_args)
  return(merged)
}