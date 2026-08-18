#' ABC Plot style
#'
#' Constructs a validated list of cosmetic options for customizing ABC
#' plots produced by \code{\link{cABC_analysis}}. Pass the result to
#' \code{cABC_analysis(..., plot_args = cABC_plot_control(...))}.
#'
#' Any option left at \code{NULL} falls back to the package default.
#'
#' @param LineWidth Numeric. Line width of the main ABC curve.
#' @param ShowUniform Logical. Whether to draw the uniform reference curve.
#' @param ShowBoundary Logical. Whether to draw lines indicating boundaries.
#' @param Plot_title Character. Title of the plot.
#' @param UniformColor Color of the uniform reference curve. Any valid R
#'   color specification, e.g. a color name (\code{"green"}), a hex code
#'   (\code{"#00FF00"}), or a palette index (\code{colors()[26]}).
#' @param IdentityColor Color of the identity line (y = x). Any valid R
#'   color specification (see \code{UniformColor}).
#' @param EquilibriumColor Color of the equilibrium diagonal (y = 1 - x).
#'   Any valid R color specification (see \code{UniformColor}).
#' @param CurveColor Color of the ABC curve itself. Any valid R color
#'   specification (see \code{UniformColor}).
#' @param PointColor Color of individual data points (shown when n < 20).
#'   Any valid R color specification (see \code{UniformColor}).
#' @param ABoundaryColor Color of the A|B boundary star and label. Any
#'   valid R color specification (see \code{UniformColor}).
#' @param BBoundaryColor Color of the B boundary star. Any valid R color
#'   specification (see \code{UniformColor}).
#' @param CBoundaryColor Color of the B|C boundary star and label. Any
#'   valid R color specification (see \code{UniformColor}).
#' @param BoundaryLineColor Color of the orthogonal boundary guide lines.
#'   Any valid R color specification (see \code{UniformColor}).
#' @param LabelColor Color of the A/B/C set labels and counts. Any valid
#'   R color specification (see \code{UniformColor}).
#' @param LegendTextSize Numeric. Font size of the legend text.
#' @param Theme A ggplot2 theme object (e.g. \code{ggplot2::theme_minimal()}).
#'  Default is \code{ggplot2::theme_light()} with slight custom tweaks.
#'  Default is \code{ggplot2::theme_light()} with slight custom tweaks.
#'
#' @return A named list of plot options, suitable for the \code{plot_args}
#'   argument of \code{\link{cABC_analysis}}.
#'
#' @examples
#' opts <- cABC_plot_style(LineWidth = 2, CurveColor = "darkred")
#' data("SwissInhabitants")
#' abc <- cABC_analysis(SwissInhabitants, PlotIt = TRUE, plotArgs = opts)
#'
#' @seealso \code{\link{cABC_analysis}}
#' @export
cABC_plot_style <- function(LineWidth         = NULL,
                              ShowUniform       = NULL,
                              ShowBoundary      = NULL,
                              Plot_title        = NULL,
                              UniformColor      = NULL,
                              IdentityColor     = NULL,
                              EquilibriumColor  = NULL,
                              CurveColor        = NULL,
                              PointColor        = NULL,
                              ABoundaryColor    = NULL,
                              BBoundaryColor    = NULL,
                              CBoundaryColor    = NULL,
                              BoundaryLineColor = NULL,
                              LabelColor        = NULL,
                              LegendTextSize    = NULL,
                              Theme             = NULL) {
  
  # Collect only the options the user actually supplied
  supplied <- as.list(environment())
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]
  
  # Merge onto defaults + warn on anything unexpected (defensive; formals
  # already restrict names
  cABC_resolve_plot_args(supplied)
}