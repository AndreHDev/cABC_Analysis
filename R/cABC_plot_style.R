#' ABC Plot style
#'
#' Constructs a validated list of cosmetic options for customizing ABC
#' plots produced by \code{\link{cABC_analysis}}. Pass the result to
#' \code{cABC_analysis(..., plot_args = cABC_plot_style(...))}.
#'
#' Any option left at \code{NULL} falls back to the package default 
#' \code{\link{cABC_default_plot_args}}.
#'
#' @param LineWidth Numeric. Line width of the main ABC curve. Default
#'   is \code{1.25}.
#' @param ShowUniform Logical. Whether to draw the uniform reference
#'   curve. Default is \code{TRUE}.
#' @param ShowBoundary Logical. Whether to draw lines indicating
#'   boundaries. Default is \code{TRUE}.
#' @param Plot_title Character. Title of the plot. Default is
#'   \code{"ABC plot"}.
#' @param UniformColor Color of the uniform reference curve. Any valid R
#'   color specification, e.g. a color name (\code{"green"}), a hex code
#'   (\code{"#00FF00"}), or a palette index (\code{colors()[26]}).
#'   Default is \code{"green"}.
#' @param IdentityColor Color of the identity line (y = x). Any valid R
#'   color specification (see \code{UniformColor}). Default is
#'   \code{grDevices::colors()[452]}.
#' @param EquilibriumColor Color of the equilibrium diagonal (y = 1 - x).
#'   Any valid R color specification (see \code{UniformColor}). Default
#'   is \code{grDevices::colors()[175]}.
#' @param CurveColor Color of the ABC curve itself. Any valid R color
#'   specification (see \code{UniformColor}). Default is \code{"blue"}.
#' @param PointColor Color of individual data points (shown when n < 20).
#'   Any valid R color specification (see \code{UniformColor}). Default
#'   is \code{"blue"}.
#' @param ABoundaryColor Color of the A|B boundary star and label. Any
#'   valid R color specification (see \code{UniformColor}). Default is
#'   \code{"red"}.
#' @param BBoundaryColor Color of the B boundary star. Any valid R color
#'   specification (see \code{UniformColor}). Default is \code{"green"}.
#' @param CBoundaryColor Color of the B|C boundary star and label. Any
#'   valid R color specification (see \code{UniformColor}). Default is
#'   \code{"blue"}.
#' @param BoundaryLineColor Color of the orthogonal boundary guide lines.
#'   Any valid R color specification (see \code{UniformColor}). Default
#'   is \code{"red"}.
#' @param LabelColor Color of the A/B/C set labels and counts. Any valid
#'   R color specification (see \code{UniformColor}). Default is
#'   \code{"red"}.
#' @param LegendTextSize Numeric. Font size of the legend text. Default
#'   is \code{4}.
#' @param LegendY Numeric, or \code{NULL} (default). Y-position where
#'   the legend starts. If \code{NULL}, this is calculated automatically
#'   based on the A/C boundary positions to avoid overlapping the plot.
#' @param LegendX Numeric. X-position of the legend. Default is
#'   \code{0.80}.
#' @param LegendSpacing Numeric. Vertical spacing between legend entries.
#'   Default is \code{0.05}.
#' @param Theme A ggplot2 theme object (e.g. \code{ggplot2::theme_minimal()}).
#'   Default is \code{ggplot2::theme_light()} with a few structural
#'   tweaks applied on top (no panel grid, a bordered panel, and a
#'   centered plot title).
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
                              LegendX           = NULL,
                              LegendY           = NULL,
                              LegendSpacing     = NULL,
                              Theme             = NULL) {
  
  # Collect only the options the user actually supplied
  supplied <- as.list(environment())
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]
  
  # Merge onto defaults + warn on anything unexpected (defensive; formals
  # already restrict names
  cABC_resolve_plot_args(supplied)
}