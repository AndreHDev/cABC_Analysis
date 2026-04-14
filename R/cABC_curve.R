#' cABC Curve Computation
#'
#' Computes ABC curve (cumulative effort vs yield) with spline interpolation and derivatives.
#' 
#' @importFrom stats spline splinefun
#' @importFrom utils head tail
#'
#' @param Data Numeric vector/matrix. First column used if matrix. Only positive values used.
#' @param p Optional x-values for spline interpolation. Default: finer grid for large datasets.
#'
#' @return List containing:
#'   Curve: Data frame with Effor (x) and Yield (y) of interpolated curve
#'   Slope: Data frame with p (x-values) and dABC (first derivative)
#' 
#' @author: MT 11/2014
#' 1.Editor: MT 01/2015
#' 2.Editor: FL
#' 3.Editor: MT 11/2017: Doku neu
#' 4.Editor: AH 01/2026: Refactor, Spline changed to hyman for monotone property
#' @noRd
cABC_curve <- function(Data, p) {
  
  # Data cleaned before
  CleanedData <- Data
  
  # Set interpolation grid
  rows <- length(CleanedData)
  if(missing(p)) {
    if(rows < 101) {
      p <- seq(0, 1, by = 0.01)
    } else {
      p <- seq(0, 1, by = 0.001)
    }
  }
  
  # Compute empirical cumulative curve
  sorted <- sort(CleanedData, decreasing = TRUE, na.last = TRUE)
  Anteil <- sorted
  y <- cumsum(Anteil) / tail(cumsum(Anteil), 1)
  x <- (1:rows) / rows
  
  # Ensure curve passes through (0,0) and (1,1)
  if(head(y, 1) > 0) {
    x <- c(0, x)
    y <- c(0, y)
  }
  if(tail(x, 1) < 1) {
    x <- c(x, 1)
    y <- c(y, 1)
  }
  
  # Spline interpolation (using splinefun from stats)
  V <- spline(x, y, xout = p, method = "hyman")
  Effort <- V$x
  Yield <- V$y
  
  # Cap Yield at 1 (interpolation error handling)
  inds <- which(Yield >= 1)
  if(length(inds) > 0) {
    ind1 <- min(inds)
    if(ind1 < length(Yield)) {
      Yield[ind1:length(Yield)] <- 1
    }
  }
  
  # Compute first derivative using splinefun
  n <- length(Effort)
  curve_fun <- splinefun(Effort, Yield)
  dABC <- curve_fun((1:n)/n, deriv = 1)
  
  # Return structured result
  list(
    Curve = data.frame(Effort = Effort, Yield = Yield),
    Slope = data.frame(p = p, dABC = dABC)
  )
}
