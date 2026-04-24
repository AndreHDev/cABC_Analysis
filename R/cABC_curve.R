#' cABC Curve Computation
#'
#' Computes cumulative percentage of largest data (effort) 
#' and cumulative percentages of sum of largest Data (yield) with monotone hyman 
#' spline interpolation used to generate in between points.
#' 
#' @importFrom stats spline splinefun
#' @importFrom utils head tail
#'
#' @param Data Numeric vector/matrix. First column used if matrix. Only positive values used.
#' @param p Optional x-values for spline interpolation. Default: finer grid for large datasets.
#'
#' @return List containing:
#'   Curve: Data frame with Effort (x) and Yield (y) of interpolated curve
#'   Slope: Data frame with p (x-values) and cABC (first derivative)
#' 
#' @keywords internal
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
  portion <- sorted
  y <- cumsum(portion) / tail(cumsum(portion), 1)
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
  # Use hyman as monotone spline
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
