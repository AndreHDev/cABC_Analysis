#' ABC Classification
#'
#' Divides a numeric dataset into three classes (A, B, and C) using
#' ABC analysis. The classification is based on geometric properties
#' of the ABC curve and identifies regions of high, balanced, and
#' low efficiency.
#' Class interpretation:
#' \tabular{ll}{
#'   A: \tab Low effort, high yield (Pareto items) \cr
#'   B: \tab Balanced effort and yield \cr
#'   C: \tab High effort, low yield (submarginal items)
#' }
#' 
#' @param Data Positive numeric vector which is not uniformly distributed.
#'   If matrix or dataframe then the first column will be used.
#'
#' @param PlotIt Logical. If \code{TRUE}, an ABC plot is generated.
#'
#' @param useGGPlot Logical, default \code{TRUE}. If \code{TRUE} a ggplot2
#'   plot is produced; if \code{FALSE} a base-R plot is produced. Only
#'   relevant when \code{PlotIt = TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{Aind, Bind, Cind}{Integer vectors of indices (into the original
#'     \code{Data}) for items assigned to classes A, B, and C respectively.
#'     In special-case returns (single point or all-identical), only
#'     \code{Aind} is populated; \code{Bind} and \code{Cind} are
#'     \code{integer(0)}.}
#'   \item{ABexchanged}{Logical; \code{TRUE} if the Pareto point and
#'     Break-even point were swapped to maintain coordinate logic (i.e. the
#'     Break-even point was to the left of the Pareto point on the curve).}
#'   \item{A, B, C}{\code{c(x, y)} coordinates for the Pareto point (A),
#'     the Break-even point (B), and the Submarginal point (C).
#'     \code{NULL} in special-case returns.}
#'   \item{smallestAData}{Cumulative yield at the boundary of Class A.
#'     \code{NULL} in special-case returns.}
#'   \item{smallestBData}{Cumulative yield at the boundary of Class B.
#'     \code{NULL} in special-case returns.}
#'   \item{AlimitIndInInterpolation}{Index of the A boundary in the
#'     interpolated \code{[p, ABC]} curve. \code{NULL} in special-case
#'     returns.}
#'   \item{BlimitIndInInterpolation}{Index of the C boundary in the
#'     interpolated \code{[p, ABC]} curve. \code{NULL} in special-case
#'     returns.}
#'   \item{p}{Numeric vector of effort values (x-axis) of the interpolation
#'     curve. \code{NULL} in special-case returns.}
#'   \item{ABC}{Numeric vector of yield values (y-axis) of the interpolation
#'     curve. \code{NULL} in special-case returns.}
#'   \item{ABLimit}{Data value closest to the threshold separating Class A
#'     from Class B. \code{NULL} in special-case returns.}
#'   \item{BCLimit}{Data value closest to the threshold separating Class B
#'     from Class C. \code{NULL} in special-case returns.}
#' }
#' 
#' @details
#' Data cleaning: Before classification, non-numeric values and
#' \code{NA}s are coerced to \code{0}, negative values are set to \code{0}.
#' A warning is issued when items are dropped. If a matrix or data frame is 
#' supplied, only the first column is used.
#'
#' Degenerate inputs (single point, all-identical values, very small datasets)
#' are caught before curve fitting, see \code{\link{cABC_handle_specials}} for
#' the full behaviour. Boundary duplicate values that span two classes after
#' classification are resolved by \code{\link{cABC_postprocess_classes}}.
#' In both cases a warning is issued when a special case is triggered.
#' 
#' @examples
#' data("SwissInhabitants")
#' abc <- cABC_analysis(SwissInhabitants, PlotIt = TRUE)
#'
#' # Extract the data belonging to each class
#' A <- abc$Aind; B <- abc$Bind; C <- abc$Cind
#' Agroup <- SwissInhabitants[A]
#' Bgroup <- SwissInhabitants[B]
#' Cgroup <- SwissInhabitants[C]
#'
#' @author André Himmelspach (01/2026)
#' @export
cABC_analysis <- function(Data, PlotIt=FALSE, useGGPlot=TRUE) {
  
  # === DATA CLEANING WITH NAME PRESERVATION ===
  Data_orig_names <- names(Data)
  if(!is.vector(Data)){
    if(is.matrix(Data) || is.data.frame(Data)){
      if(ncol(Data) > 1) {
        warning('Using only first column of data')
        Data <- as.vector(Data[,1])
        Data_orig_names <- rownames(Data)[,1]
      } else Data <- as.vector(Data)
    }
  }

  # Force numeric, negative values are set to 0, remove zeros
  # Suppressing the NA warning, as NAs are set 0, about which the user is
  # later warned about
  Data <- suppressWarnings(as.numeric(unname(Data)))
  Data[!is.finite(Data)] <- 0
  Data[Data < 0] <- 0
  
  n_original <- length(Data)
  n_used <- sum(Data > 0)
  
  if(n_used == 0) stop("No positive values remain after cleaning")
  
  if(n_used < n_original) {
    warning(sprintf('Only %d of %d items are larger then 0.', n_used, n_original))
  }
  
  # PRESERVE NAMES
  if(is.null(Data_orig_names) || length(Data_orig_names) != length(Data)) {
    Data_orig_names <- paste0("Item", seq_along(Data))
  }
  names(Data) <- Data_orig_names
  
  # === SPECIAL CASE HANDLING ===
  special_case <- cABC_handle_specials(Data)
  if (!is.null(special_case)) return(special_case)
  
  # === ABC CURVE ===
  ABCcurvedata <- suppressWarnings(cABC_curve(Data))         
  curve <- ABCcurvedata$Curve
  Effort <- curve[, 'Effort']
  Yield <- curve[, 'Yield']
  
  # === ABC BOUNDARY POINTS ===
  curve_matrix <- cbind(Effort, Yield)
  # Find Point A (Pareto Point): 
  # This is the point on the curve closest to the ideal top-left corner (0,1)
  dist_pareto <- numeric(length(Effort))
  point <- c(0,1)
  for(i in 1:length(Effort)) {
    dist_pareto[i] <- sum((point - curve[i, ])^2)
  }

  A_idx <- which.min(dist_pareto)
  A_point <- curve_matrix[A_idx, ]
  
  # Find Point B (Break-even Point):
  # Identify where the curve's slope is 1 (where marginal effort equals marginal yield)
  slope <- ABCcurvedata$Slope[, 'dABC']
  B_idx_candidates <- which.min(abs(slope - 1))
  B_idx <- max(B_idx_candidates)
  # If multiple points exist, take the rightmost one
  B_point_candidate <- curve_matrix[B_idx, ]
  
  # Check if Point A is bigger then point B, if not switch them around to keep logic
  ABexchanged <- Effort[B_idx] < Effort[A_idx]
  if(ABexchanged) {
    B_point <- A_point
    A_point <- B_point_candidate
  } else {
    B_point <- B_point_candidate
  }
  
  # Find the submarginal point
  # This is the point closest to the "Juren" point (Bx, 1).
  dist_B_line <- colSums((t(curve_matrix) - c(B_point[1], 1))^2)
  C_idx <- which.min(dist_B_line)
  C_point <- curve_matrix[C_idx, ]
  
  # === CLASSIFICATION ===
  sorted_data <- sort(Data, decreasing = TRUE)
  
  # Points on ABC line
  x_vals <- (1:length(sorted_data)) / length(sorted_data)
  y_vals <- cumsum(sorted_data) / sum(sorted_data)
  
  # Determine what points belong to what group by comparing x values of given 
  # data (Effort) with specific limit points
  A_mask <- x_vals <= A_point[1]
  C_mask <- x_vals <= C_point[1]
  
  # Apply masks to data
  sorted_indices <- order(Data, decreasing = TRUE)
  Aind <- sorted_indices[A_mask]
  Bind <- sorted_indices[C_mask & !A_mask]
  Cind <- sorted_indices[!C_mask]
  
  # For backward compatibility and edge case handling also calculate:
  ABLimit <- sort(Data, decreasing = T)[round(A_point[1]*length(Data))]
  BCLimit <- sort(Data, decreasing = T)[round(C_point[1]*length(Data))]
  
  # === SPECIAL CASE HANDLING REGARDING SETS ===
  processed <- cABC_postprocess_classes(Aind, Bind, Cind, Data, sorted_data, ABLimit, BCLimit)
  Aind <- processed$Aind
  Bind <- processed$Bind
  Cind <- processed$Cind
  
  # Preserve names
  names(Aind) <- names(Data)[Aind]
  names(Bind) <- names(Data)[Bind]
  names(Cind) <- names(Data)[Cind]
  
  # Calculate class sizes for plot
  nA <- length(Aind)
  nB <- length(Bind)
  nC <- length(Cind)
  
  # Preserve names
  names(Aind) <- names(Data)[Aind]
  names(Bind) <- names(Data)[Bind]
  names(Cind) <- names(Data)[Cind]
  
  # === PLOT (if needed) ===
  if(PlotIt) {
    bounds <- list(A = A_point, B = B_point, C = C_point)
    set_counts <- list(nA = nA, nB = nB, nC = nC)
    if(useGGPlot){
      abc <- cABC_plotGG(ABCcurvedata, Data, Boundaries = bounds, set_counts, x_vals, y_vals,
                         ShowUniform=TRUE)
      print(abc)
    }else{
      abc <- cABC_plot(ABCcurvedata, Data, Boundaries = bounds, set_counts, x_vals, y_vals,
                       ShowUniform=TRUE)
    }
  }
  
  list(
    Aind = Aind, Bind = Bind, Cind = Cind,
    ABexchanged = ABexchanged,
    A = A_point, B = B_point, C = C_point,
    smallestAData = if(ABexchanged) Yield[B_idx] else Yield[A_idx],
    smallestBData = Yield[C_idx],
    AlimitIndInInterpolation = if(ABexchanged) B_idx else A_idx,
    BlimitIndInInterpolation = C_idx,
    p = Effort,
    ABC = Yield,
    ABLimit = ABLimit, 
    BCLimit = BCLimit
  )
}
