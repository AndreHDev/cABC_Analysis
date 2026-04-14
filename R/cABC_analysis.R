#' ABC Classification
#'
#' Divides a numeric dataset into three classes (A, B, and C) using
#' ABC analysis. The classification is based on geometric properties
#' of the ABC curve and identifies regions of high, balanced, and
#' low efficiency.
#'
#' Class interpretation:
#' A = low effort, high yield 
#' B = balanced effort and yield
#' C = high effort, low yield 
#'
#' @param Data Positive numeric vector which is not uniformly distributed.
#'   If matrix or dataframe then the first column will be used.
#'
#' @param ABCcurvedata Legacy: only for internal usage, Optional list returned by 
#'    ABCcurve(), containing the interpolated ABC curve.
#'    
#' @param PlotIt Logical. If TRUE, an ABC plot is generated.
#' 
#' @param useGGPlot Logical, default TRUE. Sets if base R or gg plot should be used
#' 
#' @return A list containing: 
#'  Aind, Bind, Cind: Indices of the items partitioned into classes A, B, and C.
#'    If the A class is empty, the B class is set as the A class.
#'  ABexchanged: Logical; TRUE if the Pareto point and Break-even point were 
#'    swapped to maintain coordinate logic.
#'  A, B, C: The c(x, y) coordinates for the Pareto point (A), 
#'    the Break-even point (B), and the Submarginal point (C).
#'  smallestAData: The cumulative yield at the boundary of Class A.
#'  smallestBData: The cumulative yield at the boundary of Class B.
#'  AlimitIndInInterpolation: The index of the A boundary in the [p, ABC] curve.
#'  BlimitIndInInterpolation: The index of the C boundary in the [p, ABC] curve.
#'  p: Numeric vector of the effort (x-axis) of the interpolation curve.
#'  ABC: Numeric vector of the yield (y-axis) of the interpolation curve.
#'  ABLimit: The data point closes to the value threshold separating Class A from Class B.
#'  BCLimit: The data point closes to the value threshold separating Class B from Class C.
#'
#' @author AH (01/2026): Refactor, Group assignment fixed, monotone Spline
#' @export
cABC_analysis <- function(Data, ABCcurvedata, PlotIt=FALSE, useGGPlot=TRUE) {
  
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
    warning(sprintf('Only %d of %d items are positive.', n_used, n_original))
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
