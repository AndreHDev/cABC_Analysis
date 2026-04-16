#' cABC_plotGG - ggplot2 version matching base R cABC_plot
#'
#' @importFrom grDevices colors
#' @param CurveData Data about the ABC Curve as returned by ABC_curve
#' 
#' @param CleanData Clean original input data.
#' 
#' @param Boundaries A list with numeric vectors A, B, and C,
#'   each of length 2, giving the x/y coordinates of the ABC boundaries.
#'   
#' @param Set_counts A list with elements nA, nB, nC giving the
#'   number of observations in sets A, B, and C.
#'   
#' @param x_vals Numeric vector of x coordinates of original data points.
#' 
#' @param y_vals Numeric vector of y coordinates of original data points.
#' 
#' @param LineWidth Numeric. Line width for the ABC curve. Default is 3.
#' 
#' @param ShowUniform Logical. If TRUE (default), the uniform reference curve is
#'   drawn in addition to the identity and ABC curves.
#'   
#' @param Plot_title Character string. Title of the plot. Default is "ABC plot".
#' 
#' @details
#' The plot always uses a square coordinate system with both axes ranging from 0 to 1.
#' The diagonal y = 1 - x (equilibrium line) and the identity line y = x
#' are drawn as references. ABC set boundaries (A|B and B|C) are visualized with
#' stars and orthogonal boundary lines. 
#' Shows individual points if they are less then 20.
#' 
#' @return ggplot2 object
#' @noRd
cABC_plotGG <- function(CurveData, CleanData, Boundaries, Set_counts, x_vals, y_vals,
                        LineWidth = 1.25, ShowUniform = TRUE, Plot_title = 'ABC plot') {
  
  Effort <- CurveData$Curve[, 'Effort']
  Yield <- CurveData$Curve[, 'Yield']
  cleaned_data <- CleanData
  
  # Calculate uniform curve
  p_unif <- seq(0, 1, by = 0.01)
  if(!is.null(cleaned_data) && length(cleaned_data) > 0) {
    A <- min(cleaned_data, na.rm = TRUE)
    MaxX <- max(cleaned_data, na.rm = TRUE)
    if(A == MaxX) { A <- 0; MaxX <- 1 }
  } else {
    A <- 0; MaxX <- 1
  }
  B <- MaxX - A
  ABC_uniform <- (-0.5 * B * p_unif^2 + MaxX * p_unif) / (A + 0.5 * B)
  
  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 1), 
                       name = 'fraction of data') +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                       name = 'fraction of sum of largest data') +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(title = Plot_title)
  
  # Add uniform curve if requested
  if(ShowUniform) {
    uniform_df <- data.frame(x = p_unif, y = ABC_uniform)
    # Clip values to 0 and 1
    uniform_df$x <- pmin(pmax(uniform_df$x, 0), 1)
    uniform_df$y <- pmin(pmax(uniform_df$y, 0), 1)
    p <- p + ggplot2::geom_line(data = uniform_df, ggplot2::aes(x = x, y = y), 
                       color = 'green', linewidth = 1)
  }

  # Identity line
  p <- p + ggplot2::geom_line(data = data.frame(x = c(0, 1), y = c(0, 1)),
                              ggplot2::aes(x = x, y = y), color = colors()[452], linewidth = 0.1)
  
  # Equilibrium diagonal
  p <- p + ggplot2::geom_line(data = data.frame(x = c(0, 1), y = c(1, 0)),
                    ggplot2::aes(x = x, y = y), linetype = 'dashed', 
                    color = colors()[175], linewidth = 1)
  
  # ABC curve
  abc_df <- data.frame(x = Effort, y = Yield)
  p <- p + ggplot2::geom_line(data = abc_df, ggplot2::aes(x = x, y = y), 
                     color = 'blue', linewidth = LineWidth)
  
  # Show data points if less than 20
  if(length(x_vals) < 20) {
    points_df <- data.frame(x = x_vals, y = y_vals)
    p <- p + ggplot2::geom_point(data = points_df, ggplot2::aes(x = x, y = y),
                        shape = 1, size = 3, color = 'blue', stroke = 1.5)
  }
  
  # Boundary stars
  p <- p + 
    ggplot2::geom_point(data = data.frame(x = Boundaries$A[1], y = Boundaries$A[2]),
             ggplot2::aes(x = x, y = y), shape = 8, size = 3, color = 'red', stroke = 1.5) +
    ggplot2::geom_point(data = data.frame(x = Boundaries$B[1], y = Boundaries$B[2]),
             ggplot2::aes(x = x, y = y), shape = 8, size = 3, color = 'green', stroke = 1.5) +
    ggplot2::geom_point(data = data.frame(x = Boundaries$C[1], y = Boundaries$C[2]),
             ggplot2::aes(x = x, y = y), shape = 8, size = 3, color = 'blue', stroke = 1.5)
  
  # Boundary lines
  boundary_segments <- data.frame(
    x = c(0, Boundaries$A[1], 0, Boundaries$C[1]),
    y = c(Boundaries$A[2], 0, Boundaries$C[2], 0),
    xend = c(Boundaries$A[1], Boundaries$A[1], Boundaries$C[1], Boundaries$C[1]),
    yend = c(Boundaries$A[2], Boundaries$A[2], Boundaries$C[2], Boundaries$C[2])
  )
  p <- p + ggplot2::geom_segment(data = boundary_segments, 
                    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                    color = 'red', linewidth = 1)
  
  # A|B and B|C labels
  if(abs(Boundaries$A[1] - Boundaries$C[1]) > 0.1) {
    p <- p + 
      ggplot2::annotate("text", x = Boundaries$A[1], y = Boundaries$A[2], 
               label = 'A|B', color = 'red', size = 3.5, hjust = -0.2, vjust = 1.5) +
      ggplot2::annotate("text", x = Boundaries$C[1], y = Boundaries$C[2], 
               label = 'B|C', color = 'red', size = 3.5, hjust = -0.2, vjust = 1.5)
  } else {
    p <- p + 
      ggplot2::annotate("text", x = Boundaries$A[1] - 0.05, y = Boundaries$A[2] - 0.03, 
               label = 'A|B', color = 'red', size = 3.5) +
      ggplot2::annotate("text", x = Boundaries$C[1] + 0.025, y = Boundaries$C[2] - 0.025, 
               label = 'B|C', color = 'red', size = 3.5)
  }
  
  # A, B, C set labels with counts
  y_label <- Boundaries$A[2] / 4
  y_count <- y_label - 0.05
  
  # Calculate positions with boundary checks to prevent out-of-bounds labels
  pos_A <- Boundaries$A[1] / 2
  pos_B <- (Boundaries$A[1] + Boundaries$C[1]) / 2
  pos_C <- min((Boundaries$A[1] + Boundaries$C[1]) / 2 + max(abs(Boundaries$A[1] - Boundaries$C[1]), 0.1), 0.95)
  pos_C_count <- min((Boundaries$A[1] + Boundaries$C[1]) / 2 + max(abs(Boundaries$A[1] - Boundaries$C[1]), 0.1) + 0.02, 0.97)
  
  p <- p + 
    # Set A
    ggplot2::annotate("text", x = pos_A, y = y_label,
             label = 'A', color = 'red', size = 9, fontface = 'bold') +
    ggplot2::annotate("text", x = pos_A, y = y_count,
             label = paste0('n=', Set_counts$nA), color = 'black', size = 2.8) +
    # Set B
    ggplot2::annotate("text", x = pos_B, y = y_label,
             label = 'B', color = 'red', size = 7.5, fontface = 'bold') +
    ggplot2::annotate("text", x = pos_B, y = y_count,
             label = paste0('n=', Set_counts$nB), color = 'black', size = 2.8) +
    # Set C
    ggplot2::annotate("text", x = pos_C, y = y_label, 
             label = 'C', color = 'red', size = 6.5, fontface = 'bold') +
    ggplot2::annotate("text", x = pos_C_count, y = y_count, 
             label = paste0('n=', Set_counts$nC), color = 'black', size = 2.8)
  
  # Legend position
  legend_y_start <- if(((Boundaries$A[1] + Boundaries$C[1]) / 2 + 
                  max(abs(Boundaries$A[1] - Boundaries$C[1]), 0.1) + 0.02) < 0.8) {
    0.3
  } else {
    0.5
  }
  
  # Legend annotations (mimicking base R legend)
  legend_spacing <- 0.05
  label_size <- 5
  legend_x <- 0.85
  
  p <- p +
    ggplot2::annotate("text", x = legend_x, y = legend_y_start, 
             label = 'set limits', color = 'red', size = label_size, hjust = 0, fontface = 'italic') +
    ggplot2::annotate("text", x = legend_x, y = legend_y_start - legend_spacing, 
             label = 'data', color = 'blue', size = label_size, hjust = 0, fontface = 'italic') +
    ggplot2::annotate("text", x = legend_x, y = legend_y_start - 2 * legend_spacing, 
             label = 'uniform', color = 'green', size = label_size, hjust = 0, fontface = 'italic') +
    ggplot2::annotate("text", x = legend_x, y = legend_y_start - 3 * legend_spacing, 
             label = 'identity', color = colors()[452], size = label_size, hjust = 0, fontface = 'italic')
  
  # Theme
  p <- p + ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = 'grey', fill = NA, linewidth = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  return(p)
}