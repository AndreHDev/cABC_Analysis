#' Handle Special Cases Before ABC Classification
#'
#' Checks for degenerate input conditions that would make a standard ABC
#' analysis undefined or unreliable, and returns an early result or warning
#' where appropriate. This function is called by \code{cABC_analysis} before
#' the ABC curve is computed.
#'
#' The following special cases are handled:
#' \describe{
#'   \item{Single data point}{If only one positive value remains after
#'     cleaning, it is assigned to Class A and a warning is issued. The
#'     returned list has \code{Aind = 1} and all other fields empty/NULL.}
#'   \item{All values identical}{If every data point has the same value, all
#'     items are considered equally important. They are all assigned to Class A
#'     (not Class B as the warning text historically stated) and a warning is
#'     issued. The returned list has \code{Aind} set to all indices and all
#'     other fields empty/NULL.}
#'   \item{Very small dataset}{If three or fewer positive values remain after
#'     cleaning, a warning is issued that the ABC classification may be
#'     unstable, but processing continues normally and \code{NULL} is returned
#'     so that \code{cABC_analysis} proceeds with the standard algorithm.}
#' }
#'
#' @param Data Named numeric vector of positive values (already cleaned by
#'   \code{cABC_analysis}: no NAs, no non-positives, names preserved).
#'
#' @return \code{NULL} if no special case applies (normal processing should
#'   continue). Otherwise a named list with the same structure as the return
#'   value of \code{cABC_analysis}, where only \code{Aind} (and optionally
#'   \code{Bind}, \code{Cind}) are populated and all curve-related fields are
#'   \code{NULL} or empty.
cABC_handle_specials <- function(Data) {
  
  empty_result <- list(
    Aind = integer(0),
    Bind = integer(0),
    Cind = integer(0),
    ABexchanged = FALSE,
    A = NULL,
    B = NULL,
    C = NULL,
    smallestAData = NULL,
    smallestBData = NULL,
    AlimitIndInInterpolation = NULL,
    BlimitIndInInterpolation = NULL,
    p = NULL,
    ABC = NULL,
    ABLimit = NULL,
    BCLimit = NULL
  )
  
  n <- length(Data)
  # === Single data point ===
  # Set the single point to the A set
  if (n == 1) {
    warning("Only one data point remains after filtering. It has been assigned 
            to Class A per default.")
    empty_result$Aind <- 1
    names(empty_result$Aind) <- names(Data)
    return(empty_result)
  }

  # === Identical values ===
  # Assign all data points to the A set
  unique_vals <- unique(Data)
  if (length(unique_vals) == 1) {
    warning("All data values are identical, all data points are considered 
            equally important, assigning all items to Class B.")
    
    empty_result$Aind <- seq_along(Data)
    names(empty_result$Aind) <- names(Data)
    return(empty_result)
  }
  
  # === Very small dataset ===
  if (n <= 3) {
    warning("Extremly small dataset (< 3 values) after filtering. ABC classification may be unstable.
            Consider checking the plot.")
  }
  
  # No special case, continue normal ABC
  return(NULL)
}

#' Post-process ABC Classes to Resolve Boundary Duplicates
#'
#' After the initial class assignment in \code{cABC_analysis}, it is possible
#' for data points with the same value to be split across two or even all three
#' classes (A, B, C) because the geometric boundary cuts through a run of
#' identical values. This function detects such duplicates and consolidates all
#' occurrences of an ambiguous value into a single class using a deterministic
#' tie-breaking strategy.
#'
#' \strong{Tie-breaking rules:}
#' \enumerate{
#'   \item The class that contains the \emph{most} occurrences of the duplicate
#'     value wins outright.
#'   \item If all three classes are tied, the duplicate value is compared to
#'     both boundary limits. It is assigned to whichever boundary
#'     (\code{ABLimit} or \code{BCLimit}) it is closest to, then placed in the
#'     class above that boundary (i.e. closer to AB → A if \code{dup_val >=
#'     ABLimit}, else B; closer to BC → B if \code{dup_val >= BCLimit}, else
#'     C). If equidistant from both boundaries it is assigned to B.
#'   \item If exactly two classes are tied, the pair determines the rule:
#'     \itemize{
#'       \item \strong{A vs B}: compare to \code{ABLimit}; \code{>= ABLimit}
#'         → A, otherwise → B.
#'       \item \strong{B vs C}: compare to \code{BCLimit}; \code{>= BCLimit}
#'         → B, otherwise → C.
#'       \item \strong{A vs C}: always assign to A, since the value was already
#'         deemed important enough to appear in the top class.
#'     }
#' }
#'
#' A warning is issued whenever at least one duplicate boundary value is found,
#' prompting the user to inspect the data and the ABC plot.
#'
#' @param Aind Integer vector of indices currently assigned to Class A.
#' @param Bind Integer vector of indices currently assigned to Class B.
#' @param Cind Integer vector of indices currently assigned to Class C.
#' @param Data Named numeric vector of the (unsorted) input data, as cleaned by
#'   \code{cABC_analysis}.
#' @param sorted_data Numeric vector; \code{Data} sorted in decreasing order
#'   (used internally for boundary reference).
#' @param ABLimit Numeric scalar; the data value closest to the A/B boundary
#'   threshold (as computed in \code{cABC_analysis}).
#' @param BCLimit Numeric scalar; the data value closest to the B/C boundary
#'   threshold (as computed in \code{cABC_analysis}).
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{Aind}{Sorted integer vector of indices for Class A after deduplication.}
#'   \item{Bind}{Sorted integer vector of indices for Class B after deduplication.}
#'   \item{Cind}{Sorted integer vector of indices for Class C after deduplication.}
#' }
cABC_postprocess_classes <- function(Aind, Bind, Cind, Data, sorted_data, ABLimit , BCLimit) {
  # === Handle duplicates across class boundaries ===
  A_values <- Data[Aind]
  B_values <- Data[Bind]
  C_values <- Data[Cind]
  
  # Build a data.frame of value to class membership
  val_class <- data.frame(
    value = c(A_values, B_values, C_values),
    class = c(
      rep("A", length(A_values)),
      rep("B", length(B_values)),
      rep("C", length(C_values))
    )
  )
  
  # Find values that occur in more than one class
  duplicate_values <- with(
    val_class,
    tapply(class, value, function(x) length(unique(x)) > 1)
  )
  
  duplicate_values <- names(duplicate_values)[duplicate_values]
  
  if (length(duplicate_values) > 0) {
    warning(sprintf("Found %d duplicate value(s) spanning multiple classes.
    Reassigning all occurrences to the class with the most instances or based on
    distance to boundary if tied. Consider checking data and plot to confirm data 
    is suitable for ABC analysis.",length(duplicate_values)))
    
    for (dup_val in duplicate_values) {
      # Count occurrences in each class
      in_A <- Aind[Data[Aind] == dup_val]
      in_B <- Bind[Data[Bind] == dup_val]
      in_C <- Cind[Data[Cind] == dup_val]
      
      count_A <- length(in_A)
      count_B <- length(in_B)
      count_C <- length(in_C)
      
      # Determine which class has the most occurrences
      max_count <- max(count_A, count_B, count_C)
      classes_with_max <- c(
        if (count_A == max_count) "A" else NULL,
        if (count_B == max_count) "B" else NULL,
        if (count_C == max_count) "C" else NULL
      )

      # If there's a tie, use old, legacy boundary point proximity to decide
      if (length(classes_with_max) > 1) {

        # Duplicate value is in 3 classes the same amount of times
        if (length(classes_with_max) == 3) {
          dist_to_AB <- abs(dup_val - ABLimit)
          dist_to_BC <- abs(dup_val - BCLimit)
          
          if (dist_to_AB < dist_to_BC) {
            # Closer to AB boundary, so it's either A or B
            target_class <- if (dup_val >= ABLimit) "A" else "B"
          } else if (dist_to_BC < dist_to_AB) {
            # Closer to BC boundary, so it's either B or C
            target_class <- if (dup_val >= BCLimit) "B" else "C"
          } else {
            # Equidistant to both boundaries: assign to B (middle class)
            target_class <- "B"
          }
        } else if ("A" %in% classes_with_max && "B" %in% classes_with_max) {
          # AB boundary conflict: compare duplicate value to ABLimit
          target_class <- if (dup_val >= ABLimit) "A" else "B"
        } else if ("B" %in% classes_with_max && "C" %in% classes_with_max) {
          # BC boundary conflict: compare duplicate value to BCLimit
          target_class <- if (dup_val >= BCLimit) "B" else "C"
        } else if ("A" %in% classes_with_max && "C" %in% classes_with_max) {
          # AC conflict, most likely important, as one of values was deemed important
          target_class <- "A"
        }
      } else {
        target_class <- classes_with_max[1]
      }
      
      # Collect all indices with this duplicate value
      all_dup_indices <- c(in_A, in_B, in_C)
      
      # Remove from all classes
      Aind <- setdiff(Aind, all_dup_indices)
      Bind <- setdiff(Bind, all_dup_indices)
      Cind <- setdiff(Cind, all_dup_indices)
      
      # Assign all to the target class
      if (target_class == "A") {
        Aind <- c(Aind, all_dup_indices)
      } else if (target_class == "B") {
        Bind <- c(Bind, all_dup_indices)
      } else {
        Cind <- c(Cind, all_dup_indices)
      }
    }
    
    # Re-sort indices to maintain order
    Aind <- sort(Aind)
    Bind <- sort(Bind)
    Cind <- sort(Cind)
  }
  
  list(Aind = Aind, Bind = Bind, Cind = Cind)
}
