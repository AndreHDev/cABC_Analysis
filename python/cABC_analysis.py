#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: André Himmelspach

ABC Analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.spatial import distance
import warnings


def ABC_clean_data(data):
    """
    Clean the input data by retaining only positive numeric values.

    This function filters the data to include only finite, positive numeric
    values. Non-numeric inputs are warned about, and invalid values are set
    to zero. It preserves the index if the input is a pandas Series.

    Parameters
    ----------
    data : array-like or pandas Series
        The input data to be cleaned.

    Returns
    -------
    pandas Series
        A cleaned pandas Series containing only positive numeric values.
        Raises ValueError if no positive values remain.

    Notes
    -----
    - Converts lists or arrays to pandas Series.
    - Coerces non-numeric values to NaN and sets them to 0.
    - Warns if some values are filtered out.
    """
    if data is None:
        raise ValueError("No positive values remain after cleaning")

    # Handle DataFrame or 2D array
    if isinstance(data, pd.DataFrame):
        if data.shape[1] > 1:
            warnings.warn('Using only first column of data')
        data = data.iloc[:, 0]
    elif isinstance(data, np.ndarray) and data.ndim > 1:
        if data.shape[1] > 1:
            warnings.warn('Using only first column of data')
        data = data[:, 0]

    # Try to preserve index if Series
    if isinstance(data, pd.Series):
        ser = data.copy()
    else:
        ser = pd.Series(data)

    # Force numeric, suppress warnings for NAs
    ser = pd.to_numeric(ser, errors='coerce')
    ser[~np.isfinite(ser)] = 0
    ser[ser < 0] = 0

    n_original = len(ser)
    n_used = int((ser > 0).sum())
    if n_used == 0:
        raise ValueError("No positive values remain after cleaning")
    if n_used < n_original:
        warnings.warn(f"Only {n_used} of {n_original} items are larger then 0.")

    return ser


def cABC_handle_specials(Data):
    """
    Handle special cases in the data before performing full ABC analysis.

    This function checks for edge cases such as empty data, single data points,
    identical values, or very small datasets, and returns predefined results
    to avoid unstable classifications.

    Parameters
    ----------
    Data : pandas Series
        The cleaned data with preserved index.

    Returns
    -------
    dict or None
        A dictionary with empty or default results for special cases, or None
        if no special case is detected and full analysis should proceed.

    Notes
    -----
    - For single points or identical values, assigns all to class A.
    - Warns for small datasets.
    """
    empty_result = {
        'Aind': [], 'Bind': [], 'Cind': [], 'ABexchanged': False,
        'A': None, 'B': None, 'C': None, 'smallestAData': None,
        'smallestBData': None, 'AlimitIndInInterpolation': None,
        'BlimitIndInInterpolation': None, 'p': None, 'ABC': None,
        'ABLimit': None, 'BCLimit': None
    }

    n = len(Data)
    if n == 0:
        return empty_result

    # Single data point -> assign to A
    if n == 1:
        warnings.warn("Only one data point remains after filtering. It has been assigned to Class A per default.")
        empty_result['Aind'] = [Data.index[0]]
        return empty_result

    # Identical values -> assign all to A
    unique_vals = Data.unique()
    if len(unique_vals) == 1:
        warnings.warn("All data values are identical, all data points are considered equally important, assigning all items to Class A.")
        empty_result['Aind'] = list(Data.index)
        return empty_result

    # Very small dataset warning
    if n <= 3:
        warnings.warn("Extremly small dataset (< 3 values) after filtering. ABC classification may be unstable. Consider checking the plot.")

    # No special case
    return None


def cABC_postprocess_classes(Aind, Bind, Cind, Data, ABLimit, BCLimit):
    """    
    Post-process ABC Classes to Resolve Boundary Duplicates

    After the initial class assignment in cABC_analysis, it is possible
    for data points with the same value to be split across two or even all three
    classes (A, B, C) because the geometric boundary cuts through a run of
    identical values. This function detects such duplicates and consolidates all
    occurrences of an ambiguous value into a single class using a deterministic
    tie-breaking strategy.

    Tie-breaking rules:
    1. The class that contains the most occurrences of the duplicate
       value wins outright.
    2. If all three classes are tied, the duplicate value is compared to
       both boundary limits. It is assigned to whichever boundary
       (ABLimit or BCLimit) it is closest to, then placed in the
       class above that boundary (i.e. closer to AB → A if dup_val >=
       ABLimit, else B; closer to BC → B if dup_val >= BCLimit, else
       C). If equidistant from both boundaries it is assigned to B.
    3. If exactly two classes are tied, the pair determines the rule:
       - A vs B: compare to ABLimit; >= ABLimit → A, otherwise → B.
       - B vs C: compare to BCLimit; >= BCLimit → B, otherwise → C.
       - A vs C: always assign to A, since the value was already
         deemed important enough to appear in the top class.

    A warning is issued whenever at least one duplicate boundary value is found,
    prompting the user to inspect the data and the ABC plot.

    Parameters
    ----------
    Aind, Bind, Cind : list
        Lists of indices (labels) for each class.
    Data : pandas Series
        The original cleaned data.
    ABLimit, BCLimit : float
        Value limits between A/B and B/C classes.

    Returns
    -------
    dict
        Updated class indices: {'Aind': [...], 'Bind': [...], 'Cind': [...]}
    """
    val_map = {}
    for idx in Aind:
        val = Data.loc[idx]
        val_map.setdefault(val, {'A': [], 'B': [], 'C': []})['A'].append(idx)
    for idx in Bind:
        val = Data.loc[idx]
        val_map.setdefault(val, {'A': [], 'B': [], 'C': []})['B'].append(idx)
    for idx in Cind:
        val = Data.loc[idx]
        val_map.setdefault(val, {'A': [], 'B': [], 'C': []})['C'].append(idx)

    duplicate_values = [v for v, d in val_map.items() if sum([len(d['A'])>0, len(d['B'])>0, len(d['C'])>0]) > 1]

    if len(duplicate_values) > 0:
        warnings.warn(f"Found {len(duplicate_values)} duplicate value(s) spanning multiple classes. Reassigning occurrences based on majority or proximity to boundaries.")

        for dup_val in duplicate_values:
            groups = val_map[dup_val]
            in_A = groups['A']
            in_B = groups['B']
            in_C = groups['C']

            count_A = len(in_A)
            count_B = len(in_B)
            count_C = len(in_C)

            max_count = max(count_A, count_B, count_C)
            classes_with_max = []
            if count_A == max_count and max_count > 0:
                classes_with_max.append('A')
            if count_B == max_count and max_count > 0:
                classes_with_max.append('B')
            if count_C == max_count and max_count > 0:
                classes_with_max.append('C')

            if len(classes_with_max) > 1:
                if len(classes_with_max) == 3:
                    dist_to_AB = abs(dup_val - ABLimit)
                    dist_to_BC = abs(dup_val - BCLimit)
                    if dist_to_AB < dist_to_BC:
                        target_class = 'A' if dup_val >= ABLimit else 'B'
                    elif dist_to_BC < dist_to_AB:
                        target_class = 'B' if dup_val >= BCLimit else 'C'
                    else:
                        target_class = 'B'
                elif 'A' in classes_with_max and 'B' in classes_with_max:
                    target_class = 'A' if dup_val >= ABLimit else 'B'
                elif 'B' in classes_with_max and 'C' in classes_with_max:
                    target_class = 'B' if dup_val >= BCLimit else 'C'
                else:
                    target_class = 'A'
            else:
                target_class = classes_with_max[0]

            # Collect all indices with this duplicate value
            all_dup_indices = in_A + in_B + in_C

            # Remove from all classes
            Aind = [i for i in Aind if i not in all_dup_indices]
            Bind = [i for i in Bind if i not in all_dup_indices]
            Cind = [i for i in Cind if i not in all_dup_indices]

            # Assign all to target
            if target_class == 'A':
                Aind.extend(all_dup_indices)
            elif target_class == 'B':
                Bind.extend(all_dup_indices)
            else:
                Cind.extend(all_dup_indices)

        # Re-sort indices to maintain original order
        index_order = list(Data.index)
        def order_key(i):
            try:
                return index_order.index(i)
            except ValueError:
                return 0

        Aind = sorted(set(Aind), key=order_key)
        Bind = sorted(set(Bind), key=order_key)
        Cind = sorted(set(Cind), key=order_key)

    return {'Aind': Aind, 'Bind': Bind, 'Cind': Cind}

"""
Internal helper for hyman interpolation
"""
def hyman_slopes(x, y):
    # Initial slopes via finite differences
    n = len(x)
    h = np.diff(x)
    delta = np.diff(y) / h
    
    m = np.zeros(n)
    m[0] = delta[0]
    m[-1] = delta[-1]
    for i in range(1, n-1):
        m[i] = (delta[i-1] + delta[i]) / 2
    
    # Hyman monotonicity correction
    for i in range(n-1):
        if delta[i] == 0:
            m[i] = 0
            m[i+1] = 0
        else:
            alpha = m[i] / delta[i]
            beta = m[i+1] / delta[i]
            if alpha < 0: m[i] = 0
            if beta < 0: m[i+1] = 0
            if alpha**2 + beta**2 > 9:
                tau = 3 / np.sqrt(alpha**2 + beta**2)
                m[i] = tau * alpha * delta[i]
                m[i+1] = tau * beta * delta[i]
    
    return m

def cABC_curve(CleanedData):
    """
    Create the ABC curve with effort, yield, and derivatives.

    Computes the cumulative yield curve from sorted data, interpolates it,
    and calculates the derivative for identifying inflection points.

    Parameters
    ----------
    CleanedData : pandas Series or DataFrame
        The cleaned data. If Series, converted to DataFrame with 'value' column.

    Returns
    -------
    pandas DataFrame
        DataFrame with columns 'effort', 'yield', and 'dABC' (derivative).

    Notes
    -----
    - Sorts data in descending order.
    - Uses PCHIP interpolation for smoothness.
    - Ensures curve starts at (0,0) and ends at (1,1).
    """
    # If user passed a pandas Series, convert to DataFrame with column 'value'
    if isinstance(CleanedData, pd.Series):
        CleanedData = pd.DataFrame({'value': CleanedData})

    CleanedData_sorted = CleanedData.sort_values(by="value", ascending=False, inplace=False)

    Contrib = CleanedData_sorted["value"].array
    y = np.cumsum(Contrib)
    y = y/y[-1]
    y = y[~np.isinf(y)]
    x = np.arange(1, len(y)+1)/len(y)
    if np.min(y) > 0:
        y = np.insert(y, 0, 0, axis=0)
        x = np.insert(x, 0, 0, axis=0)
    if np.max(y) < 1:
        y = np.append(y, 1)
        x = np.append(x, 1)

    slopes = hyman_slopes(x, y)
    f = CubicHermiteSpline(x, y, slopes)
    Effort = np.linspace(0, 1, num=101, endpoint=True)
    Yield = f(Effort)
    if max(Yield) > 1:
        inds = np.where(Yield > 1)[0].tolist()[0]
        if inds < len(Yield):
            Yield[inds:len(Yield)] = 1
    slopes2 = hyman_slopes(Effort, Yield) 
    f = CubicHermiteSpline(Effort, Yield, slopes2)
    
    n = len(Effort)
    eval_points = np.arange(1, n + 1) / n   # to mirrors R's (1:n)/n
    dABC = f.derivative()(eval_points)

    return pd.DataFrame({"effort": Effort, "yield": Yield, "dABC": dABC})


def ABC_calc(CleanedData, ABCcurveData):
    """
    Calculate ABC class boundaries and classifications.

    Uses the ABC curve to find Pareto and break-even points, then determines
    class limits and assigns items to A, B, C based on positions.

    Parameters
    ----------
    CleanedData : pandas Series or DataFrame
        The cleaned data.
    ABCcurveData : pandas DataFrame
        The ABC curve data from ABC_curve().

    Returns
    -------
    dict
        Results including class indices, points, limits, etc.

    Notes
    -----
    - Identifies A, B, C points using distance to ideal curves.
    - Interpolates value limits.
    - Classifies based on cumulative effort positions.
    """
    # If user passed a pandas Series, convert to DataFrame with column 'value'
    if isinstance(CleanedData, pd.Series):
        CleanedData = pd.DataFrame({'value': CleanedData})

    CleanedData_sorted = CleanedData.sort_values(by="value", ascending=False, inplace=False)
    curve = ABCcurveData[["effort", "yield"]]

    # Calculate boundary points A, B, C 
    point = [[0.0, 1.0]]
    distPareto = distance.cdist(curve.to_numpy(), point, "euclidean")
    ParetoPointInd = np.where(distPareto == distPareto.min())[0][0]
    ParetoPoint = curve.iloc[[ParetoPointInd]]

    derivative = abs(ABCcurveData["dABC"] - 1)
    breakEvenInd = np.where(derivative == derivative.min())[0][0]
    breakEvenPoint = curve.iloc[[breakEvenInd]]

    if curve["effort"].iloc[breakEvenInd] < curve["effort"].iloc[ParetoPointInd]:
        ABexchanged = True
        JurenInd = breakEvenInd
        Bx = curve["effort"].iloc[ParetoPointInd]
        A, B = breakEvenPoint, ParetoPoint
    else:
        ABexchanged = False
        JurenInd = ParetoPointInd
        Bx = curve["effort"].iloc[breakEvenInd]
        A, B = ParetoPoint, breakEvenPoint

    Juren = [[Bx, 1.0]]
    distBx = distance.cdist(curve.to_numpy(), Juren, "euclidean")
    B_limit = np.where(distBx == distBx.min())[0][0]
    C = curve.iloc[[B_limit]]

    # Calculate value limits (for both methods)
    f = PchipInterpolator(np.linspace(1, 100, num=len(CleanedData_sorted["value"]), endpoint=True),
                          CleanedData_sorted["value"])
    interpolatedInverseEcdf = f(np.linspace(1, 100, num=1000, endpoint=True))
    ABlimit = interpolatedInverseEcdf[round(A.values.tolist()[0][0] * 1000)]
    BClimit = interpolatedInverseEcdf[round(C.values.tolist()[0][0] * 1000)]

    # Curve-based classification (positions along sorted data)
    x_vals = np.arange(1, len(CleanedData_sorted) + 1) / len(CleanedData_sorted)
    px_A = A.values.tolist()[0][0]
    px_C = C.values.tolist()[0][0]

    idx_left_A = x_vals <= px_A
    idx_left_C = x_vals <= px_C

    sorted_labels = list(CleanedData_sorted.index)
    A_positions = np.where(idx_left_A)[0]
    B_positions = np.where(idx_left_C & ~idx_left_A)[0]
    C_positions = np.where(~idx_left_C)[0]

    Aind = [sorted_labels[i] for i in A_positions]
    Bind = [sorted_labels[i] for i in B_positions]
    Cind = [sorted_labels[i] for i in C_positions]

    smallestAData = curve["yield"].iloc[JurenInd]
    smallestBData = curve["yield"].iloc[B_limit]

    return {
        "Aind": Aind, "Bind": Bind, "Cind": Cind, "ABexchanged": ABexchanged,
        "A": A, "B": B, "C": C, "smallestAData": smallestAData,
        "smallestBData": smallestBData, "AlimitIndInInterpolation": JurenInd,
        "BlimitIndInInterpolation": B_limit, "p": curve["effort"], "ABC": curve["yield"],
        "ABLimit": ABlimit, "BCLimit": BClimit
    }


def ABC_plot(cABCresults, CleanedData, ax=None):
    """
    Plot the ABC analysis results.

    Creates a visualization of the ABC curve, boundaries, and class annotations.

    Parameters
    ----------
    cABCresults : dict
        Results from ABC_calc().
    CleanedData : pandas Series or DataFrame
        The cleaned data.
    ax : matplotlib Axes, optional
        Axes to plot on; if None, uses current axes.

    Returns
    -------
    matplotlib Axes
        The axes object with the plot.

    Notes
    -----
    - Plots data points if dataset is small (<20).
    - Includes reference lines for identity and uniform distribution.
    - Annotates class sizes.
    """
    # If user passed a pandas Series, convert to DataFrame with column 'value'
    if isinstance(CleanedData, pd.Series):
        CleanedData = pd.DataFrame({'value': CleanedData})

    CleanedData_sorted = CleanedData.sort_values(by="value", ascending=False, inplace=False)

    # Data points (blue circles)
    Contrib = CleanedData_sorted["value"].array
    y = np.cumsum(Contrib)
    y = y/y[-1]
    y = y[~np.isinf(y)]
    x = np.arange(1, len(y)+1)/len(y)

    # Reference curves
    pIdent = np.linspace(0, 1, 100)
    A_min = CleanedData_sorted["value"].min()
    MaxX = CleanedData_sorted["value"].max()
    if A_min == MaxX:
        A_min = 0
        MaxX = 1
    Bmax = MaxX - A_min
    ABCuniform = (-0.5 * Bmax * pIdent**2 + MaxX * pIdent)/(A_min + 0.5 * Bmax)

    ax = ax or plt.gca()

    # Data points: only draw individual dots when dataset is small (<20)
    if len(x) < 20:
        ax.scatter(x, y, facecolors='none', edgecolors='blue', s=50, linewidth=1.5)

    # ABC curve
    ax.plot(cABCresults["p"], cABCresults["ABC"], color="dodgerblue", linewidth=3, label="ABC curve")

    # Boundary lines
    A_x, A_y = cABCresults["A"].values.tolist()[0]
    B_x, B_y = cABCresults["B"].values.tolist()[0]
    C_x, C_y = cABCresults["C"].values.tolist()[0]

    ax.plot([A_x, A_x], [0, A_y], color="salmon", linewidth=2)
    ax.plot([0, A_x], [A_y, A_y], color="salmon", linewidth=2)
    ax.plot([C_x, C_x], [0, C_y], color="salmon", linewidth=2)
    ax.plot([0, C_x], [C_y, C_y], color="salmon", linewidth=2)

    # ABC boundary stars
    ax.scatter(A_x, A_y, marker='*', s=120, color='red', linewidth=1.5, zorder=10)
    ax.scatter(B_x, B_y, marker='*', s=120, color='green', linewidth=1.5, zorder=10)
    ax.scatter(C_x, C_y, marker='*', s=120, color='blue', linewidth=1.5, zorder=10)

    # Reference lines
    ax.plot(pIdent, pIdent, color="magenta", linestyle="dashed", linewidth=1, label="Identity")
    ax.plot(pIdent, ABCuniform, color="green", linestyle="dotted", linewidth=1, label="Uniform")

    # Set annotations
    ax.text(0.5 * A_x, .1, f"Set A:\nn = {len(cABCresults['Aind'])}",
            ha='center', va='bottom', size=10, color='blue', weight='bold')
    ax.text(0.5 * (A_x + C_x), .1, f"Set B:\nn = {len(cABCresults['Bind'])}",
            ha='center', va='bottom', size=9, weight='semibold')
    ax.text(0.5 * (1 + C_x), .1, f"Set C:\nn = {len(cABCresults['Cind'])}",
            ha='center', va='bottom', size=9, weight='semibold')


    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Cumulative Effort")
    ax.set_ylabel("Cumulative Yield")
    ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98))

    return ax


def cABC_analysis(data, plot_it=False, ax=None):
    """
    ABC Classification

    Divides a numeric dataset into three classes (A, B, and C) using
    ABC analysis. The classification is based on geometric properties
    of the ABC curve and identifies regions of high, balanced, and
    low efficiency.

    Class interpretation:
    - A: Low effort, high yield (Pareto items)
    - B: Balanced effort and yield
    - C: High effort, low yield (submarginal items)

    Parameters
    ----------
    data : array-like or pandas Series/DataFrame
        Positive numeric vector which is not uniformly distributed.
        If DataFrame, the first column will be used.
    plot_it : bool, default False
        If True, an ABC plot is generated using matplotlib.
    ax : matplotlib Axes, optional
        Axes for plotting; if None and plot_it=True, creates new figure.

    Returns
    -------
    dict or None
        A dictionary containing:
        - Aind, Bind, Cind: Lists of indices (into the original data) for items
          assigned to classes A, B, and C respectively. In special-case returns
          (single point or all-identical), only Aind is populated; Bind and Cind
          are empty lists.
        - ABexchanged: bool; True if the Pareto point and Break-even point were
          swapped to maintain coordinate logic (i.e. the Break-even point was to
          the left of the Pareto point on the curve).
        - A, B, C: Coordinates [x, y] for the Pareto point (A), the Break-even
          point (B), and the Submarginal point (C). None in special-case returns.
        - smallestAData: Cumulative yield at the boundary of Class A. None in
          special-case returns.
        - smallestBData: Cumulative yield at the boundary of Class B. None in
          special-case returns.
        - AlimitIndInInterpolation: Index of the A boundary in the interpolated
          [p, ABC] curve. None in special-case returns.
        - BlimitIndInInterpolation: Index of the C boundary in the interpolated
          [p, ABC] curve. None in special-case returns.
        - p: Numeric array of effort values (x-axis) of the interpolation curve.
          None in special-case returns.
        - ABC: Numeric array of yield values (y-axis) of the interpolation curve.
          None in special-case returns.
        - ABLimit: Data value closest to the threshold separating Class A from
          Class B. None in special-case returns.
        - BCLimit: Data value closest to the threshold separating Class B from
          Class C. None in special-case returns.

    Notes
    -----
    Data cleaning: Before classification, non-numeric values and NAs are coerced
    to 0, negative values are set to 0. A warning is issued when items are converted.
    If a DataFrame is supplied, only the first column is used.

    Degenerate inputs (single point, all-identical values, very small datasets)
    are caught before curve fitting. Boundary duplicate values that span two
    classes after classification are resolved. In both cases a warning is issued
    when a special case is triggered.
    """
    CleanedData = ABC_clean_data(data)

    if len(CleanedData) == 0:
        print("No valid data after cleaning")
        return None
    # Handle special cases before computing curve
    special = cABC_handle_specials(CleanedData)
    if special is not None:
        return special

    ABCcurveData = cABC_curve(CleanedData)
    cABCresults = ABC_calc(CleanedData, ABCcurveData)

    # Postprocess classes to resolve duplicate values spanning classes
    ABLimit = cABCresults.get('ABLimit', None)
    BCLimit = cABCresults.get('BCLimit', None)

    # Ensure we have label-list A/B/C to pass to postprocessor
    Aind = cABCresults['Aind'] if isinstance(cABCresults['Aind'], list) else list(cABCresults['Aind'])
    Bind = cABCresults['Bind'] if isinstance(cABCresults['Bind'], list) else list(cABCresults['Bind'])
    Cind = cABCresults['Cind'] if isinstance(cABCresults['Cind'], list) else list(cABCresults['Cind'])

    processed = cABC_postprocess_classes(Aind, Bind, Cind, CleanedData, ABLimit, BCLimit)
    cABCresults['Aind'] = processed['Aind']
    cABCresults['Bind'] = processed['Bind']
    cABCresults['Cind'] = processed['Cind']

    print(f"Result: (nA={len(cABCresults['Aind'])}, nB={len(cABCresults['Bind'])}, nC={len(cABCresults['Cind'])})")

    if plot_it:
        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 10))
        ABC_plot(cABCresults, CleanedData, ax)
        plt.tight_layout()

    return cABCresults

