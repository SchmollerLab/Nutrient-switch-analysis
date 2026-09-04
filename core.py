import itertools

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu, skew, shapiro, normaltest, chi2_contingency
from IPython.display import display



def _run_pairwise_test(xa, xb, test):
    """Run the selected pairwise test. Returns p-value."""
    test_name = str(test).strip().lower()
    if test_name in ("welch", "welch_t", "welch_ttest"):
        _, p = ttest_ind(xa, xb, equal_var=False)
    elif test_name in ("mwu", "mannwhitney", "mannwhitneyu"):
        _, p = mannwhitneyu(xa, xb, alternative="two-sided")
    elif test_name in ("bs", "bootstrap"):
        p = sig_test_for_bs(xa, xb)
    else:
        raise ValueError(
            f"Unsupported test '{test}'. Use 'mwu' (Mann-Whitney U), 'welch' (Welch's t-test), "
            f"or 'bs' (bootstrap CI)."
        )
    return p


def _normalize_suppress_labels(suppress_labels):
    """Normalize suppression labels into a set for fast membership checks."""
    if suppress_labels is None:
        return set()
    return {str(lbl) for lbl in suppress_labels}


strain_labels = {"MMY": "WT", "WHI5DEL": "whi5Δ", "BCK2DEL": "bck2Δ", "WHI5_BCK2_DD": "DD"}

def run_stats(plotting_df, value_cols, comparisons=None, group_col="strain", 
              group_labels=strain_labels, title=None,
              ):
    """
    Print distribution summaries and pairwise statistical tests.

    Parameters
    ----------
    plotting_df : pd.DataFrame
    value_cols : str or list of str
        Column name(s) to measure. E.g. "cell_vol_fl" or ["phase_length_mins", "phase_volume_at_end"].
    comparisons : list of tuples (a, b)
        Each tuple defines one pairwise comparison between two group values.
        Example: [("MMY", "WHI5DEL"), ("MMY", "BCK2DEL")]
    group_col : str, optional
        Column to group by. Defaults to "strain".
    group_labels : dict, optional
        Mapping from group value to display label. Falls back to the raw value.
    title : str, optional
        Title to display for the statistics output.
    tests : list of str, optional
        Compatibility for comparisons
    """
    if group_labels is None:
        group_labels = {}

    if isinstance(value_cols, str):
        value_cols = [value_cols]
        
    if isinstance(group_col, str):
        group_cols = [group_col]
    else:
        group_cols = group_col
        
    if comparisons is None:
        groups = []
        for group_col in group_cols:
            groups.append(plotting_df.loc[:, group_col].dropna().unique())
        # we now have all unique entries to the groups, we should now create all possible combinations of the group entries
        possible_group_combinations = list(itertools.product(*groups))
        _comparisons = list(itertools.combinations(possible_group_combinations, 2))
    else:
        _comparisons = comparisons
        
    for value_col in value_cols:
                
        all_groups = set([g for pair in _comparisons for g in pair])

        print("\n" + "=" * 80)
        if title:
            print(title)
        print(f"Group: '{group_cols}'  |  Value: '{value_col}'")

        # Per-group distribution summary
        summary_rows = []
        for group in all_groups:
            values = group if isinstance(group, tuple) else (group,)
            mask = plotting_df[group_cols].eq(values).all(axis=1)
            x = plotting_df.loc[mask, value_col].dropna().values
            
            if len(x) < 3:
                print(f"    Not enough data points ({len(x)}) for group '{group}', skipping.")
                continue
            q1, q3 = np.percentile(x, [25, 75])
            iqr = q3 - q1
            n_outliers = int(np.sum((x < q1 - 1.5 * iqr) | (x > q3 + 1.5 * iqr)))
            # Shapiro-Wilk for n<=5000, D'Agostino-Pearson for larger samples
            if len(x) <= 5000:
                norm_stat, norm_p = shapiro(x)
                norm_test = "Shapiro-Wilk"
            else:
                norm_stat, norm_p = normaltest(x)
                norm_test = "D'Agostino-Pearson"
            summary_rows.append({
                group_col: group_labels.get(group, group),
                "n": len(x),
                "mean": round(np.mean(x), 1),
                "median": round(np.median(x), 1),
                "skewness": round(skew(x), 2),
                "n_outliers (IQR)": n_outliers,
                "normality_test": norm_test,
                "normality_p": round(norm_p, 10),
                "normal (p>0.05)": norm_p > 0.05,
            })
        display(pd.DataFrame(summary_rows))

        # Pairwise tests for each requested comparison
        test_rows = []
        for a, b in _comparisons:
            values_a = a if isinstance(a, tuple) else (a,)
            values_b = b if isinstance(b, tuple) else (b,)
            mask_a = plotting_df[group_cols].eq(values_a).all(axis=1)
            mask_b = plotting_df[group_cols].eq(values_b).all(axis=1)
            xa = plotting_df.loc[mask_a, value_col].dropna().values
            xb = plotting_df.loc[mask_b, value_col].dropna().values

            if len(xa) < 3 or len(xb) < 3:
                print(f"    Not enough data points for comparison '{a}' vs '{b}', skipping.")
                continue

            tw_stat, tw_p = ttest_ind(xa, xb, equal_var=False)
            t_stat, t_p = ttest_ind(xa, xb, equal_var=True)
            u_stat, u_p = mannwhitneyu(xa, xb, alternative="two-sided")

            label_a = group_labels.get(a, a)
            label_b = group_labels.get(b, b)
            test_rows.append({
                "comparison": f"{label_a} vs {label_b}",
                "n_1": len(xa),
                "n_2": len(xb),
                "Welch_t_p": round(tw_p, 4),
                "t_p": round(t_p, 4),
                "MannWhitneyU_p": round(u_p, 4),
                "tw_sig (p<0.05)": tw_p < 0.05,
                "t_sig (p<0.05)": t_p < 0.05,
                "MWU_sig (p<0.05)": u_p < 0.05,
                # check if the tests make different conclusions about significance
                "Any difference": (tw_p < 0.05) != (t_p < 0.05) or (tw_p < 0.05) != (u_p < 0.05) or (t_p < 0.05) != (u_p < 0.05),
            })

        display(pd.DataFrame(test_rows))


def _whisker_top(arr):
    """Upper whisker value (Q3 + 1.5*IQR, clipped to data max)."""
    if len(arr) == 0:
        return np.nan
    q1, q3 = np.percentile(arr, [25, 75])
    upper = q3 + 1.5 * (q3 - q1)
    above = arr[arr <= upper]
    return above.max() if len(above) > 0 else q3


def _plot_element_top(ax, x_center, tol):
    """
    Return the highest rendered y value (bar top or error bar cap) within
    *tol* x-units of *x_center*.  Falls back to None if nothing is found.
    This lets bracket placement use what is actually drawn rather than the
    raw data distribution, which matters for bar charts where the bars show
    means ± SE but the underlying data spreads much wider.
    """
    best = None
    # Bar patches
    for patch in ax.patches:
        if hasattr(patch, "get_x") and hasattr(patch, "get_width"):
            bx = patch.get_x() + patch.get_width() / 2
            if abs(bx - x_center) < tol:
                h = patch.get_height() if hasattr(patch, "get_height") else None
                if h is None:
                    continue
                best = h if best is None else max(best, h)
        else:
            # Path-based artists (e.g. seaborn boxplot PathPatch) don't expose get_x/get_width.
            # Convert vertices to data coordinates and infer center and top from extents.
            try:
                verts = patch.get_path().vertices
                verts_disp = patch.get_transform().transform(verts)
                verts_data = ax.transData.inverted().transform(verts_disp)
                x_vals = verts_data[:, 0]
                y_vals = verts_data[:, 1]
                bx = (np.nanmin(x_vals) + np.nanmax(x_vals)) / 2
                if abs(bx - x_center) < tol:
                    h = np.nanmax(y_vals)
                    best = h if best is None else max(best, h)
            except Exception:
                # Skip artists that do not expose a usable path/transform.
                continue
    # Error bar lines (vertical stems and caps are all Line2D objects)
    for line in ax.lines:
        xdata = line.get_xdata()
        ydata = line.get_ydata()
        if len(xdata) == 0:
            continue
        if np.any(np.abs(np.asarray(xdata, dtype=float) - x_center) < tol):
            top = float(np.nanmax(ydata))
            best = top if best is None else max(best, top)
    return best

def _p_to_stars(p):
    if p < 0.001:  return "***"
    if p < 0.01:   return "**"
    if p < 0.05:   return "*"
    return "ns"


def annotate_significance(ax, df, x_col, y_col, hue_col, x_order, hue_order,
                          min_n=5, star_fontsize=18, other_fontsize=12,
                          line_color="black", line_width=1,
                          y_gap=0.02, text_gap=0.01, star_gap=-0.02,
                          show_n=True, n_fontsize=10, n_offset=0.005,
                          test="mwu", suppress_labels=None,
                          reserve_suppressed=True):
    """
    Draw Mann-Whitney U significance annotations on a grouped boxplot axis.

    A horizontal line is drawn between each pair of hue groups for every
    x-category. No vertical tick marks are added. Stars are rendered at a
    larger font size than 'ns' / 'nd' labels.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis that already contains the boxplot.
    df : pd.DataFrame
        The data used for plotting.  Must contain x_col, y_col, hue_col.
    x_col : str
        Column mapped to the x-axis (bin labels, etc.).
    y_col : str
        Column mapped to the y-axis (measured values).
    hue_col : str
        Column mapped to hue (group labels after any renaming).
    x_order : list
        Ordered list of x-category values.
    hue_order : list of length 2
        The two hue groups to compare (left, right).
    min_n : int, optional
        Minimum number of data points per group to run the test. Default 5.
    star_fontsize : int, optional
        Font size for significant star labels. Default 14.
    other_fontsize : int, optional
        Font size for 'ns' and 'nd' labels. Default 9.
    line_color : str, optional
        Colour of the horizontal annotation line. Default 'black'.
    line_width : float, optional
        Line width of the annotation line. Default 1.
    y_gap : float, optional
        Gap above the highest whisker as a fraction of the current y-range.
        Default 0.02.
    text_gap : float, optional
        Additional vertical offset of the text above the line, as a fraction
        of the y-range.  Positive moves text upward.  Default 0.002.
    star_gap : float, optional
        Additional vertical offset of star labels above the line, as a fraction
        of the y-range.  Positive moves text upward.  Default 0.01.
    show_n : bool, optional
        If True, print n for each group below each box. Default True.
    n_fontsize : int, optional
        Font size for n labels. Default 8.
    suppress_labels : iterable of str, optional
        Labels to skip drawing (e.g. {"***"}). Default None.
    reserve_suppressed : bool, optional
        If True, suppressed labels still reserve vertical layout space.
        Default True.
    test : str, optional
        Statistical test to use for pairwise comparisons. Default "mwu".
    """


    n_hue = len(hue_order)
    dodge = 0.8 / n_hue   # seaborn default total box width = 0.8
    suppress_labels = _normalize_suppress_labels(suppress_labels)

    y_lo, y_hi = ax.get_ylim()
    y_range = y_hi - y_lo
    # reserve head room for annotations
    ax.set_ylim(y_lo, y_hi + 0.18 * y_range)
    y_lo, y_hi = ax.get_ylim()
    y_range = y_hi - y_lo

    for j, b in enumerate(x_order):
        xa = df.loc[(df[x_col] == b) & (df[hue_col] == hue_order[0]), y_col].dropna().values
        xb = df.loc[(df[x_col] == b) & (df[hue_col] == hue_order[1]), y_col].dropna().values

        if len(xa) < min_n or len(xb) < min_n:
            label = "nd"
        else:
            p = _run_pairwise_test(xa, xb, test)
            label = _p_to_stars(p)

        suppressed = label in suppress_labels

        # Use rendered bar/error-bar tops; fall back to whisker if not found
        plot_top_a = _plot_element_top(ax, j - dodge / 2, tol=dodge * 0.45)
        plot_top_b = _plot_element_top(ax, j + dodge / 2, tol=dodge * 0.45)
        if plot_top_a is not None and plot_top_b is not None:
            y_top = max(plot_top_a, plot_top_b)
        else:
            y_top = max(
                _whisker_top(xa) if len(xa) > 0 else y_lo,
                _whisker_top(xb) if len(xb) > 0 else y_lo,
            )
        y_line = y_top + y_gap * y_range
        x_left  = j - dodge / 2
        x_right = j + dodge / 2

        if suppressed and not reserve_suppressed:
            continue

        if not suppressed:
            ax.hlines(y_line, x_left, x_right,
                  colors=line_color, linewidth=line_width, zorder=5)

            fs = star_fontsize if label not in ("ns", "nd") else other_fontsize
            y_text = y_line + (star_gap if label not in ("ns", "nd") else text_gap) * y_range
            ax.text(j, y_text, label,
                ha="center", va="bottom", fontsize=fs, zorder=5)

        if show_n:
            y_lo_cur = ax.get_ylim()[0]

            ax.text(j - dodge / 2, y_lo_cur - n_offset * y_range, f"n={len(xa)}",
                    ha="center", va="top", fontsize=n_fontsize,
                    color="dimgray", zorder=5)
            ax.text(j + dodge / 2, y_lo_cur - n_offset * y_range, f"n={len(xb)}",
                    ha="center", va="top", fontsize=n_fontsize,
                    color="dimgray", zorder=5)


def annotate_significance_vs_ref(ax, df, x_col, y_col, x_order, ref_group,
                                  min_n=5, star_fontsize=18, other_fontsize=12,
                                  line_color="black", line_width=1,
                                  y_gap=0.05, bracket_step=0.07, test_buffer=0.02,
                                  text_gap=0.01, star_gap=-0.02,
                                  show_n=True, n_fontsize=10, n_offset=0.005,
                                  test="mwu", suppress_labels=None,
                                  reserve_suppressed=True,
                                  exclude_groups=None):
    """
    Draw significance bracket annotations on a plain (no-hue) boxplot,
    comparing a reference group to every other group.

    Brackets are positioned above each pair's whisker tops. Collision
    detection is done in a fast single pass: if a newly-added bracket
    overlaps an existing bracket in x-span and is within ``test_buffer``
    (fraction of y-range), it is pushed up by ``bracket_step``.
    No vertical tick marks are drawn.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis that already contains the boxplot.
    df : pd.DataFrame
        The data used for plotting.  Must contain x_col and y_col.
    x_col : str
        Column mapped to the x-axis (strain labels, etc.).
    y_col : str
        Column mapped to the y-axis (measured values).
    x_order : list
        Ordered list of x-category values (must include ref_group).
    ref_group : str
        The x-category to compare all others against (e.g. WT).
    min_n : int, optional
        Minimum n per group to run the test. Default 5.
    star_fontsize : int, optional
        Font size for star labels. Default 18.
    other_fontsize : int, optional
        Font size for 'ns' / 'nd' labels. Default 12.
    line_color : str, optional
        Colour of bracket lines. Default 'black'.
    line_width : float, optional
        Line width. Default 1.
    y_gap : float, optional
        Gap above the highest whisker as a fraction of the original y-range
        for the first bracket. Default 0.05.
    bracket_step : float, optional
        Vertical offset applied when a collision is detected, as a fraction
        of the original y-range. Default 0.07.
    test_buffer : float, optional
        Collision margin between nearby bracket y-levels (fraction of the
        original y-range). Default 0.02.
    text_gap : float, optional
        Vertical offset of 'ns' text above bracket line. Default 0.01.
    star_gap : float, optional
        Vertical offset of star text above bracket line. Default -0.02.
    show_n : bool, optional
        Print n below each box. Default True.
    n_fontsize : int, optional
        Font size for n labels. Default 10.
    n_offset : float, optional
        How far below the y-axis minimum n labels are placed. Default 0.005.
    test : str, optional
        Statistical test to use for pairwise comparisons. Supported values:
        - "mwu" (Mann-Whitney U test)
        - "welch" (Welch's t-test, unequal variances)
        - "bs" (bootstrap CI, uses sig_test_for_bs)
        Default "mwu".
    suppress_labels : iterable of str, optional
        Labels to skip drawing (e.g. {"***"}). Default None.
    reserve_suppressed : bool, optional
        If True, suppressed labels still reserve vertical layout space.
        Default True.
    """

    def _p_to_stars(p):
        if p < 0.001:  return "***"
        if p < 0.01:   return "**"
        if p < 0.05:   return "*"
        return "ns"

    ref_idx = x_order.index(ref_group)
    xa = df.loc[df[x_col] == ref_group, y_col].dropna().values
    suppress_labels = _normalize_suppress_labels(suppress_labels)

    y_lo, y_hi = ax.get_ylim()
    y_range = y_hi - y_lo
    min_spacing = bracket_step * y_range
    collision_margin = test_buffer * y_range

    _exclude = set(exclude_groups) if exclude_groups is not None else set()

    # --- Build bracket descriptors -------------------------------------------
    brackets = []
    for j, group in enumerate(x_order):
        if group == ref_group or group in _exclude:
            continue
        xb = df.loc[df[x_col] == group, y_col].dropna().values

        if len(xa) < min_n or len(xb) < min_n:
            label = "nd"
        else:
            p = _run_pairwise_test(xa, xb, test)
            label = _p_to_stars(p)

        suppressed = label in suppress_labels

        x_ref_c = float(ref_idx)
        x_grp_c = float(j)
        plot_top_ref = _plot_element_top(ax, x_ref_c, tol=0.4)
        plot_top_grp = _plot_element_top(ax, x_grp_c, tol=0.4)
        if plot_top_ref is not None and plot_top_grp is not None:
            y_top = max(plot_top_ref, plot_top_grp)
        else:
            y_top = max(
                _whisker_top(xa) if len(xa) > 0 else y_lo,
                _whisker_top(xb) if len(xb) > 0 else y_lo,
            )
        if suppressed and not reserve_suppressed:
            continue

        brackets.append({
            'j': j,
            'x_left':  min(ref_idx, j),
            'x_right': max(ref_idx, j),
            'label':   label,
            'suppressed': suppressed,
            'y_natural': y_top + y_gap * y_range,
        })

    # Respect plotting order: left-most bracket first, then moving right.
    brackets.sort(key=lambda b: (b['x_left'], b['x_right']))

    # --- Single-pass collision detection --------------------------------------
    placed = []   # list of (x_left, x_right, y_final)
    prev_y = None
    for b in brackets:
        y_final = b['y_natural']

        # If the new test is too close to any overlapping prior test,
        # place it above those by one bracket step.
        floors = []
        for (px_l, px_r, py) in placed:
            x_overlap = (b['x_left'] <= px_r) and (b['x_right'] >= px_l)
            if x_overlap and abs(y_final - py) <= collision_margin:
                floors.append(py + min_spacing)
        if floors:
            y_final = max(y_final, max(floors))

        # Enforce strict plotting-order stacking: every next bracket is above
        # the previous bracket, independent of overlap checks.
        if prev_y is not None:
            y_final = max(y_final, prev_y + min_spacing)

        b['y_final'] = y_final
        placed.append((b['x_left'], b['x_right'], y_final))
        prev_y = y_final

    # Expand y-axis to fit all brackets + text headroom
    if brackets:
        max_y = max(b['y_final'] for b in brackets)
        needed = max_y + (0.12 + abs(text_gap)) * y_range
        if needed > y_hi:
            ax.set_ylim(y_lo, needed)
            y_lo, y_hi = ax.get_ylim()
            y_range = y_hi - y_lo

    # --- Draw ----------------------------------------------------------------
    for b in brackets:
        y_line = b['y_final']
        label  = b['label']
        x_mid  = (ref_idx + b['j']) / 2

        if b.get('suppressed', False):
            continue

        ax.hlines(y_line, b['x_left'], b['x_right'],
                  colors=line_color, linewidth=line_width, zorder=5)

        fs = star_fontsize if label not in ("ns", "nd") else other_fontsize
        y_text = y_line + (star_gap if label not in ("ns", "nd") else text_gap) * y_range
        ax.text(x_mid, y_text, label,
                ha="center", va="bottom", fontsize=fs, zorder=5)

    if show_n:
        y_lo_cur = ax.get_ylim()[0]
        for j, group in enumerate(x_order):
            x_vals = df.loc[df[x_col] == group, y_col].dropna().values
            ax.text(j, y_lo_cur - n_offset * y_range, f"n={len(x_vals)}",
                    ha="center", va="top", fontsize=n_fontsize,
                    color="dimgray", zorder=5)


def annotate_significance_hue_vs_ref(ax, df, x_col, y_col, hue_col, x_order, hue_order, ref_hue,
                                      min_n=5, star_fontsize=18, other_fontsize=12,
                                      line_color="black", line_width=1,
                                      y_gap=0.05, bracket_step=0.07, test_buffer=0.02,
                                      text_gap=0.01, star_gap=-0.02,
                                      show_n=True, n_fontsize=10, n_offset=0.005,
                                      box_width=0.8, test="mwu", suppress_labels=None,
                                      reserve_suppressed=True):
    """
    Draw MWU significance bracket annotations on a hue-grouped boxplot,
    comparing a reference hue group to every other hue group, for each
    x-category.

    Brackets are drawn between the reference hue box and each other hue box.
    Collision detection is done per x-category: brackets are sorted by distance
    from the reference (nearest first) and stacked upward when they overlap.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis that already contains the boxplot.
    df : pd.DataFrame
        The data used for plotting. Must contain x_col, y_col, hue_col.
    x_col : str
        Column mapped to the x-axis.
    y_col : str
        Column mapped to the y-axis (measured values).
    hue_col : str
        Column mapped to hue (strain labels, etc.).
    x_order : list
        Ordered list of x-category values.
    hue_order : list
        Ordered list of hue values, matching seaborn hue_order.
    ref_hue : str
        The hue group to compare all others against (e.g. 'MMY').
    min_n : int, optional
        Minimum n per group to run the test. Default 5.
    star_fontsize : int, optional
        Font size for star labels. Default 18.
    other_fontsize : int, optional
        Font size for 'ns' / 'nd' labels. Default 12.
    line_color : str, optional
        Colour of bracket lines. Default 'black'.
    line_width : float, optional
        Line width. Default 1.
    y_gap : float, optional
        Gap above the highest whisker as a fraction of the y-range for the
        lowest bracket. Default 0.05.
    bracket_step : float, optional
        Vertical step when a collision is detected, as a fraction of the
        y-range. Default 0.07.
    test_buffer : float, optional
        Collision margin (fraction of y-range). Default 0.02.
    text_gap : float, optional
        Vertical offset of 'ns' text above the bracket line. Default 0.01.
    star_gap : float, optional
        Vertical offset of star text above the bracket line. Default -0.02.
    show_n : bool, optional
        Print n below each box. Default True.
    n_fontsize : int, optional
        Font size for n labels. Default 10.
    n_offset : float, optional
        Distance below the y-axis minimum for n labels. Default 0.005.
    box_width : float, optional
        Total width occupied by all hue boxes at one x-tick (seaborn default
        is 0.8). Default 0.8.
    test : str, optional
        Statistical test to use for pairwise comparisons. Supported values:
        - "mwu" (Mann-Whitney U test)
        - "welch" (Welch's t-test, unequal variances)
        - "bs" (bootstrap CI, uses sig_test_for_bs)
        Default "mwu".
    suppress_labels : iterable of str, optional
        Labels to skip drawing (e.g. {"***"}). Default None.
    reserve_suppressed : bool, optional
        If True, suppressed labels still reserve vertical layout space.
        Default True.
    """

    def _p_to_stars(p):
        if p < 0.001:  return "***"
        if p < 0.01:   return "**"
        if p < 0.05:   return "*"
        return "ns"

    n_hue = len(hue_order)
    dodge = box_width / n_hue
    ref_hue_idx = hue_order.index(ref_hue)
    suppress_labels = _normalize_suppress_labels(suppress_labels)

    y_lo, y_hi = ax.get_ylim()
    y_range = y_hi - y_lo
    min_spacing = bracket_step * y_range
    collision_margin = test_buffer * y_range

    for xi, x_val in enumerate(x_order):
        xa = df.loc[(df[x_col] == x_val) & (df[hue_col] == ref_hue), y_col].dropna().values
        x_ref_center = xi + (ref_hue_idx - (n_hue - 1) / 2) * dodge

        # --- Build bracket descriptors -------------------------------------------
        brackets = []
        for hi, hue_val in enumerate(hue_order):
            if hue_val == ref_hue:
                continue
            xb = df.loc[(df[x_col] == x_val) & (df[hue_col] == hue_val), y_col].dropna().values

            if len(xa) < min_n or len(xb) < min_n:
                label = "nd"
            else:
                p = _run_pairwise_test(xa, xb, test)
                label = _p_to_stars(p)

            suppressed = label in suppress_labels

            x_hue_center = xi + (hi - (n_hue - 1) / 2) * dodge

            # Prefer reading the rendered bar/error-bar top so bracket
            # placement works for bar charts (mean ± SE) as well as box plots.
            plot_top_ref  = _plot_element_top(ax, x_ref_center,  tol=dodge * 0.45)
            plot_top_hue  = _plot_element_top(ax, x_hue_center,  tol=dodge * 0.45)
            if plot_top_ref is not None and plot_top_hue is not None:
                y_top = max(plot_top_ref, plot_top_hue)
            else:
                y_top = max(
                    _whisker_top(xa) if len(xa) > 0 else y_lo,
                    _whisker_top(xb) if len(xb) > 0 else y_lo,
                )
            if suppressed and not reserve_suppressed:
                continue

            brackets.append({
                'hi': hi,
                'x_left':  min(x_ref_center, x_hue_center),
                'x_right': max(x_ref_center, x_hue_center),
                'x_mid':   (x_ref_center + x_hue_center) / 2,
                'label':   label,
                'suppressed': suppressed,
                'y_natural': y_top + y_gap * y_range,
            })

        # Sort: narrowest bracket (closest hue to ref) first so it sits lowest.
        brackets.sort(key=lambda b: b['x_right'] - b['x_left'])

        # --- Single-pass collision detection -------------------------------------
        placed = []
        prev_y = None
        for b in brackets:
            y_final = b['y_natural']

            floors = []
            for (px_l, px_r, py) in placed:
                x_overlap = (b['x_left'] <= px_r) and (b['x_right'] >= px_l)
                if x_overlap and abs(y_final - py) <= collision_margin:
                    floors.append(py + min_spacing)
            if floors:
                y_final = max(y_final, max(floors))

            if prev_y is not None:
                y_final = max(y_final, prev_y + min_spacing)

            b['y_final'] = y_final
            placed.append((b['x_left'], b['x_right'], y_final))
            prev_y = y_final

        # Expand y-axis to fit all brackets + text headroom
        if brackets:
            max_y = max(b['y_final'] for b in brackets)
            needed = max_y + (0.12 + abs(text_gap)) * y_range
            if needed > ax.get_ylim()[1]:
                ax.set_ylim(y_lo, needed)
                y_lo, y_hi = ax.get_ylim()
                y_range = y_hi - y_lo

        # --- Draw ----------------------------------------------------------------
        for b in brackets:
            y_line = b['y_final']
            label  = b['label']

            if b.get('suppressed', False):
                continue

            ax.hlines(y_line, b['x_left'], b['x_right'],
                      colors=line_color, linewidth=line_width, zorder=5)

            fs = star_fontsize if label not in ("ns", "nd") else other_fontsize
            y_text = y_line + (star_gap if label not in ("ns", "nd") else text_gap) * y_range
            ax.text(b['x_mid'], y_text, label,
                    ha="center", va="bottom", fontsize=fs, zorder=5)

    if show_n:
        y_lo_cur = ax.get_ylim()[0]
        for xi, x_val in enumerate(x_order):
            for hi, hue_val in enumerate(hue_order):
                x_center = xi + (hi - (n_hue - 1) / 2) * dodge
                x_vals = df.loc[(df[x_col] == x_val) & (df[hue_col] == hue_val), y_col].dropna().values
                ax.text(x_center, y_lo_cur - n_offset * y_range, f"n={len(x_vals)}",
                        ha="center", va="top", fontsize=n_fontsize,
                        color="dimgray", zorder=5)


# ── Backward-compatible aliases ───────────────────────────────────────────────
def annotate_mwu(ax, df, x_col, y_col, hue_col, x_order, hue_order, **kwargs):
    kwargs.setdefault("test", "mwu")
    if len(hue_order) <= 2:
        return annotate_significance(ax, df, x_col, y_col, hue_col, x_order, hue_order, **kwargs)

    # For 3+ hue groups, compare all groups against a reference hue.
    # Priority: explicit ref_hue in kwargs, then WT/MMY if present, else first hue.
    ref_hue = kwargs.pop("ref_hue", None)
    if ref_hue is None:
        for candidate in ("WT", "MMY"):
            if candidate in hue_order:
                ref_hue = candidate
                break
    if ref_hue is None:
        ref_hue = hue_order[0]

    return annotate_significance_hue_vs_ref(
        ax,
        df,
        x_col,
        y_col,
        hue_col,
        x_order,
        hue_order,
        ref_hue,
        **kwargs,
    )

def annotate_mwu_vs_ref(ax, df, x_col, y_col, x_order, ref_group, **kwargs):
    kwargs.setdefault("test", "mwu")
    return annotate_significance_vs_ref(ax, df, x_col, y_col, x_order, ref_group, **kwargs)

def annotate_mwu_hue_vs_ref(ax, df, x_col, y_col, hue_col, x_order, hue_order, ref_hue, **kwargs):
    kwargs.setdefault("test", "mwu")
    return annotate_significance_hue_vs_ref(ax, df, x_col, y_col, hue_col, x_order, hue_order, ref_hue, **kwargs)


def annotate_chi2_proportions_vs_ref(
    ax,
    df,
    x_col,
    outcome_col,
    x_order,
    ref_group,
    success_values=(1, True),
    failure_values=(0, False),
    y_text=1.03,
    y_lim=1.12,
    star_fontsize=16,
    other_fontsize=11,
):
    """
    Annotate chi-square significance stars on a proportion bar plot.

    Compares each group in ``x_order`` to ``ref_group`` using a 2x2 contingency
    table built from ``outcome_col`` (success vs failure counts).

    Notes
    -----
    - This helper is intended for plots where x positions are categorical
      indices matching ``x_order`` (e.g. seaborn histplot/barplot by category).
    - Stars are rendered above the non-reference groups.
    """

    def _p_to_stars(p):
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    if ref_group not in x_order:
        raise ValueError(f"ref_group '{ref_group}' is not present in x_order")

    # Ensure headroom for text while preserving larger existing limits.
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, max(y1, y_lim))

    ref_counts = df.loc[df[x_col] == ref_group, outcome_col].value_counts()
    ref_success = int(sum(ref_counts.get(v, 0) for v in success_values))
    ref_failure = int(sum(ref_counts.get(v, 0) for v in failure_values))

    if (ref_success + ref_failure) == 0:
        return

    for i, group in enumerate(x_order):
        if group == ref_group:
            continue

        grp_counts = df.loc[df[x_col] == group, outcome_col].value_counts()
        grp_success = int(sum(grp_counts.get(v, 0) for v in success_values))
        grp_failure = int(sum(grp_counts.get(v, 0) for v in failure_values))

        if (grp_success + grp_failure) == 0:
            label = "nd"
        else:
            _, p, _, _ = chi2_contingency(
                [[ref_success, ref_failure], [grp_success, grp_failure]]
            )
            label = _p_to_stars(p)

        fs = star_fontsize if label not in ("ns", "nd") else other_fontsize
        ax.text(i, y_text, label, ha="center", va="bottom", fontsize=fs)
        
def sig_test_for_bs(xa, xb):
    """
    Bootstrap significance test based on confidence intervals of the difference.
    Returns p-value corresponding to the narrowest CI that does NOT contain 0.
    """
    xa_minus_xb = xa - xb
    
    # Check narrowest CI (95%)
    p = None
    l95, h95 = np.percentile(xa_minus_xb, [2.5, 97.5])
    if l95 <= 0 <= h95 and p is None:
        p = 0.055 # ns, 0 is in 95% CI

    # Check 99%
    l99, h99 = np.percentile(xa_minus_xb, [0.5, 99.5])
    if l99 <= 0 <= h99 and p is None:
        p = 0.015 # 0 not in 95% but is in 99%
    
    # Check 99.9%
    l999, h999 = np.percentile(xa_minus_xb, [0.05, 99.95])
    if l999 <= 0 <= h999 and p is None:
        p = 0.0015 # 0 is in 99.9% CI, but not in 99% CI
    
    # 0 is not even in 99.9% CI, so p < 0.001
    if p is None:
        p = 0.0005
    print(f"95%: [{l95:.3f}, {h95:.3f}], 99%: [{l99:.3f}, {h99:.3f}], 99.9%: [{l999:.3f}, {h999:.3f}]")
    print(f"Pred. p-value: {p}, stars: {_p_to_stars(p)}")
    return p