"""
ddpcr_quant.py — ddPCR quantification
======================================
Implements BioRad's quantification pipeline confirmed from MetricsComputation source.

Concentration formula (MetricsComputation.ComputeSampleMetrics /
                        Well.ComputeCopiesPerDropletPerChannel):
    cpd = -ln(negatives / (positives + negatives))
    concentration (copies/µL) = cpd / (droplet_volume_nL / 1000)

Positive/negative counting (ClusteringInfo.GetPositivesAndNegativesForAllTargets):
    For standard 2-channel duplex, Cluster integer values 0–4:
        Ch1 positives = count(Cluster 2: Ch1+) + count(Cluster 3: Double+)
        Ch1 negatives = count(Cluster 1: NN)   + count(Cluster 4: Ch2+)
        Ch2 positives = count(Cluster 3: Double+) + count(Cluster 4: Ch2+)
        Ch2 negatives = count(Cluster 1: NN)      + count(Cluster 2: Ch1+)
        Cluster 0 (Gated/Filtered) is excluded from all counts.

Confidence intervals: Clopper-Pearson exact binomial on the negative fraction
p0 = negatives / total, converting bounds back to cpd = -ln(p0).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ── Cluster integer encoding (same for QLP binary and ddPCR JSON) ─────────────
#   0 = Gated/Filtered   1 = NN (Ch1−Ch2−)
#   2 = Ch1+ (Ch1+Ch2−)  3 = Double+ (Ch1+Ch2+)  4 = Ch2+ (Ch1−Ch2+)

# BioRad QX Manager display colours (approximated)
CLUSTER_COLORS: Dict[int, str] = {
    0: '#d0d0d0',   # Gated: light grey
    1: '#808080',   # NN: medium grey
    2: '#1560bd',   # Ch1+: blue
    3: '#009944',   # Double+: green
    4: '#e06000',   # Ch2+: orange
}


def cluster_labels(ch1_name: str = 'Ch1', ch2_name: str = 'Ch2') -> Dict[int, str]:
    """Human-readable cluster labels using actual channel names."""
    return {
        0: 'Gated',
        1: f'{ch1_name}−  {ch2_name}−',
        2: f'{ch1_name}+  {ch2_name}−',
        3: f'{ch1_name}+  {ch2_name}+',
        4: f'{ch1_name}−  {ch2_name}+',
    }


# ── Per-channel quantification result ─────────────────────────────────────────

@dataclass
class ChannelQuant:
    name: str
    positives: int
    negatives: int
    threshold: Optional[float] = None       # amplitude threshold used for classification
    copies_per_droplet: Optional[float] = None
    concentration: Optional[float] = None   # copies/µL
    ci_lower_95: Optional[float] = None     # 95% CI lower bound (copies/µL)
    ci_upper_95: Optional[float] = None     # 95% CI upper bound (copies/µL)

    @property
    def total_accepted(self) -> int:
        return self.positives + self.negatives

    @property
    def positive_fraction(self) -> Optional[float]:
        t = self.total_accepted
        return self.positives / t if t > 0 else None


# ── Per-well quantification result ────────────────────────────────────────────

@dataclass
class WellQuant:
    well_id: str
    channels: List[ChannelQuant]
    cluster_counts: Dict[int, int]          # cluster_id → droplet count
    droplet_volume_nL: float = 0.85
    total_events: int = 0                   # all droplets including gated
    has_biorad_clusters: bool = False       # True for QLP / ddPCR sources

    @property
    def gated_count(self) -> int:
        return self.cluster_counts.get(0, 0)

    @property
    def accepted_count(self) -> int:
        return sum(self.cluster_counts.get(k, 0) for k in [1, 2, 3, 4])

    def ch(self, index: int) -> Optional[ChannelQuant]:
        return self.channels[index] if index < len(self.channels) else None


# ── Positive/negative counting ────────────────────────────────────────────────

def pos_neg_from_cluster_counts(cluster_counts: Dict[int, int],
                                 n_channels: int = 2) -> List[Tuple[int, int]]:
    """
    Derive per-channel positive and negative counts from cluster integer counts.
    Standard 2-channel duplex only — use pos_neg_from_labels for multiplex.

    Confirmed from ClusteringInfo.GetPositivesAndNegativesForAllTargets:
        Ch1 positives = count(2: Ch1+) + count(3: Double+)
        Ch1 negatives = count(1: NN)   + count(4: Ch2+)
        Ch2 positives = count(3: Double+) + count(4: Ch2+)
        Ch2 negatives = count(1: NN)      + count(2: Ch1+)
    """
    c = {i: cluster_counts.get(i, 0) for i in range(5)}
    if n_channels >= 2:
        return [
            (c[2] + c[3], c[1] + c[4]),
            (c[3] + c[4], c[1] + c[2]),
        ]
    else:
        return [(c[2] + c[3], c[1] + c[4])]


def pos_neg_from_labels(df: pd.DataFrame,
                         channel_names: List[str]) -> List[Tuple[int, int]]:
    """
    Derive per-channel positive and negative counts using the Cluster_Label
    column, which contains dye-specific strings like 'FAM+  HEX−'.

    Works for both duplex (2-channel) and amplitude-multiplex (N-channel) assays
    because the label encodes the dye name directly.

    A droplet is positive for channel N if its label contains '{dye_name}+'
    A droplet is negative for channel N if its label contains '{dye_name}−'
      (includes both the '−' Unicode minus and the plain '-' ASCII minus)
    Gated droplets (Cluster == 0) are always excluded.

    Falls back to integer-based counting (pos_neg_from_cluster_counts) if the
    Cluster_Label column is absent or uses generic 'Ch1+'/'Ch2+' labels.
    """
    has_labels = ('Cluster_Label' in df.columns and
                  df['Cluster_Label'].notna().any())

    # Check if labels encode actual dye names (not generic Ch1+/Ch2+)
    if has_labels:
        sample_labels = df['Cluster_Label'].dropna().unique()
        generic = all(lbl in ('Gated', 'NN', 'Ch1+', 'Ch1+Ch2+', 'Ch2+',
                               'Ch1+  Ch2+', 'Ch1+  Ch2−', 'Ch1−  Ch2+',
                               'Ch1−  Ch2−', 'Filtered', 'Unassigned', 'Unknown')
                      for lbl in sample_labels)
        if generic:
            has_labels = False   # fall through to integer-based

    if not has_labels or 'Cluster' not in df.columns:
        # Fallback: integer-based for standard 2-channel duplex
        counts = df['Cluster'].value_counts().to_dict() if 'Cluster' in df.columns else {}
        counts = {int(k): int(v) for k, v in counts.items()}
        return pos_neg_from_cluster_counts(counts, len(channel_names))

    accepted = ~df['Cluster_Label'].isin(['Gated', 'Filtered', 'Unassigned'])
    if 'Cluster' in df.columns:
        accepted = accepted & (df['Cluster'] != 0)

    result = []
    for name in channel_names:
        pos_pat_u = f'{name}+'          # Unicode context e.g. 'FAM+  HEX−'
        # Count positives: label contains '{name}+'
        pos_mask = accepted & df['Cluster_Label'].str.contains(
            pos_pat_u, regex=False, na=False)
        # Count negatives: accepted and NOT positive for this channel
        neg_mask = accepted & ~df['Cluster_Label'].str.contains(
            pos_pat_u, regex=False, na=False)
        result.append((int(pos_mask.sum()), int(neg_mask.sum())))

    return result


# ── Core Poisson formulas ──────────────────────────────────────────────────────

def poisson_cpd(positives: int, negatives: int) -> Optional[float]:
    """
    Copies per droplet = -ln(negatives / total).
    Confirmed from Well.ComputeCopiesPerDropletPerChannel.
    Returns None if total == 0, inf if negatives == 0.
    """
    total = positives + negatives
    if total == 0:
        return None
    if negatives == 0:
        return float('inf')
    return -math.log(negatives / total)


def poisson_concentration(positives: int, negatives: int,
                           droplet_volume_nL: float) -> Optional[float]:
    """
    Concentration in copies/µL.
    Confirmed from MetricsComputation.ComputeSampleMetrics:
        num9 = droplet_volume_nL / 1000   (nL → µL)
        conc = cpd / num9
    """
    cpd = poisson_cpd(positives, negatives)
    if cpd is None or math.isinf(cpd) or droplet_volume_nL <= 0:
        return None
    return cpd / (droplet_volume_nL / 1000.0)


# ── Confidence intervals — BioRad's exact algorithm ──────────────────────────
#
# Confirmed from DDPcrCalculations.GetPoissonConfidenceCpdBands source, with
# the exact lookup tables from BayesianConfidenceIntervals (200 entries each).
#
# Table structure: 200 entries, paired as (lower_CI_count, upper_CI_count)
# for k = 0..99.  Entry [2k] = lower 95% CI bound on count, [2k+1] = upper.
# These are absolute counts (not fractions) — divided by (pos+neg) at call site.
#
# Three branches (exact source match):
#
# Branch 1 (pos < 100):
#   p_pos_lo = table[pos*2] / total  →  cpd_lo = -ln(1 - p_pos_lo)
#   p_pos_hi = table[pos*2+1] / total  →  cpd_hi = -ln(1 - p_pos_hi)
#
# Branch 2 (neg < 100):
#   lower: num6 = (1 - table[neg*2+1]/total) * z  →  cpd_lo = -ln(1-num6)
#     NOTE: z=1.96 for 95% makes num6 > 1 for any realistic total, producing
#     NaN via log(negative). This is a confirmed bug in BioRad's code.
#     For 68%, z=1.0 so this branch works. We replicate faithfully.
#   upper: num6 = 1 - table[neg*2]/total  →  cpd_hi = -ln(1-num6)   (no z!)
#
# Branch 3 (pos >= 100 AND neg >= 100): normal approximation on p_positive.
#   (exact match documented in previous implementation)
#
# The cpd - num4 / cpd + num5 lines in the source are dead code (identity
# operations) and are omitted here.

# BioRad constants
_Z_95 = 1.959963985   # c_95PctBandSigma  ≈ scipy.stats.norm.ppf(0.975)
_Z_68 = 1.0           # c_68PctBandSigma
_MAX_CPD = math.log(1_000_000)

# Exact tables from BayesianConfidenceIntervals — 200 entries each.
# Index [2k] = lower CI count, [2k+1] = upper CI count, for k=0..99.
_CI95 = [
    0.0, 2.996, 0.042, 4.776, 0.303, 6.406, 0.712, 7.951, 1.207, 9.431,
    1.758, 10.865, 2.35, 12.264, 2.974, 13.632, 3.622, 14.979, 4.292, 16.304,
    4.979, 17.614, 5.68, 18.911, 6.395, 20.193, 7.121, 21.465, 7.858, 22.726,
    8.603, 23.98, 9.356, 25.226, 10.117, 26.464, 10.885, 27.694, 11.659, 28.919,
    12.439, 30.139, 13.224, 31.353, 14.015, 32.561, 14.809, 33.767, 15.608, 34.967,
    16.412, 36.162, 17.219, 37.355, 18.03, 38.543, 18.844, 39.729, 19.661, 40.912,
    20.482, 42.09, 21.306, 43.266, 22.132, 44.439, 22.961, 45.61, 23.792, 46.779,
    24.626, 47.945, 25.462, 49.109, 26.301, 50.269, 27.142, 51.428, 27.984, 52.586,
    28.829, 53.741, 29.676, 54.893, 30.524, 56.045, 31.375, 57.194, 32.227, 58.342,
    33.08, 59.489, 33.935, 60.634, 34.792, 61.776, 35.651, 62.917, 36.51, 64.058,
    37.372, 65.196, 38.234, 66.334, 39.098, 67.47, 39.963, 68.605, 40.829, 69.739,
    41.697, 70.87, 42.566, 72.001, 43.436, 73.131, 44.307, 74.26, 45.18, 75.387,
    46.053, 76.514, 46.927, 77.64, 47.803, 78.764, 48.679, 79.888, 49.556, 81.011,
    50.434, 82.133, 51.314, 83.252, 52.194, 84.372, 53.075, 85.491, 53.957, 86.609,
    54.84, 87.726, 55.723, 88.843, 56.608, 89.958, 57.493, 91.073, 58.379, 92.187,
    59.266, 93.3, 60.153, 94.413, 61.041, 95.525, 61.93, 96.636, 62.819, 97.747,
    63.71, 98.856, 64.601, 99.964, 65.493, 101.072, 66.385, 102.18, 67.278, 103.287,
    68.172, 104.393, 69.066, 105.499, 69.961, 106.604, 70.856, 107.709, 71.752, 108.813,
    72.649, 109.916, 73.546, 111.019, 74.444, 112.121, 75.342, 113.223, 76.241, 114.324,
    77.14, 115.425, 78.04, 116.525, 78.941, 117.624, 79.842, 118.723, 80.743, 119.822,
]

_CI68 = [
    0.0, 1.369, 0.271, 2.489, 0.868, 3.845, 1.559, 5.136, 2.297, 6.389,
    3.066, 7.615, 3.857, 8.82, 4.663, 10.012, 5.484, 11.189, 6.315, 12.356,
    7.155, 13.515, 8.002, 14.667, 8.857, 15.811, 9.717, 16.951, 10.582, 18.085,
    11.453, 19.214, 12.327, 20.339, 13.205, 21.461, 14.087, 22.579, 14.972, 23.693,
    15.86, 24.805, 16.751, 25.914, 17.644, 27.02, 18.539, 28.125, 19.437, 29.227,
    20.337, 30.327, 21.239, 31.425, 22.143, 32.521, 23.048, 33.616, 23.955, 34.709,
    24.864, 35.799, 25.774, 36.889, 26.686, 37.977, 27.599, 39.064, 28.513, 40.15,
    29.429, 41.234, 30.346, 42.317, 31.264, 43.399, 32.182, 44.481, 33.102, 45.561,
    34.024, 46.638, 34.946, 47.716, 35.869, 48.793, 36.793, 49.869, 37.717, 50.945,
    38.643, 52.019, 39.569, 53.093, 40.497, 54.165, 41.425, 55.237, 42.354, 56.308,
    43.283, 57.379, 44.213, 58.449, 45.144, 59.518, 46.076, 60.586, 47.008, 61.654,
    47.941, 62.721, 48.874, 63.788, 49.808, 64.854, 50.742, 65.92, 51.677, 66.985,
    52.613, 68.049, 53.549, 69.113, 54.486, 70.175, 55.423, 71.238, 56.361, 72.3,
    57.299, 73.362, 58.237, 74.425, 59.177, 75.484, 60.116, 76.545, 61.056, 77.605,
    61.997, 78.664, 62.937, 79.724, 63.879, 80.782, 64.82, 81.841, 65.763, 82.898,
    66.705, 83.956, 67.648, 85.013, 68.591, 86.07, 69.535, 87.126, 70.479, 88.182,
    71.423, 89.238, 72.368, 90.293, 73.313, 91.348, 74.258, 92.403, 75.204, 93.457,
    76.15, 94.511, 77.096, 95.565, 78.042, 96.619, 78.989, 97.672, 79.937, 98.724,
    80.884, 99.777, 81.832, 100.829, 82.78, 101.881, 83.728, 102.933, 84.677, 103.984,
    85.626, 105.035, 86.575, 106.086, 87.525, 107.136, 88.474, 108.187, 89.424, 109.237,
]


def biorad_cpd_ci(positives: int, negatives: int,
                   alpha: float = 0.05) -> Tuple[Optional[float], Optional[float]]:
    """
    CI on copies-per-droplet using BioRad's exact algorithm.
    Exact port of DDPcrCalculations.GetPoissonConfidenceCpdBands using the
    BayesianConfidenceIntervals lookup tables verbatim.

    alpha=0.05 → 95% CI (uses _CI95 table, z=1.96)
    alpha=0.32 → 68% CI (uses _CI68 table, z=1.0)
    Returns (cpd_lower, cpd_upper).
    """
    pos, neg = float(positives), float(negatives)
    total = pos + neg
    if total == 0:
        return None, None

    # cpd = -ln(neg/total); NaN if either is NaN (not possible from int inputs)
    if neg == 0:
        cpd = float('inf')
    else:
        cpd = -math.log(neg / total)

    is_95 = (alpha <= 0.05)
    table = _CI95 if is_95 else _CI68
    z = _Z_95 if is_95 else _Z_68

    cpd_lo = float('nan')
    cpd_hi = float('nan')

    if pos < 100:
        # Branch 1: Bayesian table indexed on positives
        idx_lo = int(pos * 2.0 + 0.5)       # = pos*2 for integer pos
        idx_hi = int(pos * 2.0 + 1.0 + 0.5) # = pos*2+1 for integer pos

        p_pos_lo = table[idx_lo] / total
        p_pos_hi = table[idx_hi] / total

        if cpd == 0.0:
            cpd_lo = 0.0
        else:
            val = 1.0 - p_pos_lo
            cpd_lo = -math.log(val) if val > 0 else float('inf')

        if p_pos_hi >= 1.0:
            cpd_hi = _MAX_CPD
        else:
            val = 1.0 - p_pos_hi
            cpd_hi = -math.log(val) if val > 0 else _MAX_CPD

    elif neg < 100:
        # Branch 2: Bayesian table indexed on negatives.
        # NOTE: the *z multiplication on the lower bound is present in BioRad's
        # source and causes num6 > 1 for 95% CI with any realistic total, producing
        # NaN via log(negative). This is faithfully replicated here.
        idx_neg_hi = int(neg * 2.0 + 1.0 + 0.5)  # = neg*2+1 for integer neg
        idx_neg_lo = int(neg * 2.0 + 0.5)          # = neg*2 for integer neg

        # Lower cpd bound
        num6 = (1.0 - table[idx_neg_hi] / total) * z
        if num6 > 0.0:
            val = 1.0 - num6
            # For 95% CI, val is typically negative → log(negative) → NaN
            cpd_lo = -math.log(val) if val > 0 else float('nan')
        else:
            cpd_lo = 0.0

        # Upper cpd bound (no z multiplication — confirmed from source)
        num6 = 1.0 - table[idx_neg_lo] / total
        if num6 == 1.0:
            cpd_hi = _MAX_CPD
        else:
            val = 1.0 - num6  # = table[idx_neg_lo] / total
            cpd_hi = -math.log(val) if val > 0 else _MAX_CPD

    else:
        # Branch 3: Normal approximation on p_positive (large counts, exact match)
        p_pos = pos / total
        sigma = math.sqrt(p_pos * (1.0 - p_pos) / total)
        p_pos_lo = p_pos - z * sigma
        p_pos_hi = p_pos + z * sigma

        p_neg_hi = 1.0 - p_pos_lo
        p_neg_lo = 1.0 - p_pos_hi

        cpd_lo = -math.log(p_neg_hi) if 0 < p_neg_hi < 1 else 0.0
        cpd_hi = -math.log(p_neg_lo) if 0 < p_neg_lo < 1 else _MAX_CPD

    # Guard rails — exact match to source (NaN comparisons are false; preserved)
    if cpd_lo < 0.0:
        cpd_lo = 0.0
    if cpd_hi > _MAX_CPD:
        cpd_hi = _MAX_CPD
    if pos == 0.0:
        cpd_lo = 0.0
    if neg == 0.0:
        cpd_hi = _MAX_CPD

    return cpd_lo, cpd_hi


def concentration_ci(positives: int, negatives: int,
                      droplet_volume_nL: float,
                      alpha: float = 0.05) -> Tuple[Optional[float], Optional[float]]:
    """95% CI on concentration (copies/µL) using BioRad's algorithm."""
    if droplet_volume_nL <= 0:
        return None, None
    vol_uL = droplet_volume_nL / 1000.0
    cpd_lo, cpd_hi = biorad_cpd_ci(positives, negatives, alpha)
    ci_lo = (cpd_lo / vol_uL) if cpd_lo is not None else None
    ci_hi = (cpd_hi / vol_uL) if (cpd_hi is not None and cpd_hi < _MAX_CPD) else None
    return ci_lo, ci_hi


# ── Threshold-based cluster assignment (fallback for CSV files) ────────────────

def assign_clusters_from_thresholds(df: pd.DataFrame,
                                      channel_names: List[str],
                                      thresholds: Dict[str, float]) -> pd.Series:
    """
    Assign Cluster integers 0–4 using amplitude thresholds.
    Used as a fallback when BioRad cluster assignments are not available (CSV files).

    Returns a Series of cluster integers aligned with df's index.
    All droplets start as Cluster 1 (NN); Cluster 0 is reserved for droplets
    flagged by Gating_Flag != 0 if that column is present.
    """
    ch1_name = channel_names[0] if channel_names else 'Ch1'
    ch2_name = channel_names[1] if len(channel_names) > 1 else 'Ch2'

    def _col(name: str) -> Optional[str]:
        for candidate in [f'{name}_Amplitude', 'Ch1_Amplitude', 'Ch2_Amplitude']:
            if candidate in df.columns:
                if 'Ch1' in candidate and name == ch1_name:
                    return candidate
                if 'Ch2' in candidate and name == ch2_name:
                    return candidate
                if candidate == f'{name}_Amplitude':
                    return candidate
        return None

    ch1_col = _col(ch1_name) or ('Ch1_Amplitude' if 'Ch1_Amplitude' in df.columns else None)
    ch2_col = _col(ch2_name) or ('Ch2_Amplitude' if 'Ch2_Amplitude' in df.columns else None)

    t1 = thresholds.get(ch1_name) or thresholds.get('Ch1')
    t2 = thresholds.get(ch2_name) or thresholds.get('Ch2')

    clusters = pd.Series(1, index=df.index, dtype=int)  # default NN

    # Mark gated
    if 'Gating_Flag' in df.columns:
        clusters[df['Gating_Flag'] != 0] = 0

    if ch1_col and ch2_col and t1 is not None and t2 is not None:
        ch1_pos = df[ch1_col] >= t1
        ch2_pos = df[ch2_col] >= t2
        accepted = clusters != 0
        clusters[accepted & ~ch1_pos & ~ch2_pos] = 1
        clusters[accepted &  ch1_pos & ~ch2_pos] = 2
        clusters[accepted &  ch1_pos &  ch2_pos] = 3
        clusters[accepted & ~ch1_pos &  ch2_pos] = 4
    elif ch1_col and t1 is not None:
        ch1_pos = df[ch1_col] >= t1
        accepted = clusters != 0
        clusters[accepted & ~ch1_pos] = 1
        clusters[accepted &  ch1_pos] = 2

    return clusters


def detect_thresholds_1d(values: np.ndarray,
                          min_separation: float = 500.0) -> List[float]:
    """
    Density-based 1D threshold detection for CSV files without BioRad clusters.
    Uses log-space KDE peak detection to find valleys between populations.
    Returns a list of threshold values (may be empty if only one population).
    """
    from scipy.signal import find_peaks
    from scipy.stats import gaussian_kde

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 10:
        return []

    val_range = values.max() - values.min()
    if val_range < min_separation:
        return []

    try:
        grid = np.linspace(values.min(), values.max(), 500)
        kde = gaussian_kde(values, bw_method='scott')
        density = kde(grid)
    except Exception:
        return []

    log_density = np.log(density + 1e-10)
    log_range = log_density.max() - log_density.min()
    peaks, _ = find_peaks(log_density, prominence=0.05 * log_range)

    if len(peaks) < 2:
        return []

    thresholds = []
    for i in range(len(peaks) - 1):
        valley_slice = density[peaks[i]:peaks[i + 1]]
        local_min_idx = np.argmin(valley_slice) + peaks[i]
        raw_thresh = float(grid[local_min_idx])

        # Safety: at least 4 median-absolute-deviations above the left peak
        median = np.median(values)
        mad = np.median(np.abs(values - median))
        sigma_est = 1.4826 * mad
        safe_thresh = float(grid[peaks[i]]) + 4 * sigma_est
        thresholds.append(max(raw_thresh, safe_thresh))

    # Deduplicate thresholds that are too close together
    result = []
    for t in sorted(thresholds):
        if not result or (t - result[-1]) > min_separation:
            result.append(t)
    return result


# ── Main entry point ──────────────────────────────────────────────────────────

def get_thresholds_from_metadata(df: pd.DataFrame,
                                   channel_names: List[str]) -> Dict[str, float]:
    """
    Extract thresholds from df.attrs['well_metadata'] where available.
    Handles both QLP naming (threshold_ch1) and ddPCR naming (threshold_FAM).
    """
    meta = df.attrs.get('well_metadata', {})
    thresholds: Dict[str, float] = {}

    for i, name in enumerate(channel_names):
        ch_label = f'ch{i + 1}'
        # QLP: threshold_ch1 / threshold_ch2
        v = meta.get(f'threshold_{ch_label}')
        # ddPCR: threshold_FAM / threshold_HEX (dye name)
        if v is None:
            v = meta.get(f'threshold_{name}')
        if v is not None and np.isfinite(float(v)) and float(v) > 0:
            thresholds[name] = float(v)

    return thresholds


def compute_well_quant(df: pd.DataFrame,
                        channel_names: Optional[List[str]] = None,
                        droplet_volume_nL: float = 0.85,
                        shared_thresholds: Optional[Dict[str, float]] = None
                        ) -> WellQuant:
    """
    Compute quantification results for a single well DataFrame.

    Priority order for cluster assignment:
      1. BioRad Cluster column (QLP / ddPCR files) — use directly.
      2. Thresholds from well_metadata (QLP / ddPCR with metadata).
      3. Shared thresholds passed in (for cross-well normalisation).
      4. Auto-detected thresholds from amplitude distributions (CSV fallback).

    Parameters
    ----------
    df              Per-droplet DataFrame from a parser.
    channel_names   Override channel names; defaults to df.attrs['channel_names'].
    droplet_volume_nL  Droplet volume in nanolitres.
    shared_thresholds  Pre-computed thresholds to use if no BioRad clusters.
    """
    if channel_names is None:
        channel_names = df.attrs.get('channel_names', ['Ch1', 'Ch2'])

    well_id = df.attrs.get('well_metadata', {}).get('well_id', 'Unknown')

    # ── Resolve droplet volume ────────────────────────────────────────────
    meta_vol = df.attrs.get('well_metadata', {}).get('droplet_volume_nL')
    if meta_vol and float(meta_vol) > 0:
        droplet_volume_nL = float(meta_vol)

    # ── Determine cluster assignments ─────────────────────────────────────
    has_biorad = (
        'Cluster' in df.columns
        and df['Cluster'].dtype.kind in ('i', 'u', 'f')
        and df['Cluster'].isin([1, 2, 3, 4]).any()
    )

    if has_biorad:
        cluster_series = df['Cluster'].fillna(0).astype(int)
    else:
        # Resolve thresholds for CSV / unclassified files
        thresholds = get_thresholds_from_metadata(df, channel_names)
        if not thresholds and shared_thresholds:
            thresholds = shared_thresholds
        if not thresholds:
            # Auto-detect per channel
            thresholds = {}
            for i, name in enumerate(channel_names[:2]):
                col = f'{name}_Amplitude' if f'{name}_Amplitude' in df.columns \
                      else (f'Ch{i+1}_Amplitude' if f'Ch{i+1}_Amplitude' in df.columns else None)
                if col:
                    detected = detect_thresholds_1d(df[col].values)
                    if detected:
                        thresholds[name] = detected[0]
        cluster_series = assign_clusters_from_thresholds(df, channel_names, thresholds)

    cluster_counts = cluster_series.value_counts().to_dict()
    cluster_counts = {int(k): int(v) for k, v in cluster_counts.items()}

    # ── Count positives/negatives per channel ─────────────────────────────
    # Use label-based counting (works for both duplex and amplitude multiplex).
    # For standard duplex this is equivalent to pos_neg_from_cluster_counts.
    n_channels = len(channel_names)
    pos_neg_list = pos_neg_from_labels(df, channel_names)
    # Pad with (0,0) if fewer results than channels (shouldn't happen)
    while len(pos_neg_list) < n_channels:
        pos_neg_list.append((0, 0))

    # ── Resolve displayed thresholds (for plot annotations) ───────────────
    thresh_map = get_thresholds_from_metadata(df, channel_names)
    if not thresh_map and shared_thresholds:
        thresh_map = shared_thresholds

    # ── Build per-channel results ─────────────────────────────────────────
    channels: List[ChannelQuant] = []
    for i, name in enumerate(channel_names):
        pos, neg = pos_neg_list[i]
        thresh = thresh_map.get(name)

        cpd = poisson_cpd(pos, neg)
        conc = poisson_concentration(pos, neg, droplet_volume_nL)
        ci_lo, ci_hi = concentration_ci(pos, neg, droplet_volume_nL)

        channels.append(ChannelQuant(
            name=name,
            positives=pos,
            negatives=neg,
            threshold=thresh,
            copies_per_droplet=cpd,
            concentration=conc,
            ci_lower_95=ci_lo,
            ci_upper_95=ci_hi,
        ))

    return WellQuant(
        well_id=well_id,
        channels=channels,
        cluster_counts=cluster_counts,
        droplet_volume_nL=droplet_volume_nL,
        total_events=len(df),
        has_biorad_clusters=has_biorad,
    )


def well_quant_to_stats_row(wq: WellQuant, condition: str = '') -> Dict:
    """Flatten a WellQuant into a single CSV stats row."""
    row: Dict = {
        'Well': wq.well_id,
        'Condition': condition,
        'Total_Droplets': wq.total_events,
        'Accepted_Droplets': wq.accepted_count,
        'Gated_Droplets': wq.gated_count,
        'Droplet_Volume_nL': wq.droplet_volume_nL,
        'Source': 'BioRad' if wq.has_biorad_clusters else 'Threshold',
    }
    # Per-cluster counts
    for cid, label in {1: 'NN', 2: 'Ch1pos', 3: 'Double', 4: 'Ch2pos'}.items():
        row[f'Count_{label}'] = wq.cluster_counts.get(cid, 0)

    # Per-channel stats
    for ch in wq.channels:
        n = ch.name
        row[f'{n}_Positives'] = ch.positives
        row[f'{n}_Negatives'] = ch.negatives
        row[f'{n}_Threshold'] = ch.threshold
        row[f'{n}_CopiesPerDroplet'] = ch.copies_per_droplet
        row[f'{n}_Concentration_copies_uL'] = ch.concentration
        row[f'{n}_CI_Lower_95'] = ch.ci_lower_95
        row[f'{n}_CI_Upper_95'] = ch.ci_upper_95

    return row