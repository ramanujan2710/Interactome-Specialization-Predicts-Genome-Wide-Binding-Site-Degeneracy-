"""Shared statistical helpers used across pipeline stages."""
import numpy as np
from scipy.stats import rankdata, norm


def rankz(x):
    """Rank-transform then z-score, i.e. a normal-scores transform."""
    r = rankdata(x)
    return (r - r.mean()) / r.std()


def resid(y, X):
    """Residualize y on X (+ intercept) via ordinary least squares."""
    A = np.column_stack([X, np.ones(len(y))])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ b


def combine_probs(*channels):
    """STRING-style probabilistic-OR combination of independent evidence
    channels, each given as a 0-1000 integer score column (pandas Series or
    ndarray): 1 - prod(1 - p_i)."""
    p = np.ones(len(channels[0]), dtype=np.float64)
    for ch in channels:
        vals = ch.values if hasattr(ch, "values") else np.asarray(ch)
        p *= (1 - vals.astype(np.float64) / 1000.0)
    return 1 - p


def fisher_z_test(r1, n1, r2, n2):
    """Fisher z-test for the difference between two independent Pearson
    correlation coefficients."""
    z1, z2 = np.arctanh(r1), np.arctanh(r2)
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z = (z1 - z2) / se
    p = 2 * (1 - norm.cdf(abs(z)))
    return z, p


def post_hoc_power(r, n, alpha=0.05):
    """Post-hoc power of a two-sided Pearson correlation test at the given
    observed effect size and sample size."""
    z_r = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = norm.ppf(1 - alpha / 2)
    return 1 - norm.cdf(z_crit - z_r / se) + norm.cdf(-z_crit - z_r / se)


def winsorize(x, lo=0.05, hi=0.95):
    """Clip a pandas Series to its [lo, hi] quantile range."""
    lo_v, hi_v = x.quantile(lo), x.quantile(hi)
    return x.clip(lo_v, hi_v)
