"""
Lightweight sanity tests for the shared statistical helpers (src/tf_ppi/utils.py).
These use synthetic data only -- no external downloads required -- and are
meant as a fast regression check, not a substitute for the Stage 14
adversarial checks that run against the real analysis output.

Run with:  pytest tests/  (or: python -m pytest tests/)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pandas as pd

from tf_ppi.utils import rankz, resid, combine_probs, fisher_z_test, post_hoc_power, winsorize


def test_rankz_mean_zero_unit_std():
    x = np.array([5, 1, 3, 2, 4], dtype=float)
    z = rankz(x)
    assert abs(z.mean()) < 1e-9
    assert abs(z.std() - 1.0) < 1e-9


def test_resid_orthogonal_to_regressor():
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    y = 2 * x + rng.normal(size=200)
    r = resid(y, x)
    # residuals should be (numerically) uncorrelated with the regressor
    assert abs(np.corrcoef(r, x)[0, 1]) < 1e-8


def test_combine_probs_bounds_and_monotonicity():
    a = np.array([0, 500, 1000])
    b = np.array([0, 500, 1000])
    combined = combine_probs(a, b)
    assert (combined >= 0).all() and (combined <= 1).all()
    # combining two nonzero channels should never be less than either alone
    assert combined[1] >= 0.5
    assert combined[2] == 1.0
    assert combined[0] == 0.0


def test_combine_probs_matches_string_formula_two_channels():
    a = np.array([300])
    b = np.array([600])
    expected = 1 - (1 - 0.3) * (1 - 0.6)
    got = combine_probs(a, b)[0]
    assert abs(got - expected) < 1e-9


def test_fisher_z_test_symmetry():
    z1, p1 = fisher_z_test(0.5, 100, 0.3, 100)
    z2, p2 = fisher_z_test(0.3, 100, 0.5, 100)
    assert abs(z1 + z2) < 1e-9
    assert abs(p1 - p2) < 1e-9


def test_fisher_z_test_identical_correlations_gives_zero():
    z, p = fisher_z_test(0.4, 50, 0.4, 80)
    assert abs(z) < 1e-9
    assert abs(p - 1.0) < 1e-9


def test_post_hoc_power_increases_with_n():
    p_small = post_hoc_power(0.3, 30)
    p_large = post_hoc_power(0.3, 300)
    assert p_large > p_small


def test_winsorize_clips_extremes():
    s = pd.Series(np.arange(100, dtype=float))
    w = winsorize(s, lo=0.05, hi=0.95)
    assert w.max() <= s.quantile(0.95) + 1e-9
    assert w.min() >= s.quantile(0.05) - 1e-9
    assert len(w) == len(s)
