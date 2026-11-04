"""Helpers to load the analysis panel and derive outcomes and age grids."""
import numpy as np
import pandas as pd

from config import PROCESSED, ESP2013_FULL

OUTCOMES = {  # outcome -> elementary groups summed
    'MAL': ['MAL'],
    'BU': ['D32', 'D33', 'D42', 'D43'],
    'COMP': ['MAL', 'D32', 'D33', 'D42', 'D43'],
    'COMP_NARROW': ['MAL', 'D33', 'D43'],
    'BU_NARROW': ['D33', 'D43'],
    'CALL': ['CALL'],
    'R': ['R'],
    'ALL': ['ALL'],
}


def load_panel(path=None):
    p = pd.read_csv(path or PROCESSED / 'panel_who.csv.gz')
    for o, g in OUTCOMES.items():
        p[o] = p[g].sum(axis=1)
    return p


def collapse(p, grid, by=('geo', 'year', 'sex'), value_cols=None):
    """Collapse age groups (column 'age', lower bounds) to the lower bounds in grid."""
    grid = np.asarray(sorted(grid))
    q = p.copy()
    q['age'] = grid[np.searchsorted(grid, q['age'].to_numpy(), side='right') - 1]
    if value_cols is None:
        value_cols = [c for c in q.columns if c not in set(by) | {'age'} and pd.api.types.is_numeric_dtype(q[c])]
    return q.groupby(list(by) + ['age'], as_index=False)[value_cols].sum()


def esp_weights(grid):
    """2013 ESP weights aggregated to a grid of lower bounds (last group open)."""
    grid = sorted(grid)
    w = {}
    for a, v in ESP2013_FULL.items():
        g = max(x for x in grid if x <= a)
        w[g] = w.get(g, 0) + v
    return pd.Series(w).sort_index()


def asr(df, dcol, pcol, grid):
    """Direct age-standardised rate per 100,000 from a frame with one row per age group."""
    w = esp_weights(grid)
    r = (df.set_index('age')[dcol] / df.set_index('age')[pcol]).reindex(w.index)
    return float((r * w).sum() / w.sum() * 1e5)
