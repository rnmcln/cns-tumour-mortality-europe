"""Symmetric three-factor decomposition (Das Gupta 1993; equal to the Shapley value over orderings).

Deaths = sum over age-sex cells of size * share * rate. For a change from state 0 to state 1 the
contribution of each factor is (Das Gupta 1993, three multiplicative factors):
  A-effect = dA * [ (B0 C0 + B1 C1) / 3 + (B0 C1 + B1 C0) / 6 ], and symmetrically for B and C,
computed per cell and summed. The three effects add up exactly to the change in deaths.
Arrays may carry a leading replicate axis (bootstrap).
"""
import numpy as np


def _effect(dx, y0, y1, z0, z1):
    return dx * ((y0 * z0 + y1 * z1) / 3 + (y0 * z1 + y1 * z0) / 6)


def decompose(d0, n0, d1, n1, by_cell=False):
    """d*, n*: deaths and person-years by cell (last axis = cells). Returns dict of arrays."""
    d0, n0, d1, n1 = (np.asarray(x, float) for x in (d0, n0, d1, n1))
    s0, s1 = n0.sum(-1, keepdims=True), n1.sum(-1, keepdims=True)
    a0, a1 = s0, s1
    b0, b1 = n0 / s0, n1 / s1
    with np.errstate(invalid='ignore', divide='ignore'):
        c0 = np.where(n0 > 0, d0 / n0, 0.0)
        c1 = np.where(n1 > 0, d1 / n1, 0.0)
    size = _effect(a1 - a0, b0, b1, c0, c1)
    struct = _effect(b1 - b0, a0, a1, c0, c1)
    rate = _effect(c1 - c0, a0, a1, b0, b1)
    if by_cell:                      # cell-level effects (they sum to the totals over the last axis)
        return {'deaths_0': d0, 'deaths_1': d1, 'population_size': size, 'age_sex_structure': struct, 'rates': rate}
    size, struct, rate = size.sum(-1), struct.sum(-1), rate.sum(-1)
    return {'deaths_0': d0.sum(-1), 'deaths_1': d1.sum(-1), 'population_size': size,
            'age_sex_structure': struct, 'rates': rate}
