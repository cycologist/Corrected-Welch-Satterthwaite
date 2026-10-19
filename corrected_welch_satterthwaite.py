"""Corrected Welch-Satterthwaite effective degrees of freedom.

This is a Python translation of corrected_welch_satterthwaite.R.

The corrected estimator is

    nu_vD = (sum(w_k * s2_k))**2 / sum(w_k**2 * s2_k**2 / (nu_k + 2)) - 2

where s2_k are estimated variance components, nu_k are their degrees of
freedom, and w_k are positive weights.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from numbers import Real

__all__ = ["cws", "cws2"]


def _as_real(name: str, value: Real) -> float:
    """Validate one real-valued input and return it as a float."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number")
    value = float(value)
    if math.isnan(value):
        raise ValueError(f"{name} must not be NaN")
    return value


def _as_real_sequence(name: str, values: Sequence[Real]) -> list[float]:
    """Validate a sequence of real-valued inputs."""
    try:
        return [_as_real(name, value) for value in values]
    except TypeError as exc:
        raise TypeError(f"{name} must be a sequence of real numbers") from exc


def cws(
    s2: Sequence[Real],
    nu: Sequence[Real],
    w: Sequence[Real] | None = None,
) -> float:
    """Return the corrected Welch-Satterthwaite degrees of freedom.

    Parameters
    ----------
    s2:
        Estimated variance components. Must all be non-negative.
    nu:
        Component degrees of freedom. Must all be positive.
    w:
        Positive weights. Defaults to equal weights of 1.

    Returns
    -------
    float
        The corrected effective degrees of freedom, or NaN if the
        denominator is zero.

    Raises
    ------
    TypeError
        If an input is not numeric.
    ValueError
        If the inputs have incompatible lengths or violate the
        non-negativity/positivity constraints.
    """
    s2_values = _as_real_sequence("s2", s2)
    nu_values = _as_real_sequence("nu", nu)

    n = len(s2_values)
    if n < 1:
        raise ValueError("s2 must contain at least one element")
    if len(nu_values) != n:
        raise ValueError("s2 and nu must have the same length")

    if w is None:
        w_values = [1.0] * n
    else:
        w_values = _as_real_sequence("w", w)
        if len(w_values) != n:
            raise ValueError("s2 and w must have the same length")

    if any(value < 0 for value in s2_values):
        raise ValueError("all s2 values must be non-negative")
    if any(value <= 0 for value in nu_values):
        raise ValueError("all nu values must be positive")
    if any(value <= 0 for value in w_values):
        raise ValueError("all w values must be positive")

    denominator = math.fsum(
        wi**2 * si**2 / (nui + 2)
        for si, nui, wi in zip(s2_values, nu_values, w_values)
    )
    if denominator <= 0:
        return math.nan

    numerator = math.fsum(
        wi * si for si, wi in zip(s2_values, w_values)
    )
    return numerator**2 / denominator - 2.0


def cws2(
    s2_1: Real,
    s2_2: Real,
    nu_1: Real,
    nu_2: Real,
    w1: Real = 1.0,
    w2: Real = 1.0,
) -> float:
    """Return the corrected Welch-Satterthwaite degrees of freedom for K=2.

    The arguments are scalars; use :func:`cws` for K >= 3 or for
    vectorized input.
    """
    s2_1 = _as_real("s2_1", s2_1)
    s2_2 = _as_real("s2_2", s2_2)
    nu_1 = _as_real("nu_1", nu_1)
    nu_2 = _as_real("nu_2", nu_2)
    w1 = _as_real("w1", w1)
    w2 = _as_real("w2", w2)

    if s2_1 < 0 or s2_2 < 0:
        raise ValueError("all s2 values must be non-negative")
    if nu_1 <= 0 or nu_2 <= 0:
        raise ValueError("all nu values must be positive")
    if w1 <= 0 or w2 <= 0:
        raise ValueError("all w values must be positive")

    s2 = w1 * s2_1 + w2 * s2_2
    denominator = (
        w1**2 * s2_1**2 / (nu_1 + 2)
        + w2**2 * s2_2**2 / (nu_2 + 2)
    )
    if denominator <= 0:
        return math.nan

    return s2**2 / denominator - 2.0
