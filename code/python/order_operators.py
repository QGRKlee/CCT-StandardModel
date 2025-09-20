"""Clock operators s (order-5), sigma (order-4), u (order-3) — placeholders.

Each operator will be represented either as an 8x8 matrix in so(8)/Spin(8) representation
or via an action on shell labels with an explicit lift to matrices when required.
"""
from typing import Dict, Any
import numpy as np

def operator_s() -> np.ndarray:
    """Return a representative matrix for the order-5 operator s (placeholder)."""
    raise NotImplementedError

def operator_sigma() -> np.ndarray:
    """Return a representative matrix for the order-4 clicker sigma with sigma^2 = -I (placeholder)."""
    raise NotImplementedError

def operator_u() -> np.ndarray:
    """Return a representative matrix for the order-3 operator u (placeholder)."""
    raise NotImplementedError
