import numpy as np


# =============================================================================
# FUNKCJE METRYK JAKOŚCI
# =============================================================================
def aFactor(y_exp, y_theor):
    """Metryka aFactor – średnia wartość |y - y_est|/(y+y_est)."""
    eps = 1e-12
    return np.mean(np.abs(y_exp - y_theor) / (y_exp + y_theor + eps))

def compute_rmse(y, y_est):
    return np.sqrt(np.mean((y - y_est) ** 2))

def compute_mae(y, y_est):
    return np.mean(np.abs(y - y_est))

def compute_r2(y, y_est):
    ss_res = np.sum((y - y_est) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot != 0 else 0

def compute_pearson(y, y_est):
    if np.std(y) == 0 or np.std(y_est) == 0:
        return 0
    return np.corrcoef(y, y_est)[0, 1]