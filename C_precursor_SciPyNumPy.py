import math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from Neutron_density_SciPyNumPy import Analytic_n  

# --- Invoking the n(t) function given in ---
#     Neutron_density_SciPyNumPy.py

def make_n_func(rho_s, beta, Lambda_1, gamma_1, lambda_1, n0, dn0, source):
    def n_func(t):
        return float(Analytic_n(t, rho_s, beta, Lambda_1, gamma_1, lambda_1, n0, dn0, source))
    return n_func

# --- C(t) vectorizing---
def C_vector(ts, n_func, beta, Lambda_1, lambda_1, C0):
    ts = np.asarray(ts, dtype=float)
    f = np.array([math.exp(lambda_1*tau) * n_func(tau) for tau in ts], dtype=float)
    F = cumulative_trapezoid(f, ts, initial=0.0)  # ∫_0^t e^{λ τ} n(τ) dτ
    return np.exp(-lambda_1*ts) * (C0 + (beta/Lambda_1)*F)

if __name__ == "__main__":
    # --- Nuclear parameters 
    gamma_1  = 0.0001    # Slope ramp, gamma [1/s]  
    beta     = 0.0075    # Fraction of precursors [—] 
    lambda_1 = 0.001     # Decay constant of precursors [1/s]  
    Lambda_1 = 0.0015    # Prompt generation time [s]  
    source   = 10**8     # Source, q [n/s]  
    rho_s    = -6e-5     # b factor in at+b, also rho(0) [—] 
    # --- Initial conditions ---
    n0   = source * Lambda_1 / abs(rho_s)
    C0   = (beta / (lambda_1 * Lambda_1)) * n0
    dn0  = 0.0
    # --- Invoking n(t) ---
    n_func = make_n_func(rho_s, beta, Lambda_1, gamma_1, lambda_1, n0, dn0, source)

    # --- Mesh for the numerical integration ---
    ts = np.linspace(0.0, 1.0, 5001)  # 0 to 1 s

    # --- Evaluation ---
    Cs = C_vector(ts, n_func, beta, Lambda_1, lambda_1, C0)
    for t, c in zip(ts, Cs):
        if abs(t - round(t)) < 1e-12:  
            print(f"{t:.0f} {c:.12e}")
