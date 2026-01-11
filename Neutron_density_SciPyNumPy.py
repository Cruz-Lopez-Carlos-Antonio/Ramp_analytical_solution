# Code for the ramp insertion reactivity using SciPy/NumPy
import math
import numpy as np
from scipy.integrate import quad
import numpy.linalg as npl

# --- Input parameters (Section 5.2, Table 2) --- Based on Zhang et al. (2008, p. 748) and Palma et al. (2009, p. 1471).

gamma_1  = 0.0001    # Slope ramp, gamma [1/s]  
beta     = 0.0075    # Fraction of precursors [—] 
lambda_1 = 0.001     # Decay constant of precursors [1/s]  
Lambda_1 = 0.0015    # Prompt generation time [s]  
source   = 10**8     # Source, q [n/s]  
rho_s    = -6e-5     # b factor in at+b, also rho(0) [—]  

# ============================================================
# ================ z(t, rho_0) - Eq. (48)=====================
# ============================================================

def z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1):
    a1 = math.sqrt(gamma_1 / Lambda_1)
    return a1 * (t + (rho_0 - beta)/gamma_1 + Lambda_1*lambda_1/gamma_1)

# Stability and Control epsabs/epsrel (Section 5.3)

_QKWARGS = dict(epsabs=0.0, epsrel=1e-9, limit=200)

# ============================================================
# === Integrals - Tables 3, Eq. (49) ===
# ============================================================

def I_1(mu, zval):
    def integrand(x):
        return x**mu * np.exp(-x**2/2 + x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_2(mu, zval):
    def integrand(x):
        return x**mu * np.exp(-x**2/2 - x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_3(mu, zval):
    def integrand(x):
        return x**(mu+1) * np.exp(-x**2/2 + x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_4(mu, zval):
    def integrand(x):
        return x**(mu+1) * np.exp(-x**2/2 - x*zval)
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_5(t, rho_0, lambda_1, beta, gamma_1, Lambda_1):
    mu = lambda_1*beta/gamma_1
    A = Lambda_1/(2*gamma_1)
    B = (rho_0 - beta)/gamma_1 + t
    def integrand(x):
        return (x + lambda_1)**mu * np.exp(-(A*x**2 - B*x))
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

def I_6(t, rho_0, lambda_1, beta, gamma_1, Lambda_1):
    mu = lambda_1*beta/gamma_1
    A = Lambda_1/(2*gamma_1)
    B = (rho_0 - beta)/gamma_1 + t
    def integrand(x):
        return x*(x + lambda_1)**mu * np.exp(-(A*x**2 - B*x))
    val, _ = quad(integrand, 0.0, np.inf, **_QKWARGS)
    return float(val)

# ============================================================
# === Prefactor F ====== Table 2, 
# ============================================================

def prefactor_F(lambda_1, beta, gamma_1, Lambda_1):
    mu = lambda_1*beta/gamma_1
    return (Lambda_1/gamma_1) * (lambda_1**(-mu))

# ============================================================
# === Initial conditions - Eq. (50) ===
# ============================================================

def Constants_Ini_con(lambda_1, beta, Lambda_1, gamma_1, rho_0, n_0, dn_0, q):
    # Linear System given in Eqs. (50-54)
    zeta = z(0.0, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    mu   = lambda_1*beta/gamma_1
    Int1 = I_1(mu, zeta)
    Int2 = I_2(mu, zeta)
    Int3 = I_3(mu, zeta)
    Int4 = I_4(mu, zeta)
    Int5 = I_5(0.0, rho_0, lambda_1, beta, gamma_1, Lambda_1)
    Int6 = I_6(0.0, rho_0, lambda_1, beta, gamma_1, Lambda_1)
    F = prefactor_F(lambda_1, beta, gamma_1, Lambda_1)
    l1 = np.array([[Int1,                                  Int2],
        [math.sqrt(gamma_1/Lambda_1)*Int3, -math.sqrt(gamma_1/Lambda_1)*Int4]
    ], dtype=float)
    rhs1 = n_0 - q * F * Int5
    rhs2 = dn_0 + lambda_1 * (n_0 - q * F * Int5) - q * F * Int6
    l2 = np.array([rhs1, rhs2], dtype=float)
    #================ Normalization Proces =================================
    col_norms = npl.norm(l1, axis=0)
    l1_scaled = l1 / col_norms
    #====================== Eq. (52) ======================================= 
    #==================== Least Square Method ==============================
    scale_rhs = np.max(np.abs(l2)) if np.max(np.abs(l2)) != 0 else 1.0
    l2_scaled = l2 / scale_rhs
    A_scaled, *_ = npl.lstsq(l1_scaled, l2_scaled, rcond=None)
    #==================== Back-scaling Eq. (53) ============================
    A1, A2 = (A_scaled / col_norms) * scale_rhs
    return float(A1), float(A2)

# ============================================================
# === Analytical solution n(t) - Eq. (49) ===
# ============================================================

def Analytic_n(t, rho_0, beta, Lambda_1, gamma_1, lambda_1, n_0, dn_0, q):
    # Eq. (49) 
    A1, A2 = Constants_Ini_con(lambda_1, beta, Lambda_1, gamma_1, rho_0, n_0, dn_0, q)
    mu   = lambda_1*beta/gamma_1
    zeta = z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    F    = prefactor_F(lambda_1, beta, gamma_1, Lambda_1)

    first  = A1 * np.exp(-lambda_1 * t) * I_1(mu, zeta)
    second = A2 * np.exp(-lambda_1 * t) * I_2(mu, zeta)
    third  = q  * F * I_5(t, rho_0, lambda_1, beta, gamma_1, Lambda_1)

    sol = first + second + third
    return float(sol)

# ============================================================
# === Execution ===
# ============================================================

if __name__ == "__main__":
    n0 = source * Lambda_1 / abs(rho_s)
    times = range(0, 21)
    for k in times:
        val = Analytic_n(k, rho_s, beta, Lambda_1, gamma_1, lambda_1, n0, 0.0, source)
        print(val)
