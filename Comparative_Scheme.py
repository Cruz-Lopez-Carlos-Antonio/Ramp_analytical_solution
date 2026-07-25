import math
import time
import numpy as np
from scipy.integrate import quad
import numpy.linalg as npl
from scipy.special import gamma, gammaincc  # Palma

# ------------------------------------------------------------
# Physical parameters (common to all methods)
# ------------------------------------------------------------
a       = 1.0e-4     # ramp slope, rho(t) = a t + b
beta    = 0.0075     # delayed neutron fraction
lam     = 0.001      # precursor decay constant
Lambda1 = 0.0015     # neutron generation time
q       = 1.0e8      # external source
b       = -6.0e-5    # constant term of the ramp

# Initial conditions for Zhang et al. scenario
n0  = q * Lambda1 / abs(b)
dn0 = 0.0

# Internal notation for the analytical solution (Proposed solution)
gamma_1  = a
lambda_1 = lam
Lambda_1 = Lambda1
rho_0    = b

# Quadrature settings for the analytical integrals
_QKWARGS = dict(epsabs=0.0, epsrel=1e-9, limit=200)

# ============================================================
# 1) MIM analytical solution
# ============================================================

def z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1):
    a1 = math.sqrt(gamma_1 / Lambda_1)
    return a1 * (t + (rho_0 - beta)/gamma_1 + Lambda_1*lambda_1/gamma_1)

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

def prefactor_F(lambda_1, beta, gamma_1, Lambda_1):
    mu = lambda_1*beta/gamma_1
    return (Lambda_1/gamma_1) * (lambda_1**(-mu))

def prepare_analytic_solution(lambda_1, beta, Lambda_1, gamma_1, rho_0, n_0, dn_0, q):
    """
    Computes all time-independent quantities required by the MIM
    analytical solution. This preprocessing stage is performed once
    before evaluating n(t) over the requested time grid.
    """
    mu = lambda_1 * beta / gamma_1
    F = prefactor_F(lambda_1, beta, gamma_1, Lambda_1)

    # Auxiliary variable and analytical integrals evaluated at t = 0
    zeta_0 = z(0.0, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    Int1 = I_1(mu, zeta_0)
    Int2 = I_2(mu, zeta_0)
    Int3 = I_3(mu, zeta_0)
    Int4 = I_4(mu, zeta_0)
    Int5 = I_5(0.0, rho_0, lambda_1, beta, gamma_1, Lambda_1)
    Int6 = I_6(0.0, rho_0, lambda_1, beta, gamma_1, Lambda_1)

    # Linear system for the integration constants A1 and A2
    sqrt_factor = math.sqrt(gamma_1 / Lambda_1)
    l1 = np.array([
        [Int1,                    Int2],
        [sqrt_factor * Int3, -sqrt_factor * Int4]
    ], dtype=float)

    rhs1 = n_0 - q * F * Int5
    rhs2 = dn_0 + lambda_1 * (n_0 - q * F * Int5) - q * F * Int6
    l2 = np.array([rhs1, rhs2], dtype=float)

    # Column-wise normalization of the coefficient matrix
    col_norms = npl.norm(l1, axis=0)
    if np.any(col_norms == 0.0):
        raise np.linalg.LinAlgError(
            "At least one column of the initial-condition matrix has zero norm."
        )
    l1_scaled = l1 / col_norms

    # Normalization of the right-hand side
    max_rhs = np.max(np.abs(l2))
    if max_rhs == 0.0:
        max_rhs = 1.0
    l2_scaled = l2 / max_rhs

    # Least-squares solution of the normalized system
    A_scaled, *_ = npl.lstsq(l1_scaled, l2_scaled, rcond=None)

    # Back-scaling of the integration constants
    A1, A2 = (A_scaled / col_norms) * max_rhs

    return {
        "A1": float(A1),
        "A2": float(A2),
        "mu": float(mu),
        "F": float(F),
    }


def Analytic_n(t, rho_0, beta, Lambda_1, gamma_1, lambda_1, q, precomputed):
    """
    Evaluates the optimized MIM analytical solution using quantities
    previously computed by prepare_analytic_solution().
    """
    A1 = precomputed["A1"]
    A2 = precomputed["A2"]
    mu = precomputed["mu"]
    F = precomputed["F"]

    zeta = z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    decay_factor = math.exp(-lambda_1 * t)

    first = A1 * decay_factor * I_1(mu, zeta)
    second = A2 * decay_factor * I_2(mu, zeta)
    third = q * F * I_5(t, rho_0, lambda_1, beta, gamma_1, Lambda_1)

    return float(first + second + third)

# ============================================================
# 2) Zhang et al. analytical approximation
# ============================================================
def n_zhang(t, a, beta, Lambda1, q, b, n0):
    numerator   = beta * n0 + q * Lambda1
    denominator = beta - (a * t + b)
    return numerator / denominator

# ============================================================
# 3) Palma et al. analytical approximation
# ============================================================

# Notation mapping:
q0  = b
r   = a
k   = lam
ell = Lambda1

k1 = (k * q0 + r) / r
k2 = (beta - q0) / r      # = A3
k3 = k * q * ell / r
s  = 1.0 + k * beta / r   # Exponent appearing in Palma's analytical solution

A3 = k2

def upper_incomplete_gamma(s, z):
    # Γ(s, z) = Γ(s) * gammaincc(s, z)
    return gamma(s) * gammaincc(s, z)

A2 = -upper_incomplete_gamma(s, k * A3) \
     + (r * n0 * (k * A3)**s) / (k * q * ell * np.exp(k * A3))

# Global normalization constant determined from n(0) = n0
B_palma = n0 * (A3 ** s) / (upper_incomplete_gamma(s, k * A3) + A2)

def n_palma(t):
    t = float(t)
    x = A3 - t
    val = B_palma * math.exp(-k * t) / (x**s) * (upper_incomplete_gamma(s, k * x) + A2)
    return float(val)
    
# ============================================================
# 4) Fourth-order Runge–Kutta method
#    (one delayed-neutron precursor group and ramp reactivity)
# ============================================================

def rho_t(t, a, b):
    return a * t + b

def rhs_point_kinetics(t, y, a, b, beta, lam, Lambda1, q):
    """
    Right-hand side of the Neutron Point Kinetics Equations.

    y = [n, C]
    """
    n, C = y
    rho = rho_t(t, a, b)
    dn = (rho - beta)/Lambda1 * n + lam * C + q
    dC = beta/Lambda1 * n - lam * C
    return np.array([dn, dC], dtype=float)

def rk4_step(f, t, y, h, args):
    k1 = f(t,          y,              *args)
    k2 = f(t + 0.5*h,  y + 0.5*h*k1,   *args)
    k3 = f(t + 0.5*h,  y + 0.5*h*k2,   *args)
    k4 = f(t + h,      y + h*k3,       *args)
    return y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def rk4_solve_and_sample(times, a, b, beta, lam, Lambda1, q, n0):
    """
    Integrates the system once from t = 0 to t_max using the
    fourth-order Runge–Kutta method with a fixed step size h,
    and then returns n(t) at the requested sampling times.
    """
    times = np.asarray(times, dtype=float)
    t_max = float(times.max())

    # Integration step size; adjust it according to the desired accuracy.
    h = 1.0e-3  # For example, 0.001 s gives 20,000 steps up to t = 20 s.
    N = int(round(t_max / h))

    # Initial conditions
    C0 = beta * n0 / (lam * Lambda1)
    y = np.array([n0, C0], dtype=float)

    # Arrays used to store the complete solution for subsequent sampling
    n_full = np.empty(N+1, dtype=float)
    t_full = np.empty(N+1, dtype=float)

    n_full[0] = n0
    t_full[0] = 0.0

    args = (a, b, beta, lam, Lambda1, q)

    t = 0.0
    for k in range(1, N+1):
        y = rk4_step(rhs_point_kinetics, t, y, h, args)
        t = k * h
        n_full[k] = y[0]
        t_full[k] = t

    # Sample n(t) at the requested times, assumed to be multiples of h.
    n_out = np.empty_like(times)
    inv_h = 1.0/h
    for i, tt in enumerate(times):
        idx = int(round(tt * inv_h))
        n_out[i] = n_full[idx]

    return n_out
# ============================================================
# Generic averaged benchmark routine for scalar n(t) methods
# ============================================================

def benchmark_method(method, times, args=(), n_runs=5, warmup=True):
    """
    Benchmarks a scalar method of the form method(t, *args),
    evaluated independently at each requested time.
    """
    times = np.asarray(times, dtype=float)

    # Optional warm-up run, excluded from the timing measurements
    if warmup:
        _ = [method(t, *args) for t in times]

    elapsed_list = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        _ = [method(t, *args) for t in times]
        t1 = time.perf_counter()
        elapsed_list.append(t1 - t0)

    elapsed_array = np.array(elapsed_list)
    mean_elapsed = elapsed_array.mean()
    std_elapsed  = elapsed_array.std(ddof=1) if n_runs > 1 else 0.0
    time_per_eval = mean_elapsed / len(times)

    return mean_elapsed, std_elapsed, time_per_eval


def benchmark_rk4(times, a, b, beta, lam, Lambda1, q, n0, n_runs=5, warmup=True):
    """
    RK4-specific benchmark: integrates the complete time interval once
    and samples the solution at the requested times.
    """
    times = np.asarray(times, dtype=float)

    # Optional warm-up integration, excluded from the timing measurements
    if warmup:
        _ = rk4_solve_and_sample(times, a, b, beta, lam, Lambda1, q, n0)

    elapsed_list = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        _ = rk4_solve_and_sample(times, a, b, beta, lam, Lambda1, q, n0)
        t1 = time.perf_counter()
        elapsed_list.append(t1 - t0)

    elapsed_array = np.array(elapsed_list)
    mean_elapsed = elapsed_array.mean()
    std_elapsed  = elapsed_array.std(ddof=1) if n_runs > 1 else 0.0

    # Average total integration time divided by the number of sampled points
    time_per_eval = mean_elapsed / len(times)

    return mean_elapsed, std_elapsed, time_per_eval

# ============================================================
# Main program: comparison at t = 0, 1, ..., 20
# ============================================================

if __name__ == "__main__":

    times = np.arange(0, 21, 1, dtype=float)

    # Compute the time-independent MIM quantities only once.
    proposed_precomputed = prepare_analytic_solution(
        lambda_1,
        beta,
        Lambda_1,
        gamma_1,
        rho_0,
        n0,
        dn0,
        q,
    )

    analytic_args = (
        rho_0,
        beta,
        Lambda_1,
        gamma_1,
        lambda_1,
        q,
        proposed_precomputed,
    )
    zhang_args = (a, beta, Lambda1, q, b, n0)

    n_runs = 5

    print("=== Benchmark at t = 0, 1, 2, ..., 20 ===")
    print(
        f"Parameters: a={a}, beta={beta}, lambda={lam}, "
        f"Lambda={Lambda1}, q={q}, b={b}"
    )
    print(f"Initial n0 = {n0:.6e}")
    print()

    # --- Proposed solution ---
    mean_A, std_A, per_eval_A = benchmark_method(
        Analytic_n,
        times,
        args=analytic_args,
        n_runs=n_runs,
        warmup=True,
    )
    print("Analytical solution (Eq. (69), Proposed solution):")
    print(f"  Mean total time: {mean_A:.6f} s ± {std_A:.6f} s")
    print(f"  Mean time per evaluation: {per_eval_A:.6e} s\n")

    # --- Zhang et al. approximation ---
    mean_Z, std_Z, per_eval_Z = benchmark_method(
        n_zhang,
        times,
        args=zhang_args,
        n_runs=n_runs,
        warmup=True,
    )
    print("Zhang et al. approximate solution:")
    print(f"  Mean total time: {mean_Z:.6f} s ± {std_Z:.6f} s")
    print(f"  Mean time per evaluation: {per_eval_Z:.6e} s\n")

    # --- Palma et al. analytical approximation ---
    mean_P, std_P, per_eval_P = benchmark_method(
        n_palma,
        times,
        args=(),
        n_runs=n_runs,
        warmup=True,
    )
    print("Palma et al. analytical solution:")
    print(f"  Mean total time: {mean_P:.6f} s ± {std_P:.6f} s")
    print(f"  Mean time per evaluation: {per_eval_P:.6e} s\n")

    # --- Fourth-order Runge-Kutta method ---
    mean_R, std_R, per_eval_R = benchmark_rk4(
        times,
        a,
        b,
        beta,
        lam,
        Lambda1,
        q,
        n0,
        n_runs=n_runs,
        warmup=True,
    )
    print("Fourth-order Runge-Kutta numerical solution:")
    print(f"  Mean total time: {mean_R:.6f} s ± {std_R:.6f} s")
    print(
        "  Mean time per evaluation (over these 21 sampled points): "
        f"{per_eval_R:.6e} s\n"
    )

    # Optional comparison at selected times
    sample_times = [0, 5, 10, 15, 20, 30, 40, 60, 62, 64, 66, 68, 70, 75, 80]
    print("Comparison of n(t) at t = 0, 5, 10, 15, 20, 30, 40, 60, 62, 64, 66, 68, 70, 75, 80 :")

    nA_list = [Analytic_n(t, *analytic_args) for t in sample_times]
    nZ_list = [n_zhang(t, *zhang_args) for t in sample_times]
    nP_list = [n_palma(t) for t in sample_times]
    nR_full = rk4_solve_and_sample(
        sample_times, a, b, beta, lam, Lambda1, q, n0
    )

    for t, nA, nZ, nP, nR in zip(
        sample_times, nA_list, nZ_list, nP_list, nR_full
    ):
        # 1. Imprime la salida original que ya tenías
        print(
            f"  t = {t:4.1f} s -> Proposed solution = {nA:.6e}, "
            f"Zhang = {nZ:.6e}, Palma = {nP:.6e}, RK4 = {nR:.6e}"
        )
        
        # 2. Calcula el APE para cada método (protegiendo contra división por cero)
        ape_A = abs((nA - nR) / nR) * 100.0 if nR != 0 else 0.0
        ape_Z = abs((nZ - nR) / nR) * 100.0 if nR != 0 else 0.0
        ape_P = abs((nP - nR) / nR) * 100.0 if nR != 0 else 0.0
        
        # 3. Imprime el APE en una nueva línea indentada
        print(
            f"               -> APE (%): Proposed solution = {ape_A:.3f}%, "
            f"Zhang = {ape_Z:.3f}%, Palma = {ape_P:.3f}%\n"
        )

