import mpmath as mp

# -------------------- Precisión --------------------
mp.mp.dps = 32
mp.mp.rounding = 'nearest'

# ----------- Parameters (Zhang et al., 2008)) by default --------------------
a           = mp.mpf('0.0005')     # Slope ramp, gamma [1/s]
q           = mp.mpf('1.0e8')      # Source, q [n/s]
b           = mp.mpf('-6.0e-5')    # b factor in at+b, also ρ(0) [—]
beta        = mp.mpf('0.0075')     # Fraction of precursors [—] 
lambda_1    = mp.mpf('0.001')      # Decay constant of precursors
Lambda_1    = mp.mpf('0.0015')      # Prompt generation time [s] 

# -------------------- Initial Conditions --------------------
n0 = q * Lambda_1 / mp.fabs(b)
C0 = (beta/(Lambda_1*lambda_1)) * n0

# -------------------- Temporal Domain --------------------
t0 = mp.mpf('0.0')
tf = mp.mpf('20.0')
N  = 200_000           # Steps

# -------------------- NPKE and linear reactivity --------------------
def derivs(t, y): #System given in Eqs. (13-14)
    n, C = y
    dn = ((rho(t) - beta)/Lambda_1)*n + lambda_1*C + q
    dC = (beta/Lambda_1)*n - lambda_1*C
    return (dn, dC)

def rho(t): #Linear reactivity
    return a*t + b
# -------------------- Integrador RK4 --------------------
def rk4_step(f, t, y, h):
    n, C = y
    k1n, k1c = f(t, (n, C))
    k2n, k2c = f(t + h/2, (n + h*k1n/2, C + h*k1c/2))
    k3n, k3c = f(t + h/2, (n + h*k2n/2, C + h*k2c/2))
    k4n, k4c = f(t + h,   (n + h*k3n,   C + h*k3c))
    n_new = n + (h/6)*(k1n + 2*k2n + 2*k3n + k4n)
    c_new = C + (h/6)*(k1c + 2*k2c + 2*k3c + k4c)
    return (n_new, c_new)

def solve_rk4(f, t0, tf, y0, N):
    t0 = mp.mpf(t0); tf = mp.mpf(tf); N = int(N)
    h = (tf - t0)/N
    T = [t0 + i*h for i in range(N+1)]
    Nvals = [mp.mpf('0.0')]*(N+1)
    Cvals = [mp.mpf('0.0')]*(N+1)
    Nvals[0], Cvals[0] = y0
    y = y0
    t = t0
    for i in range(N):
        y = rk4_step(f, t, y, h)
        t = T[i+1]
        Nvals[i+1], Cvals[i+1] = y
    return T, Nvals, Cvals

# -------------------- Execution --------------------
if __name__ == "__main__":
    y0 = (n0, C0)
    T, Nvals, Cvals = solve_rk4(derivs, t0, tf, y0, N)

    # Vector of times of interest
    times = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

    h = (tf - t0)/N  # Temporal step
    print("dps =", mp.mp.dps)
    print("\n# t [s] \t n(t) \t C(t)")

    for t_req in times:
        idx = int(round((t_req - t0)/h))
        if 0 <= idx < len(T):
            print(f"{T[idx]}\t{Nvals[idx]}\t{Cvals[idx]}")
