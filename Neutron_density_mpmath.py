# Code for the ramp insertion reactivity using mpmath
import mpmath as mp
# -------------------- Precision --------------------
mp.mp.dps = 32
mp.mp.rounding = 'nearest'

# ----------- Parameters Zhang et al. (2008, p. 748) and Palma et al. (2009, p. 1471) by default --------------------
gamma_1  = mp.mpf('0.0001')     # Slope ramp, gamma [1/s]
beta     = mp.mpf('0.0075')     # Fraction of precursors [—]
lambda_1 = mp.mpf('0.001')      # Decay constant of precursors [1/s]
Lambda_1 = mp.mpf('0.0015')     # Prompt generation time [s]
source   = mp.mpf('1.0e8')      # Source, q [n/s]
rho_s    = mp.mpf('-6.0e-5')    # b factor in at+b, also ρ(0) [—]

# -------------------- z(t, rho_0), Eq. (48)--------------------
def z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1):
    a1 = mp.sqrt(gamma_1 / Lambda_1)
    return a1 * (mp.mpf(t) + (rho_0 - beta)/gamma_1 + Lambda_1*lambda_1/gamma_1)

# -------------------- Integrals Table 3 --------------------
def I_1(mu, zval):
    return mp.quad(lambda x: mp.power(x, mu) * mp.e**(-x**2/2 + x*zval), [0, mp.inf])

def I_2(mu, zval):
    return mp.quad(lambda x: mp.power(x, mu) * mp.e**(-x**2/2 - x*zval), [0, mp.inf])

def I_3(mu, zval):
    return mp.quad(lambda x: mp.power(x, mu+1) * mp.e**(-x**2/2 + x*zval), [0, mp.inf])

def I_4(mu, zval):
    return mp.quad(lambda x: mp.power(x, mu+1) * mp.e**(-x**2/2 - x*zval), [0, mp.inf])

def I_5(t, rho_0, lambda_1, beta, gamma_1, Lambda_1):
    mu = (lambda_1*beta)/gamma_1
    a = Lambda_1/(2*gamma_1)
    b = ((rho_0 - beta)/gamma_1 + mp.mpf(t))
    return mp.quad(lambda x: mp.power(x + lambda_1, mu) * mp.e**(-(a*x**2 - b*x)),\
                   [0, mp.inf])

def I_6(t, rho_0, lambda_1, beta, gamma_1, Lambda_1):
    mu = (lambda_1*beta)/gamma_1
    a = Lambda_1/(2*gamma_1)
    b = ((rho_0 - beta)/gamma_1 + mp.mpf(t))
    return mp.quad(lambda x: x*mp.power(x + lambda_1, mu) * mp.e**(-(a*x**2 - b*x)),\
                   [0, mp.inf])

# -------------------- Prefactor F given in Table 2 --------------------
def prefactor_F(lambda_1, beta, gamma_1, Lambda_1):
    mu = (lambda_1*beta)/gamma_1
    return (Lambda_1/gamma_1) * mp.power(lambda_1, -mu)

# ==================== Robust Solver 2x2, Algorithm 4 ====================
def _row_col_scale(A, b):
    #Algorithm 2
    A = mp.matrix(A)  
    b = mp.matrix(b)

    # Step 4 in Algorithm 2
    r1 = 1 / mp.nsum(lambda j: mp.fabs(A[0,j]), [0,1])\
         if mp.nsum(lambda j: mp.fabs(A[0,j]), [0,1]) != 0 else 1
    r2 = 1 / mp.nsum(lambda j: mp.fabs(A[1,j]), [0,1])\
         if mp.nsum(lambda j: mp.fabs(A[1,j]), [0,1]) != 0 else 1
    # Steps 5 in Algorithm 2
    R = mp.matrix([[r1,0],[0,r2]])
    A = R * A
    b = R * b

    #Step 6 in Algorithm 2
    c1 = 1 / (mp.fabs(A[0,0]) + mp.fabs(A[1,0])) if \
         (mp.fabs(A[0,0]) + mp.fabs(A[1,0])) != 0 else 1
    c2 = 1 / (mp.fabs(A[0,1]) + mp.fabs(A[1,1])) if \
         (mp.fabs(A[0,1]) + mp.fabs(A[1,1])) != 0 else 1

    #Step 7 in Algorithm 2
    C = mp.matrix([[c1,0],[0,c2]])
    A = A * C

    #Step 8
    def undo_col_scale(x):
        # x_orig = C * x_scaled
        return C * x
    
    #Step 9
    return A, b, undo_col_scale

def _tikhonov(A, b, tau=mp.mpf('1e-28')):#Algorithm 3
    #Step 4 in Algorithm 3
    AT = A.T
    ATA = AT * A
    ATb = AT * b
    #Step 5 in Algorithm 3
    diag_max = mp.nmax([mp.fabs(ATA[0,0]), mp.fabs(ATA[1,1])])
    all_max  = mp.nmax([mp.fabs(ATA[0,0]), mp.fabs(ATA[0,1]), \
                        mp.fabs(ATA[1,0]), mp.fabs(ATA[1,1])])
    #Step 6 in Algorithm 3
    eps = tau * (diag_max + all_max if (diag_max + all_max) != 0 else 1)
    #Step 7 in Algorithm 3
    ATAreg = mp.matrix([[ATA[0,0]+eps, ATA[0,1]],
                        [ATA[1,0],     ATA[1,1]+eps]])
    #Step 8 in Algorithm 3
    try:
        x = mp.lu_solve(ATAreg, ATb)
    #Step 9 in Algorithm 3
    except ZeroDivisionError:
        with mp.workdps(mp.mp.dps + 30):
            x = mp.lu_solve(ATAreg, ATb)
    return x

def solve2x2_robust(A, b):
    #Step 3 of Algorithm 4
    A = mp.matrix(A); b = mp.matrix(b)
    #Attempt 1 in Step 4 of Algorithm 4
    try:
        return mp.lu_solve(A, b)
    except ZeroDivisionError:
        pass

    # Attempt 2 in Step 7 of Algorithm 4
    for extra in (20, 40):
        try:
            with mp.workdps(mp.mp.dps + extra):
                return mp.lu_solve(A, b)
        except ZeroDivisionError:
            continue

    # Attempt 3 in Step 14 of Algorithm 4
    As, bs, undo = _row_col_scale(A, b)
    try:
        xs = mp.lu_solve(As, bs)
        return undo(xs)
    except ZeroDivisionError:
        pass

    # Attempt 4 Tikhonov routine 
    xs = _tikhonov(As, bs, tau=mp.mpf('1e-28'))
    return undo(xs)

# -------------------- Initial conditions, Eq. (50) --------------------
def Constants_Ini_con(lambda_1, beta, Lambda_1, gamma_1, rho_0, \
                      n_0, dn_0, q):
    zeta = z(0, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    mu   = (lambda_1*beta)/gamma_1

    Int1 = I_1(mu, zeta)
    Int2 = I_2(mu, zeta)
    Int3 = I_3(mu, zeta)
    Int4 = I_4(mu, zeta)
    Int5 = I_5(0, rho_0, lambda_1, beta, gamma_1, Lambda_1)
    Int6 = I_6(0, rho_0, lambda_1, beta, gamma_1, Lambda_1)

    F = prefactor_F(lambda_1, beta, gamma_1, Lambda_1)

    a10 = Int1
    a11 = Int2
    a20 = mp.sqrt(gamma_1/Lambda_1) * Int3
    a21 = - mp.sqrt(gamma_1/Lambda_1) * Int4

    A = mp.matrix([[a10, a11],
                   [a20, a21]])

    rhs0 = n_0 - q * F * Int5
    rhs1 = dn_0 + lambda_1 * (n_0 - q * F * Int5) - q * F * Int6
    bvec = mp.matrix([rhs0, rhs1])

    # --- Robust solver ---
    sol = solve2x2_robust(A, bvec)
    A1, A2 = sol[0], sol[1]
    return A1, A2

# -------------------- Analytical solution --------------------
def Analytic_n(t, rho_0, beta, Lambda_1, gamma_1, lambda_1, n_0, dn_0, q):
    A1, A2 = Constants_Ini_con(lambda_1, beta, Lambda_1, gamma_1, \
                               rho_0, n_0, dn_0, q)
    mu   = (lambda_1*beta)/gamma_1
    zeta = z(t, rho_0, beta, Lambda_1, gamma_1, lambda_1)
    F    = prefactor_F(lambda_1, beta, gamma_1, Lambda_1)

    First_term  = A1 * mp.e**(-lambda_1 * mp.mpf(t)) * I_1(mu, zeta)
    Second_term = A2 * mp.e**(-lambda_1 * mp.mpf(t)) * I_2(mu, zeta)
    Third_term  = q  * F * I_5(t, rho_0, lambda_1, beta, gamma_1, Lambda_1)

    Solution = First_term + Second_term + Third_term
    print(t, Solution)
    return Solution

# -------------------- Execution --------------------
if __name__ == "__main__":
    n0 = source * Lambda_1 / mp.fabs(rho_s)
    for k in range(0, 21):
        Analytic_n(k, rho_s, beta, Lambda_1, gamma_1, lambda_1, n0, \
                   mp.mpf('0.0'), source)
