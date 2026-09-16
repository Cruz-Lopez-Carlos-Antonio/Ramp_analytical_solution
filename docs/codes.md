---
layout: default
title: Python Codes
math: true
---

## Overview of the Python Scripts

The repository contains the Python 3 scripts developed for the article, which implement the analytical solution for $n(t)$ and the associated delayed neutron precursor concentration $C(t)$, both developed using the Modified Integration Method proposed in the submited paper. Aditionally, a Runge-Kutta of four order, RK4, used as reference solution is included, as well as the computational implementation of the Zhang et al. (2008) and the Palma et al. (2010) solutions. Table A contains the notation used.

**Table A:** Integral definitions and notation in the implementation.

| Integral | Definition | Notation in codes |
| :---: | :---: | :---: |
| $\bar{I}_1(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu} e^{-y^2/2+zy}\,dy$ | `I_1` |
| $\bar{I}_2(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu} e^{-y^2/2-zy}\,dy$ | `I_2` |
| $\bar{I}_3(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu+1} e^{-y^2/2+zy}\,dy$ | `I_3` |
| $\bar{I}_4(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu+1} e^{-y^2/2-zy}\,dy$ | `I_4` |
| $I_5(t),\, I_6(t)$ | $P_{0,1}(t)$ | `I_5, I_6` |

where $\mu=\lambda \beta/a$ and

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0; overflow-x: auto;">

$z(t) \;=\;\sqrt{\frac{a}{\Lambda}}\,t+\frac{\rho_0-\beta}{\sqrt{a\Lambda}}+\lambda\sqrt{\frac{\Lambda}{a}},$
</div>

and:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0; overflow-x: auto;">
$P_k(t)=\displaystyle \int_{0}^{\infty} s^{k}\,
e^{-\frac{1}{a}\!\left(\frac{\Lambda}{2}s^{2}+(\beta-b-a t)\,s\right)}
\, (s+\lambda)^{\mu}\,ds,\quad k\in\{0,1\}$
</div>

---


## 1. Neutron density $n(t)$, using SciPyNumPy (16 digit precision)
<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Ramp_analytical_solution/blob/main/Neutron_density_SciPyNumPy.py" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script implements the analytical solution of the neutron density $n(t)$ written in the following form:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0; overflow-x: auto;">
$$
n(t)=A_1 e^{-\lambda t} \bar{I}_1\big(\mu,z(t)\big)
+ A_2 e^{-\lambda t} \bar{I}_2\big(\mu,z(t)\big)
+ qF\,I_5(t),
$$
</div>

using **SciPy** and **NumPy**. The involved integrals are evaluated using `scipy.integrate.quad`.
The script also contains a linear system, derived from the initial conditions $n(0)$ and $\dot n(0)$, given by:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0; overflow-x: auto;">
$$
\underbrace{\begin{bmatrix}
\bar{I}_1(\mu,z(0)) & \bar{I}_2(\mu,z(0))\\[3pt]
\sqrt{\dfrac{a}{\Lambda}}\,\bar{I}_3(\mu,z(0)) &
-\sqrt{\dfrac{a}{\Lambda}}\,\bar{I}_4(\mu,z(0))
\end{bmatrix}}_{=:\mathbf{L}_1}
\begin{bmatrix} A_1\\ A_2 \end{bmatrix}
=
\underbrace{\begin{bmatrix}
n_0 - q\,F\,I_5(0)\\[3pt]
\dot n_0 + \lambda\big(n_0 - q\,F\,I_5(0)\big) - q\,F\,I_6(0)
\end{bmatrix}}_{=:\mathbf{L}_2}
$$
</div>

i.e.,

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0; overflow-x: auto;">
$$
\mathbf{L}_1
\begin{bmatrix}A_1 \\ A_2\end{bmatrix}
=
\mathbf{L}_2
$$
</div>

to determine the constants $A_1$ and $A_2$ via a least-squares procedure (`numpy.linalg.lstsq`), including column-wise normalization for numerical stability.

**Inputs & Numerical Parameters:**  
The script receives physical, temporal, and numerical control parameters. Specifically, the integration tolerances and the maximum number of subintervals are explicitly defined to guarantee the stability of the adaptive quadrature routine used by `scipy.integrate.quad`. 
*   `epsabs = 0.0`: Forces the algorithm to rely exclusively on relative error control.
*   `epsrel = 1e-9`: Ensures the estimated relative error in each numerical quadrature is kept below $10^{-9}$.

<div style="background:#f4f4f4; border:1px solid #ddd; border-left:4px solid #4a90e2; border-radius:4px; padding:10px; margin-bottom:15px; overflow-x:auto;">
<pre style="margin: 0; background: transparent; border: none; font-family: monospace; color: #333;">
# Physical parameters
gamma_1  = 0.0001    # Slope ramp, a [1/s]  
beta     = 0.0075    # Fraction of precursors [—] 
lambda_1 = 0.001     # Decay constant of precursors [1/s]  
Lambda_1 = 0.0015    # Prompt generation time [s]  
source   = 10**8     # External source intensity, q [n/s]  
rho_s    = -6e-5     # Initial reactivity b = rho(0) [—]  

# Stability and Control epsabs/epsrel
_QKWARGS = dict(epsabs=0.0, epsrel=1e-9, limit=200)

# Evaluation time grid
times = range(0, 21)
</pre>
</div>

**Outputs:**  
The script computes the neutron density $n(t)$ sequentially over the prescribed discrete time grid. These values are printed directly to the standard output and can be easily redirected to a text file for further plotting or analysis.

**Computational Note on Optimization:**  
While this implementation is mathematically accurate and functionally robust, it is **not** computationally optimized. Inside the main `Analytic_n` function, the routine `Constants_Ini_con` is invoked on every single time evaluation. Because the integration constants $A_1$ and $A_2$ depend solely on the initial conditions at $t=0$, recalculating them, normalizing the matrix, and solving the least-squares system at every time step introduces redundant computational overhead. An optimized version would precompute these constants once outside the time loop.

---

## 2. Neutron density $n(t)$, using mpmath (arbitrary precision)
<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Ramp_analytical_solution/blob/main/Neutron_density_mpmath.py" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script provides a **high-precision** version of the analytical solution for $n(t)$, implemented with **mpmath**.  
It includes:

- Multiprecision evaluation of the integrals $I_1,\dots,I_6$ described in the manuscript,
- A robust $2\times 2$ linear solver with row/column scaling and Tikhonov regularization,
- Control of the working precision through `mp.mp.dps`.

This implementation is used as a benchmark to assess conditioning effects and to validate the double-precision results obtained with SciPy/NumPy.

---

## 3. Delayed Precursors Density, $C(t)$, using SciPyNumPy
<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Ramp_analytical_solution/blob/main/C_precursor_SciPyNumPy.py" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script computes the delayed neutron precursor concentration $C(t)$ using the convolution formula

$$
C(t)
= C(0)\,e^{-\lambda t}
+ \frac{\beta}{\Lambda}\,e^{-\lambda t}
\int_0^t e^{\lambda\tau}\,n(\tau)\,d\tau.
$$

The integral is evaluated numerically using `scipy.integrate.cumulative_trapezoid`.  
The function $n(t)$ is imported from `Neutron_density_SciPyNumPy.py`, and the script returns a vectorized approximation of $C(t)$ over a prescribed time grid.

---

## 4. Runge-Kutta, reference solver, using mpmath
<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Ramp_analytical_solution/blob/main/Runge-Kutta%204.py" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script implements a **fourth–order Runge–Kutta (RK4)** solver in 32-digit precision (via mpmath) for the NPKE system.  

It solves simultaneously for $n(t)$ and $C(t)$ using a fine time step, and the resulting numerical solution is used as a high-accuracy reference to validate the analytical formulations and their numerical implementation.

Parameters, time step, and integration interval can be adjusted to reproduce the tables and figures reported in the manuscript.

---

## References

1. Zhang, F., Chen, W.-Z., & Gui, X.-W. (2008). Analytic method study of point-reactor kinetic equation when cold start-up. *Annals of Nuclear Energy, 35*(4), 746–749.
2. Palma, D. A., Martinez, A. S., & Gonçalves, A. C. (2009). Analytical solution of point kinetics equations for linear reactivity variation during the start-up of a nuclear reactor. *Annals of Nuclear Energy, 36*(9), 1469–1471.
