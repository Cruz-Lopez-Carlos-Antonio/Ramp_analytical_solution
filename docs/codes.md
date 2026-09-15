---
layout: default
title: Python Codes
math: true
---

## Overview of the Python Scripts

The repository contains the Python 3 scripts developed in the article, which implement the analytical solution for $n(t)$ and the associated delayed neutron precursor concentration $C(t)$, both developed with the Modified Integration Method proposed in the submited paper. Aditionally, a RK4 reference solution is included, as well as the computational implementation of the Zhang et al. (2008) and the Palma et al. (2010) solutions. Table A contains the notation used.

**Table A:** Integral definitions and notation in the implementation.

| Integral | Definition | Notation in codes |
| :---: | :---: | :---: |
| $\bar{I}_1(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu} e^{-y^2/2+zy}\,dy$ | `I_1` |
| $\bar{I}_2(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu} e^{-y^2/2-zy}\,dy$ | `I_2` |
| $\bar{I}_3(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu+1} e^{-y^2/2+zy}\,dy$ | `I_3` |
| $\bar{I}_4(\mu,z)$ | $\displaystyle \int_{0}^{\infty} y^{\mu+1} e^{-y^2/2-zy}\,dy$ | `I_4` |
| $I_5(t),\, I_6(t)$ | $P_{0,1}$ | `I_5, I_6` |

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


### 1. Neutron density $n(t)$, using SciPyNumPy (16 digit precision)
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

to determine the constants $K_1$ and $K_2$ via a least-squares procedure (`numpy.linalg.lstsq`), including column-wise normalization for numerical stability, and where the terms appearing are defined in the following Table:


### 2. Neutron density $n(t)$, using mpmath (arbitrary precision)
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

### 3. Delayed Precursors Density, $C(t)$, using SciPyNumPy
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

### 4. Runge-Kutta, reference solver, using mpmath
<div style="padding:8px; border-left:4px solid #3c6e71; margin-bottom:10px;">
  <a href="https://github.com/Cruz-Lopez-Carlos-Antonio/Ramp_analytical_solution/blob/main/Runge-Kutta%204.py" 
     target="_blank" style="font-size:16px; color:#22577a; font-weight:bold;">
     👉 Click here to view the code in a new tab
  </a>
</div>

This script implements a **fourth–order Runge–Kutta (RK4)** solver in 32-digit precision (via mpmath) for the NPKE system.  

It solves simultaneously for $n(t)$ and $C(t)$ using a fine time step, and the resulting numerical solution is used as a high-accuracy reference to validate the analytical formulations and their numerical implementation.

Parameters, time step, and integration interval can be adjusted to reproduce the tables and figures reported in the manuscript.
