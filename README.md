# Analytical Solution for the Neutron Point Kinetics Equations with Ramp Reactivity

The present repository contains the Python 3 codes associated with the development of a new analytical solution for the Neutron Point Kinetics Equations (NPKE), considering a single group of delayed neutron precursors and a linear ramp reactivity of the form  
\[
\rho(t)=a t + b.
\]

These codes accompany the manuscript *Analytical solutions of the Neutron Point Kinetics Equations under a linear reactivity ramp*, recently submitted to the **Computer Physics Communications** journal.

Unless otherwise noted, the programs are available under an MIT License.

**Authors:**  
Carlos-Antonio Cruz-López (cacl.nucl@gmail.com)  
Gilberto Espinosa-Paredes  
Juan-Luis François

---

## Overview of the Repository

The repository contains four Python scripts.  
Each one implements a different component of the analytical or numerical framework used in the manuscript.

### **1. `Neutron_density_SciPyNumPy.py`**

This script implements the analytical solution of the neutron density \( n(t) \) using SciPy and NumPy, based on the integral representation:

\[
n(t)=K_1 e^{-\lambda t}\!\int_0^\infty e^{-y^2/2 + \mathcal{E}_3(t) y}y^{\lambda\beta/a}dy
    +K_2 e^{-\lambda t}\!\int_0^\infty e^{-y^2/2 - \mathcal{E}_3(t) y}y^{\lambda\beta/a}dy
    +B_{p,1}\!\int_0^\infty e^{-\frac{1}{a}\left(\frac{\Lambda}{2}s^2+(\beta-b-at)s\right)}(s+\lambda)^{\lambda\beta/a}ds.
\]

The code evaluates the integrals numerically using `scipy.integrate.quad` and includes a numerical procedure for determining the constants \(K_1\) and \(K_2\) through a least-squares solver.

---

### **2. `Neutron_density_mpmath.py`**

This version implements the same analytical solution for \(n(t)\), but using high-precision arithmetic through **mpmath** (`mp.mp.dps = 32` by default).  
It includes:

- multiprecision versions of the integrals \(I_1\)–\(I_6\),  
- the same analytical structure for \(n(t)\),  
- and a **robust 2×2 solver** (row-column scaling + Tikhonov fallback) to guarantee numerical stability in the determination of \(K_1\) and \(K_2\).

This script is used as a benchmark for verifying numerical sensitivity and conditioning issues.

---

### **3. `C_precursor_SciPyNumPy.py`**

This program computes the delayed neutron precursor concentration \(C(t)\) using the convolution formula:

\[
C(t)
= C(0) e^{-\lambda t}
+ \frac{\beta}{\Lambda}\, e^{-\lambda t}
  \int_0^t e^{\lambda \tau}\,n(\tau)\,d\tau,
\]

where \(n(t)\) is imported directly from the analytical solution of `Neutron_density_SciPyNumPy.py`.

The integral is evaluated numerically using `scipy.integrate.cumulative_trapezoid`.  
The script produces values of \(C(t)\) over a user-defined time mesh.

---

### **4. `RK4_reference_mpmath.py`**

This code implements a **fourth-order Runge–Kutta (RK4)** integrator with 32-digit precision to solve directly the NPKE system:

\[
\frac{dn}{dt}=\frac{\rho(t)-\beta}{\Lambda}n(t)+\lambda C(t)+q,
\quad
\frac{dC}{dt}=\frac{\beta}{\Lambda}n(t)-\lambda C(t),
\]

with \(\rho(t)=at+b\).

This solver provides high-accuracy reference values to validate the analytical solution reported in the manuscript.

---

## Software Requirements

- Python 3.10+  
- SciPy ≥ 1.10  
- NumPy ≥ 1.22  
- mpmath ≥ 1.3.0  

All codes were tested under a Windows 11 environment on a 3.8 GHz desktop with 32 GB of RAM.

---

## Financial Support

The authors appreciate the financial support received from the Consejo Nacional de Humanidades, Ciencia y Tecnología (CONAHCYT), under the program *Estancias Posdoctorales por México, 2022*, with the project entitled:  
*Desarrollo de modelos fenomenológicos energéticos de orden fraccional, para la optimización y simulación en reactores nucleares de potencia*,  
by which the present development was possible.

---

If you use or adapt this code, please cite the associated article once published.
