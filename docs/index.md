---
layout: default
title: Overview
---

## Overview

This repository documents a new analytical solution for the Neutron Point Kinetics Equations (NPKE) 

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\frac{dn\left(t\right)}{dt}=\frac{\rho\left(t\right)-\beta}{\Lambda}n\left(t\right)+\lambda C\left(t\right)+q,\\
\frac{dC\left(t\right)}{dt}=\frac{\beta}{\Lambda}n\left(t\right)-\lambda C\left(t\right).\   
$$
</div>

with a linear ramp reactivity of the form

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
\rho(t) = a t + b, \quad a>0,
$$
</div>

Using the Modified Integral Method (MIM). The corresponding Python 3 codes implement:

- An analytical representation of the neutron density $$n(t)$$ using SciPy/NumPy libraries,
- An analytical representation of the neutron density $$n(t)$$ using the mpmath library with 32-digit precision by default.
- The delayed neutron precursor concentration $$C(t)$$ 
- And a high-precision Runge–Kutta 4 (RK4) reference solver.
- The analytical solution developed by Zhang et al. (2008)
- The analytical solution developed by Palma et al. (2010)

These results accompany the manuscript:  
***Analytical solutions of the Neutron Point Kinetics Equations under a linear reactivity ramp***,  
that was recently submitted to the journal *Computer Physics Communications*.

You can use the navigation bar above to explore:

<div style="background:#f8f8f8; padding:15px; border-radius:10px; border:1px solid #d0d0d0;">
  <h4 style="margin-top:0;">Analytical expression for \(n(t)\)</h4>
  $$
  n(t)=K_1 e^{-\lambda t}\int_0^\infty e^{-y^2/2 + \mathcal{E}_3(t)\,y}y^{\lambda\beta/a}\,dy
  $$
  <ul style="margin-top:10px;">
    <li>Derived from the reduced form of the NPKE</li>
    <li>Uses a change of variables and integral transforms</li>
    <li>Depends on constants \(K_1, K_2\) fixed by the initial conditions</li>
  </ul>
</div>

- **Equations** – main analytical expressions used in the article,  
- **Codes** – a summary of the Python scripts in the repository,  
- **Validation** – comparison with a high-precision RK4 benchmark,  
- **About** – authorship and financial support.

For more details of the derivation and the computational implementation, please see the manuscript. 
