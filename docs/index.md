---
layout: default
title: Overview
---

## Overview

This site documents a new analytical solution for the Neutron Point Kinetics Equations (NPKE) with a linear ramp reactivity of the form

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

- **Equations** – main analytical expressions used in the article,  
- **Codes** – a summary of the Python scripts in the repository,  
- **Validation** – comparison with a high-precision RK4 benchmark,  
- **About** – authorship and financial support.

For more details of the derivation and the computational implementation, please see the manuscript. 
