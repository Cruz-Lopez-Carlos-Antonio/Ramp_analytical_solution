---
layout: default
title: Main Equations
---

## Neutron Point Kinetics Equations

We consider the NPKE with a single group of delayed neutron precursors:

$$
\frac{dn}{dt} = \frac{\rho(t)-\beta}{\Lambda} n(t) + \lambda C(t) + q,
\qquad
\frac{dC}{dt} = \frac{\beta}{\Lambda} n(t) - \lambda C(t),
$$

with a linear ramp reactivity

$$
\rho(t) = a t + b.
$$

Here \(\Lambda\) is the prompt generation time, \(\beta\) is the total delayed neutron fraction, \(\lambda\) the precursor decay constant, and \(q\) an external source term.

---

## Integral Representation of \(n(t)\)

The analytical solution for the neutron density \(n(t)\) can be written as

$$
\begin{aligned}
n(t)=\;&
K_{1}\,e^{-\lambda t}
\int_0^\infty e^{-y^{2}/2+\mathcal{E}_{3}(t)\,y}\, y^{\lambda\beta/a}\,dy
\\[4pt]
&+K_{2}\,e^{-\lambda t}
\int_0^\infty e^{-y^{2}/2-\mathcal{E}_{3}(t)\,y}\,y^{\lambda\beta/a}\,dy
\\[4pt]
&+B_{p,1}
\int_0^\infty 
\exp\!\left[
-\frac{1}{a}\left(\frac{\Lambda}{2}s^{2}+(\beta-b-at)s\right)
\right]
\,(s+\lambda)^{\lambda\beta/a}\,ds.
\end{aligned}
$$

The functions \(\mathcal{E}_3(t)\) and \(B_{p,1}\), as well as the constants \(K_1\) and \(K_2\), are defined in the manuscript and are determined from the initial conditions \(n(0)\) and \(\dot n(0)\).

---

## Delayed Neutron Precursor Concentration

Once \(n(t)\) is known, the concentration of delayed neutron precursors is obtained from

$$
C(t)
= C(0)\,e^{-\lambda t}
+ \frac{\beta}{\Lambda}\,e^{-\lambda t}
\int_0^t e^{\lambda\tau}\,n(\tau)\,d\tau.
$$

The numerical implementation of this integral is described in the **Codes** section.

---

## RK4 Reference System

For validation purposes, a high-precision RK4 solver is applied directly to the NPKE system:

$$
\frac{dn}{dt}=\frac{\rho(t)-\beta}{\Lambda}\,n(t)+\lambda\,C(t)+q,
\qquad
\frac{dC}{dt}=\frac{\beta}{\Lambda}\,n(t)-\lambda\,C(t),
$$

with the same ramp reactivity \(\rho(t)=at+b\).  

Further details and numerical parameters are provided in the **Validation** page.
