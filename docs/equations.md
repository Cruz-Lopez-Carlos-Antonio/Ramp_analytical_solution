---
layout: default
title: Main Equations
---
# Overview of the methodology


## Neutron Point Kinetics Equations

We consider the NPKE with a single group of delayed neutron precursors and a constant source $$q$$:

<div style="background:#f7f7f7; padding:15px; border-left:4px solid #4a90e2; border-radius:6px; margin:20px 0;">
$$
n(t)=K_1 e^{-\lambda t}\int_0^\infty e^{-y^2/2 + \mathcal{E}_3(t)\,y}y^{\lambda\beta/a}\,dy
$$
</div>

$$
\frac{dn}{dt} = \frac{\rho(t)-\beta}{\Lambda} n(t) + \lambda C(t) + q,
\qquad
\frac{dC}{dt} = \frac{\beta}{\Lambda} n(t) - \lambda C(t),
$$

with a positive linear ramp reactivity

$$
\rho(t) = a t + b,  \quad a>0.
$$

Here $$\Lambda$$ is the prompt generation time, $$\beta$$ is the total delayed neutron fraction and $$\lambda$$.

---

## Modified Integration Method (MIM) and Integral Representation of $$n(t)$$

The analytical solution for the neutron density $$n(t)$$ can be written as

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

For completeness, we recall that the auxiliary function $$\mathcal{E}_3(t)$$ arises directly from the analytical reduction of the point kinetics equation with a linear reactivity. It is defined as

$$
\mathcal{E}_3(t)
=
\frac{1}{\sqrt{2a\Lambda}}
\left[
at + b - \beta + \frac{\Lambda\lambda}{2}
\right],
$$

which is obtained after completing the square in the exponent of the transformed equation for $n(t)$. The function $B_{p,1}$ appearing in the third integral term of $$n(t)$$ follows analogously from the same reduction procedure and depends linearly on the parameters $$a$$, $$b$$, $$\beta$$, $$\Lambda$$, $$\lambda$$ and $$q$$, as detailed in the associated manuscript.

---

## Delayed Neutron Precursor Concentration

Once $$n(t)$$ is known, the concentration of delayed neutron precursors is obtained from

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

with the same ramp reactivity $$\rho(t)=at+b$$.  

Further details and numerical parameters are provided in the **Validation** page.
