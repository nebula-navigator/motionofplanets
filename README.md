# Planetary Orbits Around a Star

**Sohaib Ali**  
MSc Physics and Astronomy (PandA)  
Theoretical Astrophysics Lab — Exercise 1  
Prof. Krzysztof Gęsicki  
Nicolaus Copernicus University, March 11, 2022  

---

## 📄 Abstract

This report employs Newton's laws of motion and gravitation to calculate the elliptical motion of a planet around the Sun using a numerical integrator of the Leapfrog type. Various orbits are simulated by modifying the initial conditions. All simulations were performed in Python.

---

## 🔬 Introduction

Newton's law of universal gravitation for a planet orbiting a star is:

$$
\vec{F} = -\frac{GMm}{r^2} \hat{r}
$$

To simplify calculations, we use normalized units: \( G = M = 1 \), reducing the equation of motion to:

$$
\vec{a} = -\frac{\hat{r}}{r^2}
$$

This allows us to model the motion in a dimensionless 2D space.

---

## ⚙️ Methodology

We define initial positions \((x_0, y_0)\), velocities \((v_{x0}, v_{y0})\), and compute:

- Radial distance:
  $$
  r_0 = \sqrt{x_0^2 + y_0^2}
  $$

- Acceleration components:
  $$
  a_{x0} = -\frac{x_0}{r_0^3}, \quad a_{y0} = -\frac{y_0}{r_0^3}
  $$

- Half-step velocities (Leapfrog scheme):
  $$
  v_{x0}^{1/2} = v_{x0} + a_{x0} \cdot \frac{dt}{2}, \quad v_{y0}^{1/2} = v_{y0} + a_{y0} \cdot \frac{dt}{2}
  $$

The Leapfrog method advances positions using these half-step velocities, offering second-order accuracy and time-reversibility.

---

## 🧪 Code Overview

The simulation was implemented in Python using:

- `numpy` for vector math
- `matplotlib` and `seaborn` for plotting
- `pandas` and `tabulate` for tabular outputs

Leapfrog integration was implemented in a custom function, updating position and velocity iteratively for a given number of time steps.

---

## 📊 Results

### Elliptical Orbit (dt = 0.01)

![Elliptical Orbit](./images/fig-000.png)

### Rosette Orbit (dt = 0.1)

![Rosette Orbit](./images/fig-001.png)

The orbit morphology depends sensitively on the chosen time step and initial velocity vector. At larger `dt`, the orbit deviates from a perfect ellipse due to numerical error, forming rosette-like paths.

---

## 📎 Full Lab Report

You can read the detailed lab report [here](./SohaibAli_Report1_TAL.pdf).

---

## 📚 References

1. Feynman, R. P., Leighton, R. B., & Sands, M. (2010). *The Feynman Lectures on Physics*.
2. Hut, P., Makino, J., & McMillan, S. (1995). *Building a Better Leapfrog*. ApJL, 443, L93.

