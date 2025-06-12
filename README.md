Of course. Here is the provided report converted into a single, well-formatted Markdown file. The text, equations, code, tables, and figures have been preserved and structured appropriately.

***

# THEORETICAL ASTROPHYSICS LAB - EXERCISE 1
## Planetary Orbits Around A Star

**Sohaib Ali**¹  
*MSc Physics and Astronomy (PandA)*

**Prof. Krzysztof Gęsicki**

---

**Abstract**—This report employs Newton's laws of motion and law of gravitation to calculate the elliptical motion of a planet around the Sun. These laws were used to develop a numerical integrator of leap-frog type [2]. Various orbits have been predicted by changing initial conditions for the integrator. The integrator was developed using Python.

## 1. Introduction

For a body orbiting a central mass, Newton's law of universal gravitation states:

$$
F = -G \frac{Mm}{r^2} \hat{r} \quad (1)
$$

Where F is the force vector, G is the gravitational constant, M and m are the masses of central object and orbiting body respectively, r is the distance between them, **r̂** is the unit vector pointing along the line connecting the two masses, and the negative sign indicates that the force is attractive (directed toward the center). In this exercise, to simplify things and minimize numerical errors, we use normalized units G=M=1. This reduces (1) to:

$$
a = -\frac{\hat{r}}{r^2} \quad (2)
$$

Moreover, if we know the position and velocity of a body at a specific time t, we can use Newton's laws of motion to predict the trajectory of a particle in an N-dimensional space.

## 2. Methodology

We start by defining boundary conditions for the planet. We assume a two-dimensional space with x₀ and y₀ as initial position coordinates of the planet and the corresponding velocity components vₓ₀ and vᵧ₀. We also assume that G and M in (1) are dimensionless and are equal to 1. Based on these initial positions and assumptions, I calculated r₀, vₓ₀,ₕₐₗf, vᵧ₀,ₕₐₗf, aₓ₀, aᵧ₀ using the following equations:

$$
r_0 = \sqrt{(x_0^2 + y_0^2)} \quad (3)
$$

$$
a_{x0} = -x_0/r_0^3 \quad (4)
$$

$$
a_{y0} = -y_0/r_0^3 \quad (5)
$$

$$
v_{x0,half} = v_{x0} + a_{x0}(dt/2) \quad (6)
$$

$$
v_{y0,half} = v_{y0} + a_{y0}(dt/2) \quad (7)
$$

Where r₀ is radial position at x₀ and y₀, aₓ, and aᵧ, are accelerations due gravity at x₀ and y₀. Eqs. 4 and 5 are derived from 2. This was done to simplify our calculations for the orbit. As a consequence, no gravitational interaction was considered between the Sun and the planet. The half step velocities defined by 6 and 7 are a signature element of Leapfrog integration scheme in orbital dynamics. As per leapfrog method, half-step velocities allow the positions to be updated using a velocity value that represents an average between the current and future states. Feynmann [1] discusses the need to have mid-point velocities between two time-steps in order to not "miss-out" on the information about the motion of particles while computing its motion. By computing velocities at half-time steps, the algorithm effectively centers the update of positions and velocities in time, which leads to a local truncation error of order 3 and an overall second-order accurate integration.

## 3. Code Implementation
The above equations were translated into Python and a leapfrog integrator in the form of a python function was developed [3].

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 14:29:09 2025
@author: sohaib
"""
import numpy as np
from tabulate import tabulate
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def leapfrog(dt, steps, x0, y0, vx0, vy0):
    r0 = np.sqrt(x0**2 + y0**2)  # radial position of planet at x0, y0
    print("r0 :", r0)

    ax0 = -x0 / r0**3  # acceleration component at x0
    ay0 = -y0 / r0**3  # acceleration component at y0
    print("ax(0) :", ax0)
    print("ay(0) :", ay0)

    vxhalf = vx0 + ax0 * (dt / 2)  # half integer time step velocity in x direction
    vyhalf = vy0 + ay0 * (dt / 2)  # half integer time step velocity in y direction
    print("vx(e/2) : ", vxhalf)
    print("vy(e/2) :", vyhalf)

    time = 0.0  # started at time 0
    
    # bunch of initializations for arrays to store data later
    results = []
    results.append([time, x0, y0, ax0, ay0, vx0, vy0, r0])

    # just to match the table in Feynman's book, I 
    # included an additional row in the data table to 
    # show half step velocities
    results.append([None, None, None, None, None, vxhalf, vyhalf, None])

    xcurr = x0
    ycurr = y0
    vxhalf_curr = vxhalf
    vyhalf_curr = vyhalf

    # looping and calculating parameters until time 
    # interval reaches the limit
    for n in range(steps):
        xnew = xcurr + vxhalf_curr * dt
        ynew = ycurr + vyhalf_curr * dt
        time += dt
        
        rnew = np.sqrt(xnew**2 + ynew**2)
        axnew = -xnew / rnew**3
        aynew = -ynew / rnew**3

        vxnew = vxhalf_curr + axnew * (dt / 2)
        vynew = vyhalf_curr + aynew * (dt / 2)

        results.append([time, xnew, ynew, axnew, aynew, vxnew, vynew, rnew])

        vxhalf_curr = vxnew + (dt / 2) * axnew
        vyhalf_curr = vynew + (dt / 2) * aynew
        
        results.append([None, None, None, None, None, vxhalf_curr, vyhalf_curr, None])
        
        xcurr = xnew
        ycurr = ynew

    # to print the table
    headers = ["time", "x", "y", "ax", "ay", "vx", "vy", "r"]
    print(tabulate(
        results,
        headers=headers,
        tablefmt="plain",
        floatfmt=".3f",
        missingval=""
    ))

    # since i used seaborn to plot, converting the stored
    # values in arrays to a pandas df was useful for plotting
    df = pd.DataFrame(results, columns=["time", "x", "y", "ax", "ay", "vx", "vy", "r"])
    
    ax = sns.lineplot(data=df, x="x", y="y", marker=False, legend=0, sort=False)
    sns.scatterplot(x=[0], y=[0], color="yellow", s=150, edgecolor="black", ax=ax)
    
    ax.set(xlabel='x position', ylabel='y position')
    ax.set(title=f'(dt: {dt}, steps: {steps}, x0: {x0}, y0: {y0}, vx0: {vx0}, vy0: {vy0})')
    plt.suptitle('Motion of a planet around Sun')
    plt.figure(dpi=300)
    plt.show()

    return

# 1st iteration - a coarse timestep
# dt=0.1
# steps=2000
# x0=0.500
# y0=0.000
# vx0=0
# vy0=1.630
# leapfrog(dt, steps, x0, y0, vx0, vy0)

# 2nd iteration
dt = 0.01
steps = 2000
x0 = 0.500
y0 = 0.000
vx0 = 0
vy0 = 1.630
leapfrog(dt, steps, x0, y0, vx0, vy0)
```
Using a coarse time step dt=0.100 in the first iteration leads to a Rosette-shaped orbit (Figure 1). A pseudo precession of the orbit often called numerical precession. The leapfrog method assumes that velocities at half-steps are good approximations for updating positions. However, if dt is too large, the velocity updates become less accurate. Due to a coarse time-step over a large number of steps i.e. 2000, the planet arrives at a slightly different position at the end of each iteration of loop. This leads to an open ellipse, consequently creating a precession-like effect that isn't physical. To make the large dt conform to the physical expectation, decreasing step-size dt to finer values e.g. dt=0.01 generates a perfectly elliptical orbit as shown in 2.

The values of each parameters calculated using eqs. 2, 3, 4, 5, 6 and 7 during each iteration are shown in truncated tables 1 and 2. As expected, we see more precise half-step velocities in 2, meaning that we can now model the position of the planet much more accurately within each time-step. This consequently removes the observed numerical precession and generates a perfectly elliptical orbit.

Other combination of initial conditions also yield some interesting orbits shown in figure 3. While figures 3a and 3b are obtained just by exchanging the initial velocities and positions from x to y, figures 3c represents a special orbit where increasing the velocity beyond a certain limit makes the Hamiltonian of the system positive, thereby, causing the planet to escape the system. 3d shows the limitations of our integrator that even though finer time-steps lead to more realistic solutions, in some combination of initial conditions, errors continue to add up which consequently lead to "drifts" nonetheless.

## 4. Figures






*(a) Another example of a closed orbit at y₀ = 0.5*  
*(b) Changing from a non-zero vᵧ to a non-zero vₓ*  
*(c) An example of an unbound/hyperbolic orbit by increasing velocity (positive total energy)*  
*(d) Limitations of the 2nd-order leapfrog method leading to drifts*

## 5. Simulation Data

### 5.1. Iteration 1 - Coarse time-step
The simulation was initialized with r₀ = 0.5, aₓ(0) = -4.0, aᵧ(0) = -0.0, vₓ(ε/2) = -0.02, and vᵧ(ε/2) = 1.63. Table 1 presents the time evolution data. The second row after each full step contains the half-step velocity values.

**Table 1. Iteration 1 with dt=0.1**
| Time  | x       | y       | aₓ      | aᵧ      | vₓ      | vᵧ      | r       |
| :---- | :------ | :------ | :------ | :------ | :------ | :------ | :------ |
| 0.000 | 0.500   | 0.000   | -4.000  | -0.000  | 0.000   | 1.630   | 0.500   |
|       |         |         |         |         | -0.200  | 1.630   |         |
| 0.100 | 0.480   | 0.163   | -3.685  | -1.251  | -0.384  | 1.567   | 0.507   |
|       |         |         |         |         | -0.568  | 1.505   |         |
| 0.200 | 0.423   | 0.313   | -2.897  | -2.146  | -0.713  | 1.398   | 0.527   |
|       |         |         |         |         | -0.858  | 1.290   |         |
| 0.300 | 0.337   | 0.443   | -1.958  | -2.569  | -0.956  | 1.162   | 0.556   |
|       |         |         |         |         | -1.054  | 1.033   |         |
| 0.400 | 0.232   | 0.546   | -1.112  | -2.617  | -1.110  | 0.903   | 0.593   |
|       |         |         |         |         | -1.165  | 0.772   |         |
| ...   | ...     | ...     | ...     | ...     | ...     | ...     | ...     |
| 2.100 | -1.022  | 0.002   | 0.957   | -0.002  | 0.009   | -0.797  | 1.022   |
|       |         |         |         |         | 0.057   | -0.797  |         |
*Note: The second row after each full time step contains half-step velocity values.*

### 5.2. Iteration 2 - Fine time-step
The simulation was initialized with r₀ = 0.5, aₓ(0) = -4.0, aᵧ(0) = -0.0, vₓ(ε/2) = -0.2, and vᵧ(ε/2) = 1.63. Table 2 presents the time evolution data. Note that the second row after each full time step contains the half-step velocity values.

**Table 2. Iteration 2 with dt=0.01**
| Time  | x       | y       | aₓ      | aᵧ      | vₓ      | vᵧ      | r       |
| :---- | :------ | :------ | :------ | :------ | :------ | :------ | :------ |
| 0.000 | 0.500   | 0.000   | -4.000  | -0.000  | 0.000   | 1.630   | 0.500   |
|       |         |         |         |         | -0.020  | 1.630   |         |
| 0.010 | 0.500   | 0.016   | -3.997  | -0.130  | -0.040  | 1.629   | 0.500   |
|       |         |         |         |         | -0.060  | 1.629   |         |
| 0.020 | 0.499   | 0.033   | -3.987  | -0.260  | -0.080  | 1.627   | 0.500   |
|       |         |         |         |         | -0.100  | 1.626   |         |
| 0.030 | 0.498   | 0.049   | -3.972  | -0.389  | -0.120  | 1.624   | 0.501   |
|       |         |         |         |         | -0.140  | 1.622   |         |
| 0.040 | 0.497   | 0.065   | -3.950  | -0.517  | -0.159  | 1.620   | 0.501   |
|       |         |         |         |         | -0.179  | 1.617   |         |
| ...   | ...     | ...     | ...     | ...     | ...     | ...     | ...     |
| 0.210 | 0.417   | 0.324   | -2.839  | -2.205  | -0.752  | 1.372   | 0.527   |
|       |         |         |         |         | -0.767  | 1.361   |         |
*Note: The second row after each full step contains half-step velocity values.*

---
## References
[1] Feynman, R. P., Leighton, R. B., & Sands, M. (1964). *The Feynman Lectures on Physics: Volume I, Chapter 9*. Available at: https://www.feynmanlectures.caltech.edu/I_09.html (Accessed: March 10, 2025).

[2] Hairer, E., Lubich, C., & Wanner, G. (2006). *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations*. Springer, Berlin.

[3] Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press, Cambridge.

---
¹MSc Physics and Astronomy (PandA)
