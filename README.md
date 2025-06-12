# Motion of Planets Around a Star

This repository contains:

- **Code:** `leapfrog.py` — a 2nd-order leapfrog integrator for planetary orbits.  
- **Figures:** located in `images/`  

---

## Abstract  
This study employs normalized units (G = M = 1) and Newton’s laws to simulate a planet’s motion around the Sun. A leapfrog integrator was implemented in Python to explore how varying initial conditions and time-step sizes affect orbital trajectories :contentReference[oaicite:0]{index=0}.

## 1. Introduction  
For a body orbiting a central mass, Newton’s law of gravitation gives the acceleration  
\[
\mathbf{a} = -\frac{\hat{\mathbf{r}}}{r^2}
\]  
where \(\hat{\mathbf{r}}\) is the unit vector from the planet to the star. By choosing \(G = M = 1\), we reduce the governing equations to this form. Given initial position \((x_0,y_0)\) and velocity \((v_{x0},v_{y0})\), one can predict trajectories using numerical integration :contentReference[oaicite:1]{index=1}.

## 2. Methodology  
1. **Initial conditions:**  
   - Position: \((x_0,\,y_0)\)  
   - Velocity: \((v_{x0},\,v_{y0})\)  
2. **Leapfrog scheme:**  
   - Compute radial distance  
     \[
       r_0 = \sqrt{x_0^2 + y_0^2}
     \]  
   - Compute accelerations  
     \[
       a_{x0} = -\frac{x_0}{r_0^3}
       \quad,\quad
       a_{y0} = -\frac{y_0}{r_0^3}
     \]  
   - Compute half-step velocities  
     \[
       v_{x,\tfrac12} = v_{x0} + a_{x0}\,\frac{\Delta t}{2},
       \quad
       v_{y,\tfrac12} = v_{y0} + a_{y0}\,\frac{\Delta t}{2}
     \]  
   - Update positions and velocities over each full step, appending both full-step and half-step values to the results table :contentReference[oaicite:2]{index=2}.

## 3. Code Implementation  
```python
def leapfrog(dt, steps, x0, y0, vx0, vy0):
    r0 = np.sqrt(x0**2 + y0**2)
    ax0 = -x0/r0**3;  ay0 = -y0/r0**3
    vx_half = vx0 + ax0*(dt/2);  vy_half = vy0 + ay0*(dt/2)
    time = 0.0
    results = [[time, x0, y0, ax0, ay0, vx0, vy0, r0],
               [None, None, None, None, None, vx_half, vy_half, None]]
    x, y = x0, y0
    vxh, vyh = vx_half, vy_half

    for _ in range(steps):
        x_new = x + vxh*dt
        y_new = y + vyh*dt
        time += dt
        r_new = np.sqrt(x_new**2 + y_new**2)
        ax_new = -x_new/r_new**3;  ay_new = -y_new/r_new**3
        vx_new = vxh + ax_new*(dt/2);  vy_new = vyh + ay_new*(dt/2)
        results.append([time, x_new, y_new, ax_new, ay_new, vx_new, vy_new, r_new])
        vxh = vx_new + ax_new*(dt/2);  vyh = vy_new + ay_new*(dt/2)
        results.append([None, None, None, None, None, vxh, vyh, None])
        x, y = x_new, y_new
    return results
