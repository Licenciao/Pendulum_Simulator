# Pendulum Simulator

A didactic collection of numerical simulations in **C (SDL2 / Raylib)** and **Python (Matplotlib)** designed to study and visualize the nonlinear and chaotic dynamics of pendulum systems in 2D and 3D.

---

## Project Proposal

Even simple mechanical rules can produce rich, unpredictable motion. This repository provides an accessible, hands-on way to explore how pendulum systems evolve over time—from a classic double pendulum to an arbitrary chain of $N$ coupled pendulums and an elastic double pendulum suspended by spider-silk threads that can stretch, go slack, or snap under tension.

---

## Mathematical & Numerical Foundations

Across the three main 2D models, the simulator combines analytical mechanics with numerical methods:

1. **Lagrangian Mechanics (Rigid Models):** Instead of solving constraint forces directly, the rigid-rod models use generalized angular coordinates $\theta_i$ (measured from the downward vertical). Applying the Euler–Lagrange equations yields a coupled linear system at every instant:
   $$M(\theta)\,\alpha = \text{RHS}(\theta, \omega)$$
   where $M(\theta)$ is the configuration-dependent **mass matrix**, $\alpha = \ddot{\theta}$ is the vector of angular accelerations, and $\text{RHS}$ groups centrifugal, gravitational, and viscous joint-damping torques.
2. **Newtonian Mechanics in Cartesian Space (Elastic Model):** Because the arm lengths are no longer fixed, the elastic pendulum tracks positions and velocities $(x_i, y_i, v_{xi}, v_{yi})$ directly under gravity, viscoelastic thread tension, and joint friction.
3. **4th-Order Runge–Kutta (RK4) Integration:** All three models advance the state vector $\mathbf{y}$ in time using the classic RK4 scheme:
   $$\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{\Delta t}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)$$
   For stiff configurations (large $N$ or high elastic stiffness), the time step $\Delta t$ is automatically scaled down to preserve numerical stability and energy accuracy.

---

## Physical Models & Governing Equations

### 1. Double Pendulum (`pendulo_doble_c.c`, `pendulo_doble.py`)

Two point masses $m_1, m_2$ connected by rigid massless rods of lengths $l_1, l_2$, with viscous rotational damping coefficients $b_1, b_2$ at the joints. Using $\Delta = \theta_1 - \theta_2$, the $2 \times 2$ equations of motion are:

$$
\begin{bmatrix}
(m_1 + m_2)l_1^2 & m_2 l_1 l_2 \cos\Delta \\
m_2 l_1 l_2 \cos\Delta & m_2 l_2^2
\end{bmatrix}
\begin{bmatrix}
\ddot{\theta}_1 \\
\ddot{\theta}_2
\end{bmatrix}
=
\begin{bmatrix}
-m_2 l_1 l_2 \sin\Delta\,\dot{\theta}_2^2 - (m_1 + m_2)g l_1 \sin\theta_1 - b_1\dot{\theta}_1 \\
+m_2 l_1 l_2 \sin\Delta\,\dot{\theta}_1^2 - m_2 g l_2 \sin\theta_2 - b_2\dot{\theta}_2
\end{bmatrix}
$$

- **Mechanical Energy:**
  $$E = T + V = \frac{1}{2}(m_1 + m_2)l_1^2\dot{\theta}_1^2 + \frac{1}{2}m_2 l_2^2\dot{\theta}_2^2 + m_2 l_1 l_2 \dot{\theta}_1\dot{\theta}_2\cos\Delta - (m_1 + m_2)g l_1 \cos\theta_1 - m_2 g l_2 \cos\theta_2$$

---

### 2. $N$-Link Pendulum (`pendulo_n_c.c`)

Generalizes the double pendulum to a chain of $N$ point masses (setting $N = 2$ reproduces the double pendulum exactly). Defining the cumulative hanging mass from link $i$ downward as:

$$S_i = \sum_{k=i}^{N} m_k$$

the coupled $N \times N$ system $M(\theta)\ddot{\theta} = \text{RHS}$ has entries:

$$M_{ij} = S_{\max(i,j)}\, l_i l_j \cos(\theta_i - \theta_j)$$

$$\text{RHS}_i = -\sum_{j=1}^{N} S_{\max(i,j)}\, l_i l_j \sin(\theta_i - \theta_j)\,\dot{\theta}_j^2 - S_i g l_i \sin\theta_i - b_i\dot{\theta}_i$$

- **Linear Solver:** Because $M(\theta)$ is symmetric and positive-definite, the angular accelerations $\ddot{\theta}$ are solved at each RK4 evaluation via **in-place Cholesky factorization** ($M = L L^T$), followed by forward and backward substitution.
- **Mechanical Energy:**
  $$E = \frac{1}{2}\dot{\theta}^T M(\theta)\dot{\theta} - g\sum_{i=1}^{N} S_i l_i \cos\theta_i$$

---

### 3. Elastic Double Pendulum (`pendulo_doble_elastico_c.c`)

Replaces the rigid rods with extensible *Nephila* dragline spider-silk threads of cross-sectional area $A$, Young's modulus $E \approx 10\text{ GPa}$, and tensile strength $\sigma_{\max} \approx 1.1\text{ GPa}$.

- **Stiffness, Axial Damping & Rupture Limit:**
  $$k_i = \frac{E\,A}{l_i}, \qquad c_1 = 2\zeta\sqrt{k_1(m_1 + m_2)}, \qquad c_2 = 2\zeta\sqrt{k_2 m_2}, \qquad T_{\max} = \sigma_{\max} A$$
- **Thread Tension (Tension-Only / Slack Condition):**
  For instantaneous arm vectors $\mathbf{r}_1 = (x_1, y_1)$ and $\mathbf{r}_2 = (x_2 - x_1, y_2 - y_1)$ with lengths $L_i = \|\mathbf{r}_i\|$, unit vectors $\hat{\mathbf{u}}_i = \mathbf{r}_i / L_i$, elongation $\Delta L_i = L_i - l_i$, and radial stretching rate $v_{r,i} = \dot{\mathbf{r}}_i \cdot \hat{\mathbf{u}}_i$:
  $$T_i = \begin{cases} \max\left(0,\; k_i \Delta L_i + c_i v_{r,i}\right) & \text{if } \Delta L_i > 0 \\ 0 & \text{if } \Delta L_i \le 0 \text{ (slack thread)} \end{cases}$$
  If $T_i > T_{\max}$ at any integration step, the thread snaps and the simulation reports which pendulum broke.
- **Newton's Equations of Motion:**
  Including joint viscous damping $\mathbf{F}_{t,i} = -\frac{b_i \omega_i}{L_i}\hat{\mathbf{t}}_i$ (where $\omega_i = (\mathbf{r}_i \times \dot{\mathbf{r}}_i)_z / L_i^2$ and $\hat{\mathbf{t}}_i = (-u_{iy}, u_{ix})$):
  $$m_1 \ddot{\mathbf{r}}_1 = -T_1\hat{\mathbf{u}}_1 + T_2\hat{\mathbf{u}}_2 - m_1 g\,\hat{\mathbf{j}} + \mathbf{F}_{t,1}$$
  $$m_2 (\ddot{\mathbf{r}}_1 + \ddot{\mathbf{r}}_2) = -T_2\hat{\mathbf{u}}_2 - m_2 g\,\hat{\mathbf{j}} + \mathbf{F}_{t,2}$$
- **Mechanical Energy:**
  $$E = \frac{1}{2}m_1 v_1^2 + \frac{1}{2}m_2 v_2^2 + m_1 g y_1 + m_2 g y_2 + \frac{1}{2}k_1 \max(0, \Delta L_1)^2 + \frac{1}{2}k_2 \max(0, \Delta L_2)^2$$
