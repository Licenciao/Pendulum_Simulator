# Pendulum Simulator

An interactive collection of **C** and **Python** simulations built to visualize how pendulums move, interact, and become chaotic in 2D and 3D.

## About the Project

Even simple pendulums become unpredictable when linked together. This project lets you experiment with masses, lengths, gravity, and friction to see how small changes in starting conditions completely transform the motion.

## How the Math Works

At each tiny time step $\Delta t$, the simulator calculates the forces acting on each mass to find its acceleration, and then updates velocities and positions using the **Runge–Kutta (RK4)** method:

$$y_{n+1} = y_n + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

---

## Models & Core Equations

### 1. Double Pendulum (`pendulo_doble_c.c`, `pendulo_doble.py`)
Two rigid arms of lengths $l_1, l_2$ and masses $m_1, m_2$ swinging at angles $\theta_1, \theta_2$ (with $\Delta = \theta_1 - \theta_2$ and joint friction $b_1, b_2$). Both arms pull on each other through this coupled system:

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

### 2. $N$-Pendulum (`pendulo_n_c.c`)
Extends the double pendulum to a chain of $N$ connected arms. Using $S_i = \sum_{k=i}^{N} m_k$ as the total mass hanging from link $i$ downward, the system solves $M\,\ddot{\theta} = \text{RHS}$ (via Cholesky decomposition):

$$M_{ij} = S_{\max(i,j)}\, l_i l_j \cos(\theta_i - \theta_j)$$

$$\text{RHS}_i = -\sum_{j=1}^{N} S_{\max(i,j)}\, l_i l_j \sin(\theta_i - \theta_j)\,\dot{\theta}_j^2 - S_i g l_i \sin\theta_i - b_i\dot{\theta}_i$$

### 3. Elastic Double Pendulum (`pendulo_doble_elastico_c.c`)
Replaces rigid rods with stretchy spider-silk threads (Young's modulus $E$, cross-section $A$, and breaking strength $\sigma_{\max}$). Because the arms stretch, motion is tracked in $(x, y)$ coordinates using Newton's second law ($\mathbf{F} = m\mathbf{a}$):

- **Spring stiffness & breaking limit:**
  $$k_i = \frac{E\,A}{l_i}, \qquad T_{\max} = \sigma_{\max} A$$
- **Thread tension** (threads only pull when stretched by $\Delta L_i = L_i - l_i > 0$ and go slack otherwise):
  $$T_i = \max\left(0,\; k_i \Delta L_i + c_i \dot{L}_i\right)$$
  If $T_i > T_{\max}$, the silk snaps and the simulation stops.
