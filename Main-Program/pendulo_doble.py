import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # Backend necesario para la animación en vivo
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# ============================================================
# SIMULADOR DE PÉNDULO DOBLE
# Ángulos medidos desde la vertical hacia abajo.
# Incluye fricción viscosa independiente en cada articulación.
# No requiere SciPy: usa integración RK4 implementada aquí.
# ============================================================

print("--- Simulador de Péndulo Doble ---")

# --- Condiciones iniciales ---
A1_user = float(input("Ángulo inicial del brazo 1 (grados, ej: 120): "))
A2_user = float(input("Ángulo inicial del brazo 2 (grados, ej: -10): "))

w1_user = float(input("Velocidad angular inicial 1 (grados/s, ej: 0): "))
w2_user = float(input("Velocidad angular inicial 2 (grados/s, ej: 0): "))

# --- Geometría ---
l1 = float(input("Longitud del brazo 1 (m, ej: 1.0): "))
l2 = float(input("Longitud del brazo 2 (m, ej: 1.0): "))

# --- Parámetros físicos ---
m1 = float(input("Masa 1 (kg, ej: 1.0): "))
m2 = float(input("Masa 2 (kg, ej: 1.0): "))
g = float(input("Gravedad (m/s**2, ej: 9.81): "))

# b1 y b2 son coeficientes de amortiguamiento rotacional:
# torque de fricción = -b * velocidad_angular
b1 = float(input("Fricción articulación 1 (N*m*s/rad, ej: 0.02): "))
b2 = float(input("Fricción articulación 2 (N*m*s/rad, ej: 0.02): "))

t_max = float(input("Tiempo de simulación (s, ej: 20): "))


# --- Conversión a radianes ---
theta1_0 = np.radians(A1_user)
theta2_0 = np.radians(A2_user)
omega1_0 = np.radians(w1_user)
omega2_0 = np.radians(w2_user)


# ============================================================
# ECUACIONES DEL PÉNDULO DOBLE
# Estado:
# y = [theta1, omega1, theta2, omega2]
#
# Las aceleraciones se obtienen resolviendo, en cada instante,
# el sistema acoplado:
#
#     M(q) * q_ddot = fuerzas_generalizadas
#
# Esto permite incluir masas, longitudes y fricción sin usar
# la aproximación de ángulo pequeño del péndulo simple.
# ============================================================

def derivadas(estado):
    theta1, omega1, theta2, omega2 = estado

    delta = theta1 - theta2
    c = np.cos(delta)
    s = np.sin(delta)

    # Matriz de masa del sistema
    M11 = (m1 + m2) * l1**2
    M12 = m2 * l1 * l2 * c
    M21 = M12
    M22 = m2 * l2**2

    # Términos gravitacionales, centrífugos y de fricción
    rhs1 = (
        -m2 * l1 * l2 * s * omega2**2
        -(m1 + m2) * g * l1 * np.sin(theta1)
        -b1 * omega1
    )

    rhs2 = (
        +m2 * l1 * l2 * s * omega1**2
        -m2 * g * l2 * np.sin(theta2)
        -b2 * omega2
    )

    aceleraciones = np.linalg.solve(
        np.array([[M11, M12],
                  [M21, M22]], dtype=float),
        np.array([rhs1, rhs2], dtype=float)
    )

    alpha1, alpha2 = aceleraciones

    return np.array([
        omega1,
        alpha1,
        omega2,
        alpha2
    ], dtype=float)


# --- Integrador Runge-Kutta de cuarto orden ---
def paso_rk4(estado, dt):
    k1 = derivadas(estado)
    k2 = derivadas(estado + 0.5 * dt * k1)
    k3 = derivadas(estado + 0.5 * dt * k2)
    k4 = derivadas(estado + dt * k3)

    return estado + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)


# ============================================================
# SIMULACIÓN NUMÉRICA
# ============================================================

dt = 0.005
tiempos = np.arange(0.0, t_max + dt, dt)

estados = np.zeros((len(tiempos), 4), dtype=float)
estados[0] = [theta1_0, omega1_0, theta2_0, omega2_0]

for i in range(1, len(tiempos)):
    estados[i] = paso_rk4(estados[i - 1], dt)

theta1 = estados[:, 0]
omega1 = estados[:, 1]
theta2 = estados[:, 2]
omega2 = estados[:, 3]


# ============================================================
# CONVERSIÓN DE ÁNGULOS A POSICIONES X-Y
# ============================================================

x1 = l1 * np.sin(theta1)
y1 = -l1 * np.cos(theta1)

x2 = x1 + l2 * np.sin(theta2)
y2 = y1 - l2 * np.cos(theta2)


# ============================================================
# ENERGÍA MECÁNICA
# Útil para comprobar el comportamiento:
# - sin fricción debe mantenerse aproximadamente constante;
# - con fricción debe disminuir.
# ============================================================

def energia_mecanica(th1, om1, th2, om2):
    T = (
        0.5 * (m1 + m2) * l1**2 * om1**2
        + 0.5 * m2 * l2**2 * om2**2
        + m2 * l1 * l2 * om1 * om2 * np.cos(th1 - th2)
    )

    V = (
        -(m1 + m2) * g * l1 * np.cos(th1)
        -m2 * g * l2 * np.cos(th2)
    )

    return T + V


energia = energia_mecanica(theta1, omega1, theta2, omega2)


# ============================================================
# CONFIGURACIÓN VISUAL
# ============================================================

fig, ax = plt.subplots(figsize=(7, 7))

alcance = 1.15 * (l1 + l2)
ax.set_xlim(-alcance, alcance)
ax.set_ylim(-alcance, 0.35 * alcance)
ax.set_aspect("equal")
ax.grid(True, linestyle="--", alpha=0.3)

ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")

# Brazo completo: pivote -> masa 1 -> masa 2
linea, = ax.plot(
    [], [],
    "o-",
    lw=3,
    markersize=9
)

# Trayectoria de la segunda masa
trayectoria, = ax.plot(
    [], [],
    lw=1.2,
    alpha=0.45
)

texto_tiempo = ax.text(
    0.02, 0.96, "",
    transform=ax.transAxes,
    va="top"
)


def init():
    linea.set_data([], [])
    trayectoria.set_data([], [])
    texto_tiempo.set_text("")
    return linea, trayectoria, texto_tiempo


# Mostrar aproximadamente 30-35 FPS, aunque la integración use dt pequeño
salto_frames = max(1, int(round(0.03 / dt)))
indices_animacion = range(0, len(tiempos), salto_frames)


def update(i):
    # Posición de los tres puntos: pivote, masa 1 y masa 2
    linea.set_data(
        [0.0, x1[i], x2[i]],
        [0.0, y1[i], y2[i]]
    )

    # Cola reciente de la segunda masa
    segundos_cola = 2.0
    puntos_cola = max(1, int(segundos_cola / dt))
    inicio = max(0, i - puntos_cola)

    trayectoria.set_data(
        x2[inicio:i + 1],
        y2[inicio:i + 1]
    )

    ang1 = np.degrees(theta1[i])
    ang2 = np.degrees(theta2[i])

    texto_tiempo.set_text(
        f"t = {tiempos[i]:.2f} s\n"
        f"θ1 = {ang1:.1f}°   θ2 = {ang2:.1f}°\n"
        f"E = {energia[i]:.3f} J"
    )

    ax.set_title("Péndulo doble")

    return linea, trayectoria, texto_tiempo


ani = FuncAnimation(
    fig,
    update,
    frames=indices_animacion,
    init_func=init,
    blit=True,
    interval=30,
    repeat=False
)

plt.show()
