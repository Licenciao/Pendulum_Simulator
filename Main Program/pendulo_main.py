import numpy as np
import matplotlib
matplotlib.use("TkAgg") # Backend necesario para la animación en vivo
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


    
    # --- Interacción con el Usuario ---
print("--- Simulador de Péndulo Simple ---")
A_user = float(input("Amplitud inicial (grados, ej: 45) "))
l_user = float(input("Longitud del brazo (metros, ej: 1.0) "))
t_max = float(input("Tiempo de simulación (segundos, ej: 20) "))
g = float(input("Ingrese gravedad deseada (m/s**2, ej: 9.81) "))
m = float(input("Ingrese masa deseada (kg, ej: 1.0) "))
b = float(input("Ingrese coeficiente de friccion (kg, ej: 1.0)"))
# --- Parámetros Físicos ---
A = np.radians(A_user) 
l = l_user
w = np.sqrt(g/l)

# --- Funciones de Cálculo ---
def f(A, w, t):
    return A * np.cos(w * t)


def f_con_friccion(A, w, t, b):
    # La amplitud A se multiplica por un exponencial negativo
    return A * np.exp(-b * t) * np.cos(w * t)

def theta_to_xy(theta, l):
    x = l * np.sin(theta)
    y = -l * np.cos(theta)
    return x, y

# --- Configuración Visual ---
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-l - 5, l + 5)
ax.set_ylim(-l - 5, l/2)
ax.set_aspect('equal')
ax.grid(True, linestyle='--', alpha=0.3)

# 'o-' dibuja el círculo (masa) y la línea (brazo) en un solo objeto
linea, = ax.plot([], [], 'o-', lw=3, color='royalblue', 
                 markersize=10, markerfacecolor='orange', markeredgecolor='black')

def init():
    linea.set_data([], [])
    return linea,

def update(t):
    theta = f_con_friccion(A, w, t, b)
    
    x, y = theta_to_xy(theta, l)
    
    # --- Cálculo de Fuerzas ---
    # Velocidad angular: v_ang = -A * w * sin(w * t)
    v_lineal = -l * A * w * np.sin(w * t) 
    
    fg = m * g
    ft = -m * g * np.sin(theta)
    fc = m * (v_lineal**2) / l
    tension = m * g * np.cos(theta) + fc
    
    # Aplicar al gráfico
    linea.set_data([0, x], [0, y])
    
    # Ejemplo de cómo mostrarlo en el título dinámico
    ax.set_title(f"Tensión: {tension:.2f}N | F_tangencial: {ft:.2f}N | Gravedad: {g:.2f}N")
    
    return linea,

# Creamos el rango de tiempo
t_frames = np.arange(0, t_max, 0.05)

# Lanzar la animación
ani = FuncAnimation(fig, update, frames=t_frames, 
                    init_func=init, blit=True, interval=30, repeat=False)



plt.show()
