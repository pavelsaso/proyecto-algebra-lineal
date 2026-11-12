import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

L1, L2 = 1.0, 0.8
d = pd.read_csv("trayectoria.csv")

# ---------- Figura estatica con 3 graficas ----------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

ax[0].plot(d.x_obj, d.y_obj, "k--", label="objetivo")
ax[0].plot(d.x, d.y, "r.", markersize=3, label="alcanzado")
ax[0].set_aspect("equal")
ax[0].set_xlabel("x"); ax[0].set_ylabel("y")
ax[0].set_title("Trayectoria del efector")
ax[0].legend(); ax[0].grid(True)

ax[1].plot(d.t, d.q1, label="q1")
ax[1].plot(d.t, d.q2, label="q2")
ax[1].set_xlabel("parametro t (rad)"); ax[1].set_ylabel("angulo (rad)")
ax[1].set_title("Angulos de las articulaciones")
ax[1].legend(); ax[1].grid(True)

ax[2].plot(d.t, d.kappa)
ax[2].set_xlabel("parametro t (rad)"); ax[2].set_ylabel("kappa(J)")
ax[2].set_title("Condicionamiento a lo largo del recorrido")
ax[2].grid(True)

fig.tight_layout()
fig.savefig("trayectoria.png", dpi=200)

# ---------- Animacion ----------
q1 = d.q1.values
ex = L1 * np.cos(q1)
ey = L1 * np.sin(q1)
hx = d.x.values
hy = d.y.values

fig2, a = plt.subplots(figsize=(6, 6))
a.set_xlim(-2, 2); a.set_ylim(-2, 2)
a.set_aspect("equal"); a.grid(True)
a.set_title("Brazo 2R siguiendo un circulo")
a.plot(d.x_obj, d.y_obj, "k--", lw=1)
brazo, = a.plot([], [], "o-", lw=4, color="tab:blue")
rastro, = a.plot([], [], "r-", lw=1.5)

def update(i):
    brazo.set_data([0, ex[i], hx[i]], [0, ey[i], hy[i]])
    rastro.set_data(hx[:i+1], hy[:i+1])
    return brazo, rastro

ani = FuncAnimation(fig2, update, frames=len(d), interval=30, blit=True)
ani.save("brazo.gif", writer=PillowWriter(fps=30))
plt.show()