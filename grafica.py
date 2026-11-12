import pandas as pd
import matplotlib.pyplot as plt

d = pd.read_csv("kappa.csv")
d = d[d.kappa < 100]   # recorta la zona de la singularidad
plt.plot(d.q2, d.kappa)
plt.xlabel("q2 (rad)")
plt.ylabel("kappa(J)")
plt.title("Condicionamiento del Jacobiano del brazo 2R")
plt.grid(True)
plt.savefig("kappa.png", dpi=200)
plt.show()