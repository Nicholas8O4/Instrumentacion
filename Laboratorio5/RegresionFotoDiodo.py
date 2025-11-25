import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Datos experimentales
x = np.array([105, 201, 306, 403, 500, 597, 705, 801, 906, 1003, 1100, 1206, 1305, 1400, 1505], dtype=float)
y = np.array([0.553, 0.993, 1.438, 1.813, 2.172, 2.512, 2.866, 3.17, 3.488, 3.766, 4.02, 4.3, 4.56, 4.79, 5.04], dtype=float)

# Modelos
def modelo_inverso(x, A, B, C):
    return A - B / (x + C)

def modelo_exponencial(x, A, B, k):
    return A - B * np.exp(-k * x)

# Ajuste de cada modelo
popt_inv, _ = curve_fit(modelo_inverso, x, y, p0=[5, 100, 10])
popt_exp, _ = curve_fit(modelo_exponencial, x, y, p0=[5, 5, 0.001])

# Predicciones
y_inv = modelo_inverso(x, *popt_inv)
y_exp = modelo_exponencial(x, *popt_exp)

# Calcular R²
def r2(y, y_pred):
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    return 1 - ss_res/ss_tot

R2_inv = r2(y, y_inv)
R2_exp = r2(y, y_exp)

# Mostrar resultados
print("Modelo inverso: V = A - B/(x + C)")
print(f"A={popt_inv[0]:.6f}, B={popt_inv[1]:.6f}, C={popt_inv[2]:.6f}, R²={R2_inv:.6f}\n")

print("Modelo exponencial: V = A - B·exp(-k·x)")
print(f"A={popt_exp[0]:.6f}, B={popt_exp[1]:.6f}, k={popt_exp[2]:.6f}, R²={R2_exp:.6f}\n")

# Elegir mejor modelo
if R2_inv > R2_exp:
    mejor_modelo = modelo_inverso
    popt_mejor = popt_inv
    R2_mejor = R2_inv
    etiqueta = 'A - B/(x + C)'
else:
    mejor_modelo = modelo_exponencial
    popt_mejor = popt_exp
    R2_mejor = R2_exp
    etiqueta = 'A - B·exp(-k·x)'

# Graficar
x_smooth = np.linspace(min(x), max(x), 400)
y_smooth = mejor_modelo(x_smooth, *popt_mejor)

plt.figure(figsize=(8,5))
plt.plot(x, y, 'o', label='Datos experimentales')
plt.plot(x_smooth, y_smooth, '-', label=f'Ajuste {etiqueta}')
plt.xlabel('Lux')
plt.ylabel('Vcm')
plt.title('Ajuste del modelo a los datos')
plt.grid(True)

# Mostrar solo R² en el gráfico
plt.gca().text(0.05, 0.95, f'$R^2 = {R2_mejor:.6f}$',
               transform=plt.gca().transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))

plt.legend()
plt.tight_layout()
plt.show()
