import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Datos (Lux, Vcm)
x = np.array([400,800,1200,1600,2000,2400,2800,3200,3600,4000,4400,4800,5200,5600,6000], dtype=float)
y = np.array([0.0019,2.111,2.994,3.487,3.823,4.05,4.23,4.39,4.51,4.61,4.69,4.76,4.83,4.88,4.93], dtype=float)

def model(x, A, B, C):
    return A - B / (x + C)

# Intentar importar curve_fit
try:
    from scipy.optimize import curve_fit
except ImportError:
    raise ImportError(
        "scipy no está instalado. Instálalo con:\n"
        "  py -m pip install scipy\n"
        "o\n"
        "  python -m pip install scipy\n"
    )

# Estimaciones iniciales (puedes ajustar si falla la convergencia)
p0 = [5.0, 1000.0, 0.0]

# Ajuste
popt, pcov = curve_fit(model, x, y, p0=p0, maxfev=200000)
perr = np.sqrt(np.diag(pcov))
A_fit, B_fit, C_fit = popt

# Predicciones y R^2
y_pred = model(x, *popt)
ss_res = np.sum((y - y_pred)**2)
ss_tot = np.sum((y - np.mean(y))**2)
R2 = 1 - ss_res/ss_tot

# Mostrar resultados en consola
print("Ajuste: V(x) = A - B / (x + C)")
print(f"A = {A_fit:.6f} ± {perr[0]:.6f}")
print(f"B = {B_fit:.6f} ± {perr[1]:.6f}")
print(f"C = {C_fit:.6f} ± {perr[2]:.6f}")
print(f"R^2 = {R2:.6f}\n")

# Crear DataFrame y guardar CSV
df = pd.DataFrame({"Lux": x, "V_obs": y, "V_pred": y_pred, "residual": y - y_pred})
csv_name = "resultados_ajuste_lux_inv.csv"
df.to_csv(csv_name, index=False)
print(f"Guardado: {csv_name}")

# Gráfica con solo R^2 mostrado
x_smooth = np.linspace(x.min(), x.max(), 400)
y_smooth = model(x_smooth, *popt)

plt.figure(figsize=(8,5))
plt.plot(x, y, 'o', label='Datos observados')
plt.plot(x_smooth, y_smooth, '-', label='Ajuste: A - B/(x + C)')
plt.xlabel('Lux')
plt.ylabel('Vcm')
plt.title('Ajuste: V(x) = A - B/(x + C)')
plt.grid(True)

# Mostrar solo R^2 (texto en la esquina superior izquierda)
r2_text = rf"$R^2 = {R2:.6f}$"
plt.gca().text(0.05, 0.95, r2_text, transform=plt.gca().transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white', alpha=0.8))

plt.legend()
plt.tight_layout()
plt.show()
