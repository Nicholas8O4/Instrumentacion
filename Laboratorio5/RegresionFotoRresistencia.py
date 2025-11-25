import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Datos
m = np.array([200,300,400,500,600,700,800,900,1000,2000,3000,4000,6000,8000], dtype=float)
V = np.array([0.2162,2.015,2.833,3.217,3.434,3.609,3.785,3.961,4.02,4.35,4.49,4.61,4.68,4.82], dtype=float)

def model(m, A, B, C):
    return A - B / (m + C)

try:
    from scipy.optimize import curve_fit
except ImportError:
    raise ImportError(
        "scipy no está instalado. Instálalo con:\n"
        "  py -m pip install scipy\n"
        "o\n"
        "  python -m pip install scipy\n"
    )

# Estimaciones iniciales y límites útiles
A0 = V.max()
B0 = (A0 - V[0]) * (m[0] + 10)
C0 = -40
p0 = [A0, B0, C0]

min_m = m.min()
C_lower = -min_m + 1e-3
bounds_lower = [-100.0, -1e7, C_lower]
bounds_upper = [100.0, 1e7, 1e7]

# Ajuste
popt, pcov = curve_fit(model, m, V, p0=p0, bounds=(bounds_lower, bounds_upper), maxfev=20000)
A_fit, B_fit, C_fit = popt
perr = np.sqrt(np.diag(pcov))

# Predicciones y R^2
V_pred = model(m, *popt)
residuals = V - V_pred
ss_res = np.sum(residuals**2)
ss_tot = np.sum((V - np.mean(V))**2)
R2 = 1 - ss_res/ss_tot

# DataFrame con resultados
df = pd.DataFrame({
    "m_g": m,
    "V_obs": V,
    "V_pred": V_pred,
    "residual": residuals
})

# Mostrar resultados en consola
print("Parámetros ajustados:")
print(f"A = {A_fit:.6f} ± {perr[0]:.6f}")
print(f"B = {B_fit:.6f} ± {perr[1]:.6f}")
print(f"C = {C_fit:.6f} ± {perr[2]:.6f}")
print(f"R^2 = {R2:.6f}\n")

print("Tabla resumen:")
print(df.to_string(index=False))

# Guardar CSV
csv_name = "resultados_ajuste.csv"
df.to_csv(csv_name, index=False)
print(f"\nGuardado: {csv_name}")

# Gráfica (con R^2 mostrado en la gráfica)
m_smooth = np.linspace(m.min(), m.max(), 400)
V_smooth = model(m_smooth, *popt)

plt.figure(figsize=(8,5))
plt.plot(m, V, marker='o', linestyle='None', label='Datos observados')
plt.plot(m_smooth, V_smooth, linestyle='-', label='Ajuste: A - B/(m + C)')
plt.xlabel('m (g)')
plt.ylabel('V (Vcm)')
plt.title('Ajuste no lineal: V(m) = A - B/(m + C)')
plt.grid(True)

# Mostrar solo R^2 en la gráfica (sin mostrar RMSE)
r2_text = rf"$R^2 = {R2:.6f}$"
# Ubicación del cuadro: ajusta 'xy' si se superpone (ej: (0.05, 0.95) = esquina superior izquierda)
plt.gca().text(0.05, 0.95, r2_text, transform=plt.gca().transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white', alpha=0.8))

plt.legend()
plt.tight_layout()
plt.show()
