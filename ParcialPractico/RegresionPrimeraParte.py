import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import json
from caas_jupyter_tools import display_dataframe_to_user

x = np.array([-7.9, -5.9, -4.0, -1.8, -0.1, 1.8, 4.2, 6.6, 7.8, 10.1, 12.2, 14.1, 16.3, 17.9, 19.8, 22.4,
              24.4, 26.2, 28.0, 30.2, 31.8, 33.9, 36.0, 37.9, 40.0, 42.0, 44.0, 45.9, 47.9, 49.9, 51.9,
              54.0, 56.0, 58.0, 60.0, 62.3, 64.1, 66.0, 68.0, 70.0, 72.0, 73.9, 76.1, 77.9, 80.0, 82.1,
              84.1, 86.1], dtype=float)

y = np.array([0.50625,0.7365,1.0,1.536,1.78575,2.1013125,2.4615,2.7995625,2.9720625,3.1798125,
              3.33391875,3.48825,3.624375,3.723375,3.81,3.952875,4.045875,4.1131875,4.1825625,
              4.257375,4.3111875,4.3674375,4.4176875,4.46025,4.505625,4.5440625,4.576125,4.6048125,
              4.639125,4.6644375,4.6876875,4.7165625,4.73775,4.756125,4.7715,4.785,4.8,4.809375,
              4.8271875,4.835625,4.8444375,4.8556875,4.86675,4.8718125,4.88025,4.8885,4.8924375,
              4.8969375], dtype=float)

# Comprueba tamaños
assert x.size == y.size

# Modelos
xmin = x.min()
def r2_score(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred)**2)
    ss_tot = np.sum((y_true - np.mean(y_true))**2)
    return 1 - ss_res/ss_tot

def linear(x, a, b): return a*x + b
def poly2(x, a, b, c): return a*x**2 + b*x + c
def poly3(x, a, b, c, d): return a*x**3 + b*x**2 + c*x + d
def exp_sat_shift(x, a, b, c): return c + a*(1 - np.exp(-b*(x - xmin)))
def asymp_exp_shift(x, a, b, c): return c + a * np.exp(-b*(x - xmin))
def michaelis_shift(x, vmax, km, c):
    xs = x - xmin + 1e-8
    return c + vmax * xs / (km + xs)
def logistic4(x, A, K, B, M): return A + (K - A) / (1 + np.exp(-B*(x - M)))
def gompertz_shift(x, a, b, c): return a * np.exp(-b * np.exp(-c * (x - xmin)))

candidates = {
    "Linear": (linear, [0.01, 1.0]),
    "Poly2": (poly2, [1e-4, 1e-2, 1.0]),
    "Poly3": (poly3, [1e-6, 1e-3, 1e-2, 1.0]),
    "Exp_saturating_shift": (exp_sat_shift, [4.8, 0.05, 0.15]),
    "Asymp_exp_shift": (asymp_exp_shift, [4.8, 0.05, 0.15]),
    "Michaelis_shift": (michaelis_shift, [5.0, 10.0, 0.15]),
    "Logistic4": (logistic4, [0.15, 5.0, 0.2, 10.0]),
    "Gompertz_shift": (gompertz_shift, [5.0, 1.0, 0.02])
}

# Ajuste
results = []
from scipy.optimize import curve_fit
for name, (func, p0) in candidates.items():
    try:
        if name == "Logistic4":
            popt, pcov = curve_fit(func, x, y, p0=p0,
                                   bounds=([-np.inf, -np.inf, 1e-8, np.min(x)],
                                           [np.inf, np.inf, np.inf, np.max(x)]),
                                   maxfev=400000)
        elif name in ("Exp_saturating_shift","Asymp_exp_shift","Michaelis_shift"):
            popt, pcov = curve_fit(func, x, y, p0=p0,
                                   bounds=([-np.inf, 1e-8, -np.inf],[np.inf, np.inf, np.inf]),
                                   maxfev=300000)
        else:
            popt, pcov = curve_fit(func, x, y, p0=p0, maxfev=200000)
        y_pred = func(x, *popt)
        r2 = r2_score(y, y_pred)
        results.append((name, float(r2), popt))
    except Exception as e:
        results.append((name, None, None))

# Ordenar
results_sorted = sorted(results, key=lambda r: (r[1] is not None, r[1]), reverse=True)
rows = [{"Modelo": r[0], "R2": (r[1] if r[1] is not None else None), "Params": (np.round(r[2],6).tolist() if r[2] is not None else None)} for r in results_sorted]
res_df = pd.DataFrame(rows)
display_dataframe_to_user("Resultados de ajustes (ordenados por R²)", res_df)

# Guardar CSV y seleccionar mejor modelo
res_df.to_csv("/mnt/data/resumen_ajustes_dataset3.csv", index=False)
best = next((r for r in results_sorted if r[1] is not None), None)
best_name, best_r2, best_popt = best

# Nombrar parámetros
def params_named_dict(name, popt):
    if name == "Gompertz_shift":
        a, b, c = popt; return {"a": float(a), "b": float(b), "c": float(c)}
    if name == "Logistic4":
        A, K, B, M = popt; return {"A": float(A), "K": float(K), "B": float(B), "M": float(M)}
    if name == "Poly3":
        a, b, c, d = popt; return {"a3": float(a), "a2": float(b), "a1": float(c), "a0": float(d)}
    if name == "Poly2":
        a, b, c = popt; return {"a2": float(a), "a1": float(b), "a0": float(c)}
    if name == "Exp_saturating_shift":
        a, b, c = popt; return {"a": float(a), "b": float(b), "c": float(c)}
    if name == "Asymp_exp_shift":
        a, b, c = popt; return {"a": float(a), "b": float(b), "c": float(c)}
    if name == "Michaelis_shift":
        vmax, km, c = popt; return {"vmax": float(vmax), "km": float(km), "c": float(c)}
    if name == "Linear":
        a, b = popt; return {"a": float(a), "b": float(b)}
    return {"p"+str(i): float(v) for i, v in enumerate(popt)}

best_params_dict = params_named_dict(best_name, best_popt)

# Guardar coeficientes JSON/TXT
out_json = {"modelo": best_name, "R2": best_r2, "params": best_params_dict}
with open("/mnt/data/mejor_ajuste_dataset3_coeficientes.json", "w", encoding="utf-8") as f:
    json.dump(out_json, f, indent=4, ensure_ascii=False)
with open("/mnt/data/mejor_ajuste_dataset3_coeficientes.txt", "w", encoding="utf-8") as f:
    f.write(f"Modelo: {best_name}\nR2: {best_r2:.6f}\n")
    for k,v in best_params_dict.items():
        f.write(f"{k} = {v:.6f}\n")

# Mostrar en consola (a la salida visible del notebook)
print("Mejor modelo:", best_name)
print("R^2:", best_r2)
print("Parámetros nombrados:")
for k,v in best_params_dict.items():
    print(f"  {k} = {v:.6f}")

# Graficar
def format_equation(name, popt, xmin):
    if name == "Gompertz_shift":
        a,b,c = popt
        return rf"$y(T)={a:.6f}\,\exp\!\left(-{b:.6f}\,\exp\!\left(-{c:.6f}\,(T-{xmin:.2f})\right)\right)$"
    if name == "Logistic4":
        A,K,B,M = popt
        return rf"$y(T)={A:.6f}+ \dfrac{{{K:.6f}-{A:.6f}}}{{1+\exp(-{B:.6f}(T-{M:.6f}))}}$"
    if name == "Poly3":
        a,b,c,d = popt
        return rf"$y(T)={a:.6e}T^3+{b:.6e}T^2+{c:.6e}T+{d:.6e}$"
    if name == "Poly2":
        a,b,c = popt
        return rf"$y(T)={a:.6e}T^2+{b:.6e}T+{c:.6e}$"
    if name == "Exp_saturating_shift":
        a,b,c = popt
        return rf"$y(T)={c:.6f}+{a:.6f}\left(1-\exp(-{b:.6f}(T-{xmin:.2f}))\right)$"
    if name == "Asymp_exp_shift":
        a,b,c = popt
        return rf"$y(T)={c:.6f}+{a:.6f}\exp(-{b:.6f}(T-{xmin:.2f}))$"
    if name == "Michaelis_shift":
        vmax,km,c = popt
        return rf"$y(T)={c:.6f}+{vmax:.6f}\dfrac{{(T-{xmin:.2f})}}{{{km:.6f}+(T-{xmin:.2f})}}$"
    return "Ecuación no formateada"

equation = format_equation(best_name, best_popt, xmin)
r2_text = rf"$R^2 = {best_r2:.6f}$"
text_block = equation + "\n" + r2_text

fig, ax = plt.subplots(figsize=(10,5))
ax.plot(x, y, 'o', label='Datos')
xx = np.linspace(x.min(), x.max(), 600)
func = candidates[best_name][0]
yy = func(xx, *best_popt)
ax.plot(xx, yy, '-', label=f'Ajuste: {best_name}')
ax.set_xlabel("Temperatura")
ax.set_ylabel("Vc(m)")
ax.grid(True)
ax.legend(loc='lower right')

# Mostrar resumen de parámetros si la ecuación es muy larga
max_chars = 180
if len(text_block) > max_chars:
    short_lines = [f"{k}={v:.4f}" for k, v in best_params_dict.items()]
    short_text = ", ".join(short_lines)
    text_to_show = short_text + "\n" + r2_text
else:
    text_to_show = text_block

fig.text(0.02, 0.92, text_to_show, fontsize=10, va='top', ha='left',
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='black'))

plt.title(f"Mejor ajuste: {best_name} (R² = {best_r2:.6f})")
plt.tight_layout(rect=[0,0,1,0.95])
outpath = "/mnt/data/mejor_ajuste_dataset3.png"
plt.savefig(outpath, dpi=150, bbox_inches='tight', pad_inches=0.3)
plt.show()

print("Resumen CSV guardado en: /mnt/data/resumen_ajustes_dataset3.csv")
print("Coeficientes guardados en JSON/TXT en /mnt/data/")
print("Gráfica guardada en:", outpath)
