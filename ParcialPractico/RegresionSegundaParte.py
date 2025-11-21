import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import json

# --- Datos (temperaturas altas) ---
x = np.array([88.0,90.0,91.9,94.0,96.0,98.0,100.0,102.1,104.1,106.1,108.1,110.2,112.0,114.0,116.1,118.1,119.7], dtype=float)
y = np.array([4.903125,4.9111875,4.9149375,4.9190625,4.9231875,4.927875,4.93275,4.9365,4.9381875,4.939875,4.94625,4.9464375,4.9505625,4.951125,4.954125,4.9550625,4.959], dtype=float)

xmin = x.min()

# --- Modelos no polinómicos ---
def exp_saturating_shift(x, a, b, c):
    return c + a*(1 - np.exp(-b*(x - xmin)))

def asymp_exp_shift(x, a, b, c):
    return c + a * np.exp(-b*(x - xmin))

def michaelis_shift(x, vmax, km, c):
    xs = x - xmin + 1e-12
    return c + vmax * xs / (km + xs)

def logistic4(x, A, K, B, M):
    return A + (K - A) / (1 + np.exp(-B*(x - M)))

def gompertz_shift(x, a, b, c):
    return a * np.exp(-b * np.exp(-c * (x - xmin)))

def rational2(x, p, q, r):
    return (p*x + q) / (x + r)

def inverse_power(x, a, b, c):
    return c + a * (x**(-b))

candidates = {
    "Exp_saturating_shift": (exp_saturating_shift, [0.06, 0.02, 4.9]),
    "Asymp_exp_shift": (asymp_exp_shift, [0.06, 0.02, 4.9]),
    "Michaelis_shift": (michaelis_shift, [0.5, 10.0, 4.9]),
    "Logistic4": (logistic4, [4.9, 5.0, 0.05, 105.0]),
    "Gompertz_shift": (gompertz_shift, [4.9, 1.0, 0.01]),
    "Rational": (rational2, [0.01, 4.8, 80.0]),
    "Inverse_power": (inverse_power, [1.0, 1.0, 4.9])
}

def r2_score(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred)**2)
    ss_tot = np.sum((y_true - np.mean(y_true))**2)
    return 1 - ss_res/ss_tot

# --- Ajuste ---
results = []
for name, (func, p0) in candidates.items():
    try:
        if name == "Logistic4":
            popt, pcov = curve_fit(func, x, y, p0=p0,
                                   bounds=([-10, -10, 1e-8, x.min()], [10, 10, 1.0, x.max()]),
                                   maxfev=200000)
        elif name in ("Exp_saturating_shift", "Asymp_exp_shift"):
            popt, pcov = curve_fit(func, x, y, p0=p0, bounds=([-10, 1e-8, -10], [10, 1.0, 10]), maxfev=200000)
        else:
            popt, pcov = curve_fit(func, x, y, p0=p0, maxfev=200000)
        y_pred = func(x, *popt)
        r2 = r2_score(y, y_pred)
        results.append((name, float(r2), popt))
    except Exception:
        results.append((name, None, None))

# Ordenar por R2
results_sorted = sorted(results, key=lambda r: (r[1] is not None, r[1]), reverse=True)
rows = [{"Modelo": r[0], "R2": (r[1] if r[1] is not None else None), "Params": (np.round(r[2],8).tolist() if r[2] is not None else None)} for r in results_sorted]
df = pd.DataFrame(rows)
print(df.to_string(index=False))
df.to_csv("resumen_no_polinomico_dataset4.csv", index=False)

# Elegir mejor no polinómico
best = next((r for r in results_sorted if r[1] is not None), None)
if best is None:
    raise RuntimeError("No se obtuvo ningún ajuste válido.")
best_name, best_r2, best_popt = best

# --- Nombres de parámetros y ecuación numérica legible ---
def params_named_dict(name, popt):
    if name == "Exp_saturating_shift":
        a,b,c = popt; return {"a":float(a), "b":float(b), "c":float(c)}
    if name == "Asymp_exp_shift":
        a,b,c = popt; return {"a":float(a), "b":float(b), "c":float(c)}
    if name == "Michaelis_shift":
        vmax,km,c = popt; return {"vmax":float(vmax),"km":float(km),"c":float(c)}
    if name == "Logistic4":
        A,K,B,M = popt; return {"A":float(A),"K":float(K),"B":float(B),"M":float(M)}
    if name == "Gompertz_shift":
        a,b,c = popt; return {"a":float(a),"b":float(b),"c":float(c)}
    if name == "Rational":
        p,q,r = popt; return {"p":float(p),"q":float(q),"r":float(r)}
    if name == "Inverse_power":
        a,b,c = popt; return {"a":float(a),"b":float(b),"c":float(c)}
    return {f"p{i}": float(v) for i,v in enumerate(popt)}

best_params = params_named_dict(best_name, best_popt)

def equation_string(name, popt, xmin):
    # devuelve la ecuación numérica con coeficientes insertados (formato legible)
    if name == "Exp_saturating_shift":
        a,b,c = popt
        return f"y(T) = {c:.8f} + {a:.8f}*(1 - exp(-{b:.8f}*(T - {xmin:.2f})))"
    if name == "Asymp_exp_shift":
        a,b,c = popt
        return f"y(T) = {c:.8f} + {a:.8f}*exp(-{b:.8f}*(T - {xmin:.2f}))"
    if name == "Michaelis_shift":
        vmax,km,c = popt
        return f"y(T) = {c:.8f} + {vmax:.8f}*(T - {xmin:.2f})/({km:.8f} + (T - {xmin:.2f}))"
    if name == "Logistic4":
        A,K,B,M = popt
        return f"y(T) = {A:.8f} + ({K:.8f} - {A:.8f})/(1 + exp(-{B:.8f}*(T - {M:.8f})))"
    if name == "Gompertz_shift":
        a,b,c = popt
        return f"y(T) = {a:.8f} * exp(-{b:.8f} * exp(-{c:.8f}*(T - {xmin:.2f})))"
    if name == "Rational":
        p,q,r = popt
        return f"y(T) = ({p:.8f}*T + {q:.8f})/(T + {r:.8f})"
    if name == "Inverse_power":
        a,b,c = popt
        return f"y(T) = {c:.8f} + {a:.8f} * T^(-{b:.8f})"
    return "Ecuación no formateada"

# Imprimir y guardar ecuación + coeficientes
eq_str = equation_string(best_name, best_popt, xmin)
print("\n--- Mejor modelo no polinómico ---")
print("Modelo:", best_name)
print(f"R² = {best_r2:.10f}")
print("Ecuación (numérica):")
print(eq_str)
print("Coeficientes nombrados:")
for k,v in best_params.items():
    print(f"  {k} = {v:.8f}")

# Guardar a JSON y TXT
out = {"modelo": best_name, "R2": best_r2, "ecuacion": eq_str, "params": best_params}
with open("mejor_no_polinomico_dataset4_coeficientes.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=4, ensure_ascii=False)
with open("mejor_no_polinomico_dataset4_coeficientes.txt", "w", encoding="utf-8") as f:
    f.write("Mejor ajuste no polinómico\n")
    f.write(f"Modelo: {best_name}\n")
    f.write(f"R2: {best_r2:.10f}\n\n")
    f.write("Ecuación:\n")
    f.write(eq_str + "\n\n")
    f.write("Coeficientes:\n")
    for k,v in best_params.items():
        f.write(f"{k} = {v:.8f}\n")

print("\nArchivos guardados:")
print(" - resumen_no_polinomico_dataset4.csv")
print(" - mejor_no_polinomico_dataset4_coeficientes.json")
print(" - mejor_no_polinomico_dataset4_coeficientes.txt")

# Graficar mejor ajuste y mostrar ecuación en figura
xx = np.linspace(x.min(), x.max(), 400)
yy = candidates[best_name][0](xx, *best_popt)

plt.figure(figsize=(9.5,4.5))
plt.plot(x, y, 'o', label='Datos')
plt.plot(xx, yy, '-', label=f'Ajuste: {best_name}')
plt.xlabel("Temperatura")
plt.ylabel("Vc(m)")
plt.grid(True)
plt.legend(loc='lower right')
plt.title(f"Mejor no polinómico: {best_name} (R² = {best_r2:.8f})")
plt.text(0.02, 0.92, eq_str + f"\nR² = {best_r2:.8f}", transform=plt.gcf().transFigure,
         fontsize=9, va='top', bbox=dict(facecolor='white', alpha=0.9, edgecolor='black'))
plt.tight_layout(rect=[0,0,1,0.95])
plt.savefig("mejor_no_polinomico_dataset4.png", dpi=150, bbox_inches='tight')
plt.show()
print("Gráfica guardada en: mejor_no_polinomico_dataset4.png")
