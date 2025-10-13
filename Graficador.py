from __future__ import annotations
import argparse
import os
import re
from itertools import cycle
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd

# Expresión regular para números (coma o punto decimal, exponentes)
PATRON_NUMERO = re.compile(r"[-+]?\d*[\.,]?\d+(?:[eE][-+]?\d+)?")


def convertir_a_flotante(token: str) -> Optional[float]:
    try:
        return float(token.replace(',', '.'))
    except Exception:
        return None


def leer_bloques(ruta: str) -> List[Tuple[pd.DataFrame, Optional[str]]]:
    if not os.path.isfile(ruta):
        raise FileNotFoundError(ruta)

    texto = open(ruta, 'r', encoding='utf-8').read()
    bloques = [b.strip() for b in re.split(r"\n\s*\n+", texto) if b.strip()]

    series: List[Tuple[pd.DataFrame, Optional[str]]] = []
    for bloque in bloques:
        lineas = [l.strip() for l in bloque.splitlines() if l.strip()]
        if not lineas:
            continue
        encabezado = None
        if not PATRON_NUMERO.search(lineas[0]):
            encabezado = lineas[0]
            datos_lineas = lineas[1:]
        else:
            datos_lineas = lineas

        filas = []
        for i, linea in enumerate(datos_lineas):
            tokens = PATRON_NUMERO.findall(linea)
            valores = [convertir_a_flotante(t) for t in tokens]
            valores = [v for v in valores if v is not None]
            if len(valores) >= 2:
                x, y = float(valores[0]), float(valores[1])
            elif len(valores) == 1:
                x, y = float(i), float(valores[0])
            else:
                continue
            filas.append({'x': x, 'y': y})

        if not filas:
            continue
        datos = pd.DataFrame(filas).sort_values('x', kind='mergesort').reset_index(drop=True)
        series.append((datos, encabezado))
    return series


def inferir_etiquetas(series: List[Tuple[pd.DataFrame, Optional[str]]]) -> Tuple[str, str]:
    for _, encabezado in series:
        if not encabezado:
            continue
        coincidencia = re.search(r"\bvs\.?\b", encabezado, flags=re.IGNORECASE)
        if coincidencia:
            izquierda = encabezado[:coincidencia.start()].strip()
            derecha = encabezado[coincidencia.end():].strip()
            etiqueta_x = izquierda.split()[0] if izquierda else 'X'
            etiqueta_y = derecha.split()[0] if derecha else 'Y'
            return etiqueta_x, etiqueta_y
    return 'X', 'Y'


def graficar_series(series: List[Tuple[pd.DataFrame, Optional[str]]],
                    etiqueta_x: str = 'X', etiqueta_y: str = 'Y', titulo: Optional[str] = None,
                    guardar: Optional[str] = None, invertir_x: bool = False) -> None:
    if not series:
        raise ValueError('No hay series para graficar')

    figura, (ax_lineal, ax_log) = plt.subplots(1, 2, figsize=(14, 6))
    marcadores = cycle(['o', 's', '^', 'D', 'v', 'P'])
    colores = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

    # Escala lineal
    for i, (datos, encabezado) in enumerate(series, 1):
        if datos.empty:
            continue
        etiqueta = f"Serie {i} ({len(datos)})"
        if encabezado:
            etiqueta += f" — {encabezado}"
        ax_lineal.plot(datos['x'], datos['y'], marker=next(marcadores), linestyle='-', label=etiqueta, color=next(colores))

    ax_lineal.set_xlabel(etiqueta_x)
    ax_lineal.set_ylabel(etiqueta_y)
    ax_lineal.set_title('Escala lineal')
    ax_lineal.grid(True, linestyle=':', alpha=0.6)
    ax_lineal.legend()
    if invertir_x:
        ax_lineal.invert_xaxis()

    # Escala logarítmica
    marcadores = cycle(['o', 's', '^', 'D', 'v', 'P'])
    colores = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    hay_log = False
    for i, (datos, encabezado) in enumerate(series, 1):
        datos_pos = datos[(datos['x'] > 0) & (datos['y'] > 0)]
        if datos_pos.empty:
            continue
        etiqueta = f"Serie {i} ({len(datos_pos)})"
        if encabezado:
            etiqueta += f" — {encabezado}"
        ax_log.plot(datos_pos['x'], datos_pos['y'], marker=next(marcadores), linestyle='-', label=etiqueta, color=next(colores))
        hay_log = True

    ax_log.set_xlabel(etiqueta_x)
    ax_log.set_ylabel(etiqueta_y)
    ax_log.set_title('Escala log-log')
    ax_log.grid(True, which='both', linestyle=':', alpha=0.6)
    if hay_log:
        ax_log.set_xscale('log')
        ax_log.set_yscale('log')
        ax_log.legend()
        if invertir_x:
            ax_log.invert_xaxis()
    else:
        ax_log.text(0.5, 0.5, 'Sin datos válidos para log-log', ha='center', va='center', transform=ax_log.transAxes)

    if titulo:
        figura.suptitle(titulo)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if guardar:
        figura.savefig(guardar, dpi=300)
        print(f'Figura guardada en: {guardar}')
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description='Graficador simple usando pandas')
    parser.add_argument('archivo', nargs='?', default='datos.txt')
    parser.add_argument('-g', '--guardar')
    parser.add_argument('--invertir-x', action='store_true')
    parser.add_argument('--etiqueta-x')
    parser.add_argument('--etiqueta-y')
    argumentos = parser.parse_args()

    try:
        series = leer_bloques(argumentos.archivo)
    except FileNotFoundError:
        print('Archivo no encontrado:', argumentos.archivo)
        return

    if not series:
        print('No se encontraron series en el archivo.')
        return

    if argumentos.etiqueta_x or argumentos.etiqueta_y:
        etiqueta_x = argumentos.etiqueta_x or 'X'
        etiqueta_y = argumentos.etiqueta_y or 'Y'
    else:
        etiqueta_x, etiqueta_y = inferir_etiquetas(series)

    print(f"Ejes -> X: '{etiqueta_x}', Y: '{etiqueta_y}'")
    titulo = f"Gráfica ({len(series)} series) — {argumentos.archivo}"
    graficar_series(series, etiqueta_x=etiqueta_x, etiqueta_y=etiqueta_y, titulo=titulo, guardar=argumentos.guardar, invertir_x=argumentos.invertir_x)


if __name__ == '__main__':
    main()
