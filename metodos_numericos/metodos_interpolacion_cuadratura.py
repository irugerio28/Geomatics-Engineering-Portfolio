"""
Autor: Ivan Rugerio
Título: Biblioteca de Interpolación, Cuadratura e Integración Numérica
Descripción: Módulo que implementa algoritmos fundamentales de análisis numérico:
             - Interpolación: Newton (Diferencias Divididas) y Ajustes Polinómicos.
             - Integración: Regla del Trapecio, Simpson Mixto y Cuadratura Gaussiana (n=2, n=3).
             - Diferenciación: Extrapolación de Richardson.
Tecnologías: Python, NumPy, SciPy, Matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from typing import Callable, Tuple

# =============================================================================
# SECCIÓN 1: HERRAMIENTAS DE INTERPOLACIÓN
# =============================================================================

def diferencias_divididas(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    """
    Calcula la tabla de diferencias divididas de Newton.
    Retorna la tabla completa, donde la primera fila contiene los coeficientes b_i.
    """
    n = len(xs)
    tabla = np.zeros((n, n))
    tabla[:, 0] = ys
    for j in range(1, n):
        for i in range(n - j):
            tabla[i, j] = (tabla[i + 1, j - 1] - tabla[i, j - 1]) / (xs[i + j] - xs[i])
    return tabla

def evaluar_newton(xs: np.ndarray, coeficientes: np.ndarray, val: float) -> float:
    """Evalúa el polinomio de interpolación de Newton en un punto 'val'."""
    n = len(xs)
    res = coeficientes[0]
    for i in range(1, n):
        termino = coeficientes[i]
        for j in range(i):
            termino *= (val - xs[j])
        res += termino
    return res

def ajuste_polinomico(x: np.ndarray, y: np.ndarray, grado: int, x_eval: float = None) -> Tuple[np.poly1d, float]:
    """
    Realiza un ajuste polinómico (mínimos cuadrados) de grado especificado.
    Retorna el polinomio objeto y el valor evaluado (si se solicitó).
    """
    coef = np.polyfit(x, y, grado)
    polinomio = np.poly1d(coef)
    valor = polinomio(x_eval) if x_eval is not None else None
    return polinomio, valor

# =============================================================================
# SECCIÓN 2: HERRAMIENTAS DE INTEGRACIÓN (CUADRATURA)
# =============================================================================

def integrar_gauss_legendre(func: Callable, a: float, b: float, n: int) -> float:
    """
    Calcula la integral definida usando Cuadratura Gaussiana de n puntos (n=2 o n=3).
    Realiza automáticamente el cambio de variable de [a, b] a [-1, 1].
    """
    # Pesos y nodos definidos para [-1, 1]
    if n == 2:
        w = np.array([1.0, 1.0])
        x_nodes = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
    elif n == 3:
        w = np.array([5/9, 8/9, 5/9])
        x_nodes = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    else:
        raise ValueError("Solo se soportan n=2 y n=3 en esta implementación.")

    # Cambio de variable: t = ((b-a)/2) * x + (b+a)/2
    # dx = ((b-a)/2) dt
    factor_dx = (b - a) / 2.0
    suma = 0.0
    
    for i in range(n):
        t = factor_dx * x_nodes[i] + (b + a) / 2.0
        suma += w[i] * func(t)
        
    return factor_dx * suma

def regla_trapecio_multiple(func: Callable, a: float, b: float, n: int) -> float:
    """Calcula la integral usando la regla del trapecio compuesta."""
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = func(x)
    return (h / 2) * (y[0] + 2 * np.sum(y[1:n]) + y[n])

def regla_simpson_mixta_n5(func: Callable, a: float, b: float) -> float:
    """
    Calcula la integral para n=5 segmentos combinando:
    - Simpson 1/3 para los primeros 2 segmentos.
    - Simpson 3/8 para los últimos 3 segmentos.
    """
    n_total = 5
    h = (b - a) / n_total
    
    # Generar puntos: x0, x1, x2, x3, x4, x5
    x = [a + i*h for i in range(n_total + 1)]
    
    # Parte 1: Simpson 1/3 (Intervalo [x0, x2])
    I_13 = (h / 3) * (func(x[0]) + 4 * func(x[1]) + func(x[2]))
    
    # Parte 2: Simpson 3/8 (Intervalo [x2, x5])
    I_38 = (3 * h / 8) * (func(x[2]) + 3 * func(x[3]) + 3 * func(x[4]) + func(x[5]))
    
    return I_13 + I_38

# =============================================================================
# SECCIÓN 3: HERRAMIENTAS DE DIFERENCIACIÓN
# =============================================================================

def diferencia_central(func: Callable, x: float, h: float) -> float:
    """Derivada numérica usando diferencia central O(h^2)."""
    return (func(x + h) - func(x - h)) / (2 * h)

def extrapolacion_richardson(D_h1: float, D_h2: float) -> float:
    """
    Mejora la aproximación de la derivada usando Extrapolación de Richardson O(h^4).
    Asume que h2 = h1 / 2.
    Formula: D_rich = (4 * D(h/2) - D(h)) / 3
    """
    return (4 * D_h2 - D_h1) / 3

# =============================================================================
# SECCIÓN 4: DEMOSTRACIONES (CASOS DE ESTUDIO)
# =============================================================================

def demo_interpolacion():
    print("\n=== DEMO: Interpolación y Ajuste de Curvas ===")
    
    # Caso Newton
    x1 = np.array([-3, -1, 1, 3, 5], dtype=float)
    y1 = np.array([-51, -11, -11, -3, 61], dtype=float)
    D1 = diferencias_divididas(x1, y1)
    y_newton = evaluar_newton(x1, D1[0, :], 2)
    print(f"1. Newton (x=2): {y_newton:.4f}")

    # Caso Ajuste Parabólico (Trayectoria)
    x2 = np.arange(0, 90, 5, dtype=float)
    y2 = np.array([0,8.1,15.1,21.1,25.9,29.7,32.3,34.4,33.7,33.7,32.1,
                  29.9,25.4,20.5,14.4,7.3,-0.96,-10.3])
    poly2, val2 = ajuste_polinomico(x2, y2, 2, x_eval=42.3)
    print(f"2. Ajuste Parabólico (x=42.3): {val2:.4f}")
    
    # Visualización rápida
    t = np.linspace(0, 85, 100)
    plt.figure("Ajuste Parabólico")
    plt.plot(x2, y2, 'o', label='Datos')
    plt.plot(t, poly2(t), '-', label='Ajuste Cuadrático')
    plt.legend()
    plt.grid(True)
    print("   (Gráfica generada)")

def demo_cuadratura():
    print("\n=== DEMO: Cuadratura Gaussiana ===")
    
    # Integral 1: ∫[1, 5] (2t * e^(t²)) / (6t³ - 4) dt
    def g_a(t): 
        den = 6 * t**3 - 4
        return (2 * t * np.exp(t**2)) / den if den != 0 else 0

    res_n2 = integrar_gauss_legendre(g_a, 1, 5, n=2)
    res_n3 = integrar_gauss_legendre(g_a, 1, 5, n=3)
    print(f"Integral A (n=2): {res_n2:.4f}")
    print(f"Integral A (n=3): {res_n3:.4f}")

    # Integral 2: ∫[-1, 1] cos(t + π) dt
    # Valor exacto: -2*sin(1) ≈ -1.6829
    def g_c(t): return np.cos(t + np.pi)
    res_c = integrar_gauss_legendre(g_c, -1, 1, n=3)
    print(f"Integral B (n=3): {res_c:.4f} (Exacto: {-1.6829419:.4f})")

def demo_integracion_numerica():
    print("\n=== DEMO: Integración (Trapecio vs Simpson) ===")
    
    # Integral: ∫[1, 5] ... (Misma función g_a anterior)
    def func_a(x): return (2 * x * np.exp(x**2)) / (6 * x**2 - 4)
    
    valor_exacto, _ = quad(func_a, 1, 5)
    
    trap = regla_trapecio_multiple(func_a, 1, 5, n=6)
    simp = regla_simpson_mixta_n5(func_a, 1, 5)
    
    print(f"Valor SciPy (Ref): {valor_exacto:.6f}")
    print(f"Trapecio (n=6):    {trap:.6f} (Err: {abs((valor_exacto-trap)/valor_exacto)*100:.2f}%)")
    print(f"Simpson (n=5):     {simp:.6f} (Err: {abs((valor_exacto-simp)/valor_exacto)*100:.2f}%)")

def demo_richardson():
    print("\n=== DEMO: Extrapolación de Richardson ===")
    
    def f(x): return np.sin(x) - np.exp(-x)
    def f_prime_real(x): return np.cos(x) + np.exp(-x)
    
    x_val = 0.5
    h1, h2 = 0.5, 0.25
    
    D1 = diferencia_central(f, x_val, h1)
    D2 = diferencia_central(f, x_val, h2)
    D_rich = extrapolacion_richardson(D1, D2)
    real = f_prime_real(x_val)
    
    print(f"Diferencia Central (h={h1}): {D1:.5f}")
    print(f"Diferencia Central (h={h2}): {D2:.5f}")
    print(f"Richardson Extrapolado:      {D_rich:.5f}")
    print(f"Valor Real:                  {real:.5f}")

# =============================================================================
# EJECUCIÓN PRINCIPAL
# =============================================================================
if __name__ == "__main__":
    demo_interpolacion()
    demo_cuadratura()
    demo_integracion_numerica()
    demo_richardson()
    
    print("\n--- Cerrar ventanas de gráficos para finalizar ---")
    plt.show()