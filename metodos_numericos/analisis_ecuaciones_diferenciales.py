"""
Autor: Ivan Rugerio
Título: Solucionador Numérico de Ecuaciones Diferenciales (EDO)
Descripción: Biblioteca de métodos numéricos (Euler, RK4, Taylor, Milne, Diferencias Finitas) 
             para la resolución de EDOs de primer orden, orden superior y sistemas acoplados.
Tecnologías: Python, NumPy, Pandas
"""

import numpy as np
import pandas as pd
from typing import Callable

# =============================================================================
# SECCIÓN 1: MOTOR DE MÉTODOS NUMÉRICOS (SOLVERS)
# =============================================================================

def euler(f: Callable, x0: float, y0: float, h: float, x_end: float) -> pd.DataFrame:
    """Resuelve una EDO de 1er orden con el método de Euler."""
    x_values = [x0]
    y_values = [y0]
    
    x, y = x0, y0
    
    while x < x_end - 1e-9: # Corrección para punto flotante
        y = y + h * f(x, y)
        x = x + h
        x_values.append(x)
        y_values.append(y)
        
    return pd.DataFrame({'x': x_values, 'y': y_values})

def euler_gauss(f: Callable, x0: float, y0: float, h: float, x_end: float) -> pd.DataFrame:
    """Resuelve una EDO de 1er orden con el método de Euler-Gauss (Heun / RK2)."""
    x_values = [x0]
    y_values = [y0]
    
    x, y = x0, y0
    
    while x < x_end - 1e-9:
        y_pred = y + h * f(x, y)
        y = y + (h / 2) * (f(x, y) + f(x + h, y_pred))
        x = x + h
        x_values.append(x)
        y_values.append(y)
        
    return pd.DataFrame({'x': x_values, 'y': y_values})

def rk2_midpoint(f: Callable, x0: float, y0: float, h: float, x_end: float) -> pd.DataFrame:
    """Resuelve una EDO de 1er orden con el método RK2 (Punto Medio)."""
    x_values = [x0]
    y_values = [y0]
    
    x, y = x0, y0
    
    while x < x_end - 1e-9:
        k1 = f(x, y)
        k2 = f(x + h/2, y + (h/2) * k1)
        y = y + h * k2
        x = x + h
        x_values.append(x)
        y_values.append(y)
        
    return pd.DataFrame({'x': x_values, 'y': y_values})

def taylor_orden3(f_derivs: Callable, x0: float, y0: float, h: float, x_end: float) -> pd.DataFrame:
    """Resuelve una EDO de 1er orden con Series de Taylor de 3er grado."""
    x_values = [x0]
    y_values = [y0]
    
    x, y = x0, y0
    
    while x < x_end - 1e-9:
        # f_derivs debe devolver [y', y'', y''']
        y_p, y_pp, y_ppp = f_derivs(x, y)
        y = y + h * y_p + (h**2 / 2) * y_pp + (h**3 / 6) * y_ppp
        x = x + h
        x_values.append(x)
        y_values.append(y)
        
    return pd.DataFrame({'x': x_values, 'y': y_values})

def euler_system(f: Callable, t0: float, u0: np.ndarray, h: float, t_end: float) -> tuple[list, list]:
    """Resuelve un sistema de EDOs con el método de Euler."""
    t_values = [t0]
    u_values = [u0]
    
    t = t0
    u = u0.copy()
    n_steps = int(round((t_end - t0) / h))
    
    for _ in range(n_steps):
        u = u + h * f(t, u)
        t = t + h
        t_values.append(t)
        u_values.append(u)
        
    return t_values, u_values

def rk4_system(f: Callable, t0: float, u0: np.ndarray, h: float, t_end: float) -> tuple[list, list]:
    """Resuelve un sistema de EDOs con el método RK4."""
    t_values = [t0]
    u_values = [u0]
    
    t = t0
    u = u0.copy()
    n_steps = int(round((t_end - t0) / h))

    for _ in range(n_steps):
        k1 = f(t, u)
        k2 = f(t + h/2, u + (h/2) * k1)
        k3 = f(t + h/2, u + (h/2) * k2)
        k4 = f(t + h, u + h * k3)
        
        u = u + (h / 6) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + h
        t_values.append(t)
        u_values.append(u)
        
    return t_values, u_values

def taylor_orden2_system(f_derivs: Callable, t0: float, u0: np.ndarray, h: float, n_steps: int) -> tuple[list, list]:
    """Genera n_steps de un sistema de EDOs con Taylor de 2do grado (Útil para iniciar Milne)."""
    t_values = [t0]
    u_values = [u0]
    
    t = t0
    u = u0.copy()

    for _ in range(n_steps):
        u_p, u_pp = f_derivs(t, u)
        u = u + h * u_p + (h**2 / 2) * u_pp
        t = t + h
        t_values.append(t)
        u_values.append(u.copy())
        
    return t_values, u_values

def milne_system(f: Callable, t_start: list, u_start: list, h: float, t_end: float) -> tuple[list, list]:
    """Resuelve un sistema de EDOs con el método de Milne (Predictor-Corrector)."""
    if len(t_start) != 4 or len(u_start) != 4:
        raise ValueError("El método de Milne requiere 4 puntos de inicio.")

    t_values = list(t_start)
    u_values = list(u_start)
    f_values = [f(t, u) for t, u in zip(t_start, u_start)]
    
    t = t_values[-1]
    
    while t < t_end - 1e-9:
        # Predictor (Adams-Bashforth 4)
        u_pred = u_values[-1] + (h/24) * (55*f_values[-1] - 59*f_values[-2] + 37*f_values[-3] - 9*f_values[-4])
        t = t + h
        f_pred = f(t, u_pred)
        
        # Corrector (Adams-Moulton 3)
        u_corr = u_values[-1] + (h/24) * (9*f_pred + 19*f_values[-1] - 5*f_values[-2] + f_values[-3])
        
        # Actualización
        t_values.append(t)
        u_values.append(u_corr)
        f_values.append(f(t, u_corr))
        f_values = f_values[-4:] # Mantener buffer
        
    return t_values, u_values

# =============================================================================
# SECCIÓN 2: DEFINICIÓN DE PROBLEMAS DE EJEMPLO
# =============================================================================

def ejemplo_1_edo_simple():
    """Caso: 3y' - 4xy + e^x = 0"""
    print("\n=== CASO 1: EDO de 1er Orden ===")
    
    def f1(x, y): return (4 * x * y - np.exp(x)) / 3
    
    def f1_derivs_taylor(x, y):
        y_p = (4 * x * y - np.exp(x)) / 3
        y_pp = (4 * y + 4 * x * y_p - np.exp(x)) / 3
        y_ppp = (8 * y_p + 4 * x * y_pp - np.exp(x)) / 3
        return [y_p, y_pp, y_ppp]

    params = {'x0': 0.0, 'y0': 0.1, 'h': 0.05, 'x_end': 0.5}
    
    print("-> Euler:")
    print(euler(f1, **params).head())
    print("-> Taylor 3er Grado:")
    print(taylor_orden3(f1_derivs_taylor, **params).head())

def ejemplo_2_sistema_acoplado():
    """Caso: Sistema de EDOs acopladas (Milne y Euler)"""
    print("\n=== CASO 2: Sistema de EDOs Acopladas ===")
    
    def f2_system(t, u):
        y, x = u[0], u[1]
        return np.array([2 * t * y**2 - t**2 * y, 2 * t - x + y])

    t0, u0, h, t_end = 0.0, np.array([0.1, -0.2]), 0.1, 1.0
    
    print("-> Euler Multivariable:")
    t, u = euler_system(f2_system, t0, u0, h, t_end)
    df = pd.DataFrame(u, columns=['y(t)', 'x(t)'])
    df.insert(0, 't', t)
    print(df.head())

def ejemplo_3_problema_frontera():
    """Caso: Problema de Valor en la Frontera (BVP) con Diferencias Finitas"""
    print("\n=== CASO 3: Problema de Valor en la Frontera (BVP) ===")
    # y'' - y' + y = 17.5e^(-3t)
    # y(0) = 2.5, y(1) = 0.1245
    
    h = 1/3
    t_vals = [0.0, 1/3, 2/3, 1.0]
    y_0, y_3 = 2.5, 0.1245  
    
    # Coeficientes diferencias finitas
    C1 = (1 / h**2) + (1 / (2 * h))
    C2 = (-2 / h**2) + 1
    C3 = (1 / h**2) - (1 / (2 * h))

    # Sistema Ax = B
    A = np.array([[C2, C3], [C1, C2]])
    B = np.array([
        17.5 * np.exp(-3 * t_vals[1]) - C1 * y_0,
        17.5 * np.exp(-3 * t_vals[2]) - C3 * y_3
    ])

    try:
        y_int = np.linalg.solve(A, B)
        print(f"Solución nodos internos: y1={y_int[0]:.4f}, y2={y_int[1]:.4f}")
    except np.linalg.LinAlgError:
        print("Sistema singular")

# =============================================================================
# EJECUCIÓN PRINCIPAL
# =============================================================================

if __name__ == "__main__":
    print("--- INICIANDO SUITE DE PRUEBAS NUMÉRICAS ---")
    ejemplo_1_edo_simple()
    ejemplo_2_sistema_acoplado()
    ejemplo_3_problema_frontera()
    print("\n--- FIN DE LA EJECUCIÓN ---")