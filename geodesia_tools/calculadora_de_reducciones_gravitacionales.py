"""
Autor: Ivan Rugerio
Título: Calculadora de Reducciones Gravitacionales y Anomalías
Descripción: Herramienta de línea de comandos (CLI) para calcular correcciones gravimétricas 
             (Deriva, Aire Libre, Bouguer, Eotvos) y anomalías de gravedad utilizando 
             diversos modelos de referencia (GRS80, WGS84, Helmert, etc.).
Tecnologías: Python, Math
"""

import math

# --- CALCULADORA DE REDUCCIONES DE GRAVEDAD ---

print ("**** Bienvenid@ a la Calculadora de Reducciones de Gravedad ****")

# --- MENÚ PRINCIPAL ---

print("\n¿Qué deseas calcular hoy?")
print("[1] Corrección por Deriva Instrumental")
print("[2] Anomalías de Gravedad Terrestre")
print("[3] Anomalías de Gravedad Marina/Aérea")

modo = input("Selecciona una opción (1-3): ")

datos = []
g_base_i, g_base_f, t_i, t_f = 0, 0, 0, 0

# --- CAPTURA DE DATOS ---

if modo == '1':
    print("\n--- MODO: DERIVA INSTRUMENTAL ---")
    g_base_i = float(input("Lectura Base Inicial (mGal): "))
    t_i = float(input("Hora Inicial (decimal): "))
    g_base_f = float(input("Lectura Base Final (mGal): "))
    t_f = float(input("Hora Final (decimal): "))

    print("\nIntroduce tus estaciones. 'fin' para terminar.")
    while True:
        est = input("Nombre Estación (o 'fin'): ")
        if est.lower() == 'fin': break
        g_leida = float(input("Gravedad Leída: "))
        t_obs = float(input("Hora de lectura: "))
        datos.append({'id': est, 'g_raw': g_leida, 't': t_obs})

elif modo == '2':
    print("\n--- MODO: GRAVEDAD TERRESTRE ---")
    print("Escribe 'fin' en latitud para terminar.")
    while True:
        lat_str = input("\nLatitud (grados decimales) (o 'fin'): ")
        if lat_str.lower() == 'fin': break
        lon = float(input("Longitud (grados decimales): "))
        H = float(input("Altura Ortométrica (m): "))
        g_obs = float(input("Gravedad Observada (mGal): "))
        C_T = float(input("Corrección de Terreno (0 si no tienes): "))
        datos.append({'lat': float(lat_str), 'lon': lon, 'H': H, 'g': g_obs, 'CT': C_T})

elif modo == '3':
    print("\n--- MODO: GRAVEDAD DINÁMICA ---")
    print("Escribe 'fin' en latitud para terminar.")
    while True:
        lat_str = input("\nLatitud (grados decimales) (o 'fin'): ")
        if lat_str.lower() == 'fin': break
        lon = float(input("Longitud (grados decimales): "))
        g_obs = float(input("Gravedad Observada (mGal): "))
        vel = float(input("Velocidad (Nudos): "))
        rumbo = float(input("Rumbo (Azimut): "))
        datos.append({'lat': float(lat_str), 'lon': lon, 'g': g_obs, 'v': vel, 'az': rumbo})

# --- SELECCIÓN DE MÉTODO DE CALCULO DE GRAVEDAD TEORICA ---

metodo_g = '5' # Default GRS80 General
if modo == '2' or modo == '3':
    print("\nElige la fórmula de Gravedad Normal:")
    print("  [1] Helmert (1901)")
    print("  [2] Bowie (1917)")
    print("  [3] Cassini (Internacional 1930)")
    print("  [4] GRS67 (1967)")
    print("  [5] GRS80 (Fórmula General)")
    print("  [6] GRS80 (Precisión 0.1 µgal)")
    print("  [7] Somigliana (WGS84)")

    metodo_g = input("Opción (1-7): ")
    if metodo_g not in ['1','2','3','4','5','6','7']: metodo_g = '5'

# --- FUNCIONES DE CÁLCULO CIENTÍFICO ---

# 1. GRAVEDAD NORMAL
def calcular_gravedad_normal(latitud_grados, metodo):
    phi = math.radians(latitud_grados)
    sin_phi_2 = math.sin(phi)**2
    sin_2phi_2 = math.sin(2 * phi)**2

    if metodo == '1': # Helmert
        return 978030 * (1 + 0.005302 * sin_phi_2 - 0.000007 * sin_2phi_2)
    elif metodo == '2': # Bowie
        return 978039 * (1 + 0.005294 * sin_phi_2 - 0.000007 * sin_2phi_2)
    elif metodo == '3': # Internacional 1930
        return 978049 * (1 + 0.0052884 * sin_phi_2 - 0.0000059 * sin_2phi_2)
    elif metodo == '4': # GRS67
        return 978031.8 * (1 + 0.0053024 * sin_phi_2 - 0.0000058 * sin_2phi_2)
    elif metodo == '5': # GRS80 (General)
        return 978032.7 * (1 + 0.0053024 * sin_phi_2 - 0.0000058 * sin_2phi_2)
    elif metodo == '6': # GRS80 (Precisa)
        ye = 978032.67715
        b = 0.0052790414
        b1 = 0.0000232718
        b2 = 0.0000001262
        b3 = 0.0000000007
        sin_phi_4 = sin_phi_2**2
        sin_phi_6 = sin_phi_2**3
        sin_phi_8 = sin_phi_2**4
        return ye * (1 + b*sin_phi_2 + b1*sin_phi_4 + b2*sin_phi_6 + b3*sin_phi_8)
    elif metodo == '7': # Somigliana WGS84
        ge = 978032.53359
        k = 0.00193185138639
        e2 = 0.00669437999013
        return ge * ((1 + k * sin_phi_2) / math.sqrt(1 - e2 * sin_phi_2))
    else:
        return 978032.7

# 2. CORRECCIONES GEODÉSICAS
def calcular_deriv_instrumental(g_base_ini, g_base_fin, t_ini, t_fin, t_obs):
    if t_fin == t_ini: return 0.0
    tasa_deriva = (g_base_fin - g_base_ini) / (t_fin - t_ini)
    delta_t = t_obs - t_ini
    return - (tasa_deriva * delta_t)

def calcular_cor_aire_libre(H_metros):
    return 0.3086 * H_metros

def calcular_cor_bouguer(H_metros, rho_g_cm3=2.67):
    return 0.0419 * rho_g_cm3 * H_metros

def calcular_cor_curvatura(H_metros):
    return (1.464 * 10**-3 * H_metros) - (3.533 * 10**-7 * H_metros**2)

def calcular_cor_atmosferica(H_metros):
    return 0.87 - (0.0000965 * H_metros)

def calcular_cor_eotvos(latitud, velocidad_nudos, rumbo_azimut):
    phi = math.radians(latitud)
    alpha = math.radians(rumbo_azimut)
    return 7.503 * velocidad_nudos * math.sin(alpha) * math.cos(phi) + 0.004154 * velocidad_nudos**2

# --- PROCESAMIENTO Y GENERACIÓN DE REPORTE ---

print("\n" + "="*30)
print("RESULTADOS GENERADOS (Formato CSV):")
print("="*30 + "\n")

if modo == '1':
    print("Estacion,Hora,G_Cruda,Corr_Deriva,G_Corregida")
    for d in datos:
        C_D = calcular_deriv_instrumental(g_base_i, g_base_f, t_i, t_f, d['t'])
        g_corr = d['g_raw'] + C_D
        print(f"{d['id']},{d['t']},{d['g_raw']},{C_D:.4f},{g_corr:.4f}")

elif modo == '2':
    print("Latitud,Longitud,Altura_m,G_Obs,G_Normal,C_AireLibre,C_Bouguer,C_Terreno,C_Curvatura,C_Atmosferica,Anom_AireLibre,Anom_Bouguer_Simple,Anom_Bouguer_Completa")

    for d in datos:
        Gn = calcular_gravedad_normal(d['lat'], metodo_g)
        F = calcular_cor_aire_libre(d['H'])
        CB = calcular_cor_bouguer(d['H'])
        CCurv = calcular_cor_curvatura(d['H'])
        CAtm = calcular_cor_atmosferica(d['H'])

        # Cálculos de Anomalías
        Ag_AL = d['g'] - Gn + F
        Ag_BS = Ag_AL - CB
        Ag_BC = Ag_BS + d['CT'] + CCurv + CAtm

        print(f"{d['lat']},{d['lon']},{d['H']},{d['g']},{Gn:.4f},{F:.4f},{CB:.4f},{d['CT']},{CCurv:.4f},{CAtm:.4f},{Ag_AL:.4f},{Ag_BS:.4f},{Ag_BC:.4f}")

elif modo == '3':
    print("Latitud,Longitud,Velocidad_kt,Rumbo,G_Obs,G_Normal,C_Eotvos,Anomalia_AireLibre_Marina")

    for d in datos:
        Gn = calcular_gravedad_normal(d['lat'], metodo_g)
        Ce = calcular_cor_eotvos(d['lat'], d['v'], d['az'])
        Anomalia = d['g'] + Ce - Gn
        print(f"{d['lat']},{d['lon']},{d['v']},{d['az']},{d['g']},{Gn:.4f},{Ce:.4f},{Anomalia:.4f}")

print("\n--- FIN DEL PROCESO ---")