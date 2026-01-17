# 🌍 Portafolio de Ingeniería Geomática

![Python](https://img.shields.io/badge/Python-3.8%2B-blue?style=for-the-badge&logo=python&logoColor=white)
![PostgreSQL](https://img.shields.io/badge/PostgreSQL-PostGIS-336791?style=for-the-badge&logo=postgresql&logoColor=white)
![Status](https://img.shields.io/badge/Estado-Activo-success?style=for-the-badge)

Colección de algoritmos y herramientas profesionales desarrolladas para resolver problemas complejos de **Geodesia Física**, **Posicionamiento GNSS**, **Análisis Espacial** y **Cálculo Numérico**.

## 📂 Estructura del Repositorio

### 🛰️ 1. Herramientas de Geodesia (`geodesia_tools/`)
Scripts especializados en el procesamiento de datos físicos y satelitales.
* **`satelite.py`**: Algoritmo de posicionamiento de punto único (SPP) utilizando pseudodistancias y mínimos cuadrados iterativos. Incluye visualización de convergencia.
* **`calculadora_de_reducciones_gravitacionales.py`**: CLI robusta para calcular anomalías de gravedad (Aire Libre, Bouguer) y correcciones geodésicas (Eotvos, Deriva) soportando modelos GRS80 y WGS84.
* **`QC_CORS.py`**: Módulo para el control de calidad y monitoreo de Estaciones de Referencia de Operación Continua (CORS).

### 📐 2. Métodos Numéricos Avanzados (`metodos_numericos/`)
Implementación desde cero de algoritmos matemáticos fundamentales para ingeniería.
* **`analisis_ecuaciones_diferenciales.py`**: Solvers para EDOs y sistemas acoplados (Euler, Runge-Kutta 4, Milne Predictor-Corrector).
* **`metodos_interpolacion_cuadratura.py`**: Biblioteca para interpolación (Lagrange, Newton) e integración numérica de alta precisión (Cuadratura Gaussiana, Simpson Mixto).

### 🗄️ 3. Base de Datos Geoespacial (`database/`)
* **`ir_pf.sql`**: Script SQL para **PostgreSQL + PostGIS**.
  * Gestión de solicitudes financieras y préstamos.
  * **Análisis de Riesgo Espacial:** Algoritmo que intercepta ubicaciones de clientes con polígonos de riesgo en la CDMX para determinar la viabilidad de créditos automáticamente.

---

## 🚀 Instalación y Uso

1. **Clonar el repositorio:**
   ```bash
   git clone [https://github.com/irugerio28/Geomatics-Engineering-Portfolio.git](https://github.com/irugerio28/Geomatics-Engineering-Portfolio.git)
   cd Geomatics-Engineering-Portfolio
