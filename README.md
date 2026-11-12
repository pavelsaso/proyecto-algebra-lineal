# Cinemática inversa de un brazo robótico planar 2R (Álgebra Lineal)

Proyecto final de Álgebra Lineal, UDLAP, Otoño 2026.
Docente: Andreiv Domínguez Martínez.
Autores: Pavel Bautista Fuentes

## Descripción
Simulación en C++ (Eigen) de un brazo planar de dos eslabones. Se resuelve la
cinemática inversa con mínimos cuadrados amortiguados (SVD) y se comparan
Gauss-Newton con amortiguamiento fijo, Levenberg-Marquardt adaptativo y la
solución analítica. Incluye validaciones y análisis de sensibilidad.

## Requisitos
- Compilador g++ con soporte C++17 (MinGW en Windows)
- Eigen 3.4 (https://eigen.tuxfamily.org), solo encabezados
- Python 3 con numpy, pandas, matplotlib y pillow:
  `pip install numpy pandas matplotlib pillow`

## Cómo reproducir
1. Descarga Eigen y descomprímelo (por ejemplo, en `C:\C++\eigen`).
2. Compila y ejecuta:
   ```
   g++ -std=c++17 -O2 -I C:\C++\eigen main.cpp -o brazo.exe
   .\brazo.exe > resultados.txt
   ```
3. Genera las figuras:
   ```
   python grafica.py
   python animar.py
   python figs.py
   ```

## Archivos
| Archivo | Contenido |
|---|---|
| `main.cpp` | Código completo: modelo, métodos, validaciones y exportación de datos |
| `grafica.py` | Gráfica de κ(J) vs q2 (`kappa.png`) |
| `animar.py` | Figura de la trayectoria 1 y animación (`brazo.gif`) |
| `figs.py` | Esquema, mapa de convergencia, mapa de calor, trayectoria 2 y fotogramas |
| `kappa.csv`, `trayectoria.csv`, `trayectoria2.csv`, `mapa.csv` | Datos simulados generados por el programa |
| `resultados.txt` | Salida del programa |

## Resultados principales
- Levenberg-Marquardt: 480/480 objetivos convergen; GN con λ fijo: 425/480.
- Error de posición máximo menor que 1e-6; tiempo medio de alrededor de 1 µs.
- Jacobiano validado con diferencias finitas (error 3.4e-10).

## Datos y licencias
Todos los datos son simulados; no hay datos personales. Eigen se distribuye
bajo licencia MPL 2.0. Semilla del generador aleatorio: 42.