# Método de Newton para la Inversa Matricial 

Este repositorio contiene la implementación algorítmica y el estudio teórico del **Método de Newton aplicado a la inversión de matrices**, desarrollado como proyecto de la asignatura Resolución Numérica (Doble Grado en Matemáticas e Ingeniería Informática, Universitat Politècnica de València).

**Calificación obtenida:** 10/10 

##  Descripción del Proyecto
El trabajo aborda la resolución de la ecuación matricial no lineal $F(X) = X^{-1} - A = 0$ mediante métodos iterativos. Se ha programado en MATLAB y estructurado en varias fases de análisis matemático y computacional:

* **Desarrollo y Convergencia:** Deducción de la expresión iterativa y demostración analítica de su convergencia cuadrática.
* **Análisis Dinámico:** Estudio de la estabilidad del método (puntos fijos superatractores y atractores) y representación gráfica de cuencas de atracción (planos dinámicos).
* **Inversa de Moore-Penrose:** Extensión del algoritmo para el cálculo de la pseudoinversa en matrices no cuadradas, garantizando el cumplimiento de las condiciones de Moore-Penrose.
* **Modificaciones Multipaso:** Análisis de alternativas computacionales como el método de Traub y métodos congelando la evaluación.

##  Estructura del Repositorio
* `/codigo`: Contiene los scripts desarrollados en MATLAB.
  * `Newton_matrcial.m`: Implementación del método estándar para matrices cuadradas.
  * `Newton_Pseudoinversa.m`: Variación para matrices no cuadradas (Moore-Penrose).
  * `trabajo_resonum.xml`: Archivo matlab con pruebas hechas sobre el código anterior. 
* `/docs`: Contiene la memoria completa del proyecto en formato PDF (`Trabajo_Resonum.pdf`) con todas las deducciones matemáticas y tablas de resultados.

##  Uso rápido de los scripts (MATLAB)

**¡Importante! Condición de convergencia:** 
Para que el método converja a la solución y no diverja, la aproximación inicial $X_0$ debe cumplir obligatoriamente la condición necesaria $||I - A X_0|| < 1$. 

Una forma efectiva de garantizar esto (tal y como se demuestra en el estudio) es definir la estimación inicial como $X_0 = \beta A^t$, utilizando un escalar pequeño como $\beta = \frac{1}{||A||^2}$.

Para probar el algoritmo principal, puedes ejecutar el siguiente fragmento en MATLAB:

```matlab
% 1. Definimos la matriz A que queremos invertir
A = [4 1; 1 3]; 

% 2. Calculamos una aproximación inicial X0 válida para garantizar convergencia
beta = 1 / (norm(A)^2);
X0 = beta * A'; 

% 3. Ejecutamos el método de Newton (máx 100 iteraciones, tolerancia 1e-10)
[iter, SOL, incre, I, tcl, tcc, increl, ACOC, dif_iter] = Newton_matrcial(A, X0, 100, 1e-10);

% Mostramos la solución
disp('Inversa calculada:');
disp(SOL);
