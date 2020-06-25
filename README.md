# Anotacion <i>Antirrhinum Majus</i>
Código utilizado en el método para la anotación del genoma de <i>Antirrhinum majus</i> a partir de los de <i>Solanum lycopersicum</i> y <i>Arabidopsis Thaliana</i> desarrollado en mi Trabajo de Fin de Grado en la Universidad de Murcia durante el curso 2019/2020.

## Scripts
- <b>filtroPorEvalor.py</b>: Se encarga de leer la entrada de los BLAST y filtrar todos aquellos cuyo e-value supere un umbral.
- <b>filtroPorPID.py</b>: Se encarga de leer los datos producidos por el anterior para aplicar un filtro de manera que se eliminen los que no superen un determinado porcentaje de identidad. A continuación aplica otro filtro para quedarse solo con los resultados recíprocos. Finalmente escoge para cada gen la pareja con la que tenga mejor e-value.
- <b>obtenerAnotaciones.py</b>: A partir de los datos producidos por el anterior, obtiene las anotaciones asociadas a los genes ortólogos mediante consultas a Biomart.
- <b>seleccion.py</b>: Utilizando las anotaciones que hemos recuperado, selecciona cuáles de ellas escogemos para cada gen. Además, contiene la funcionalidad necesaria para medir el rendimiento del método aplicándoselo a los genes ya anotados.

## Dependencias
Los scripts utilizan alguna o varias de las siguientes librerías:
- <b>Pandas</b>: Utilizada para el manejo de grandes estructuras de datos.
- <b>Pybiomart</b>: Utilizada para la comunicación con Biomart.
- <b>Goatools</b>: Utilizada para el tratamiento de términos GO.

## Ficheros de partida
- <b>anotacionGO.tsv</b>: Las anotaciones que ya existen para los genes de <i>Antirrhinum majus</i>.
- <b>genes_no_anotados.tsv</b>: Los genes de <i>Antirrhinum majus</i> para los que no existe anotación.

## Resultados
En este directorio está el fichero final con todas las anotaciones obtenidas para <i>Antirrhinum majus</i> utilizando nuestro método.
