# SIDLE: Caracterización de la microbiota en muestras de tejido parafinado (FFPE)

Este repositorio contiene el código para el análisis bioinformático de muestras de tejido tumoral de cáncer de mama conservadas mediante fijación en formalina e inclusión en parafina (FFPE) obtenidas en el instituto de genética médica y molecular (INGEMM) en el Hospital de la Paz, Madrid, España.

Las muestras FFPE presentan un ADN altamente fragmentado, generalmente con longitudes de entre 100 y 300 pb. Esta característica dificulta el empleo de métodos convencionales de secuenciación basados en regiones largas (*single long-region framework*).

Para abordar esta limitación, se empleó **SMURF** (*Short MUltiple Regions Framework*), un algoritmo diseñado para combinar la información obtenida mediante la amplificación y secuenciación de múltiples regiones del gen 16S rRNA, permitiendo reconstruir perfiles taxonómicos de la microbiota a partir de fragmentos cortos de ADN.

La implementación de este enfoque se realizó mediante **SIDLE** (*SMURF Implementation Done to acceLerate Efficiency*), un plugin integrado en el entorno bioinformático **QIIME 2**.

Documentación oficial: https://q2-sidle.readthedocs.io/en/latest/

## Entorno de trabajo y herramientas

El análisis se desarrolló utilizando las siguientes herramientas:

- **Docker:** creación de un entorno de ejecución aislado y reproducible, con versiones compatibles de las herramientas y dependencias necesarias.
- **QIIME 2:** plataforma bioinformática empleada para el procesamiento y análisis de los datos de secuenciación de microbiota.
- **SIDLE:** plugin de QIIME 2 utilizado para implementar el enfoque SMURF.
- **Visual Studio Code:** entorno de desarrollo utilizado para la edición y gestión del código.
