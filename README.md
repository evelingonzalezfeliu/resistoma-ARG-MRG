# Resistoma ARG–MRG

Repositorio de los trabajos SLURM y scripts usados para detectar, depurar e integrar genes de resistencia a antimicrobianos (ARG), genes asociados a resistencia a metales (MRG) y elementos genéticos móviles en una colección de aislados bacterianos.

Los genomas, bases de datos, contenedores y resultados completos no se almacenan en Git. El repositorio conserva los comandos reproducibles, los formatos de entrada y los scripts de análisis.

## Estructura

```text
resistoma-ARG-MRG/
├── slurm/                  # Trabajos ejecutados en el clúster
├── scripts/
│   ├── preprocessing/     # Conversión de GBFF a FNA, FAA y GFF3
│   └── results/           # QC, integración, colocalización y figuras
├── inputs/                 # Formatos y lista de muestras de ejemplo
├── README.md
└── LICENSE
```

## Entradas esperadas

La estructura de trabajo usada en el clúster fue:

```text
PROJECT_ROOT/
├── lista.txt
├── gbff/
│   └── SAMPLE.gbff
├── secuencias_extraidas/
│   ├── fna/SAMPLE.fna
│   ├── faa/SAMPLE.faa
│   └── gff/SAMPLE.gff
└── analysis/
```

Los formatos requeridos por cada paso se detallan en [inputs/README.md](inputs/README.md).

## Paso a paso de los trabajos SLURM

Esta es la sección principal para reproducir las búsquedas ejecutadas en el clúster. Las rutas `BASE` de los archivos SLURM deben apuntar al directorio correcto del proyecto.

Antes de comenzar:

```bash
export PROJECT_ROOT=/mnt/beegfs/home/efeliu/work2026/Unidad_MLatorre
export REPO=/ruta/al/clone/resistoma-ARG-MRG
cd "$PROJECT_ROOT/analysis"
mkdir -p logs containers
```

Los comandos siguientes se ejecutan desde `$PROJECT_ROOT/analysis`. `lista.txt` debe contener un identificador de aislado por línea, sin encabezado, líneas vacías ni duplicados. Los arrays se construyeron para 100 aislados.

### Resumen de los SLURM

| Orden | Script | Función | Entrada principal | Salida principal | Tipo |
|---:|---|---|---|---|---|
| 1 | `01_amrfinder_singularity.sbash` | Descarga el contenedor AMRFinderPlus 4.2.7 | Imagen Docker de NCBI | `analysis/containers/amrfinder_4.2.7_db2026-08-07.1.sif` | Preparación |
| 2 | `02_descargar_amrfinder_db.slurm` | Descarga y verifica la base AMRFinderPlus | Contenedor AMRFinderPlus | `analysis/amrfinder_db/` | Preparación |
| 3 | `descargar_bacmet.slurm` | Descarga BacMet2 EXP, mapping y checksums | Servidor de BacMet | `analysis/bacmet_db/` | Preparación |
| 4 | `descargar_mobsuite_db.slurm` | Inicializa la base local de MOB-suite | Contenedor MOB-suite | `analysis/mobsuite_db/` | Preparación |
| 5 | `03_ejecutar_amrfinder.slurm` | Ejecuta AMRFinderPlus usando solamente ADN | FNA | `analysis/resultados_amrfinder/` | Opcional/alternativo |
| 6 | `04_amrfinder_combinado_array.slurm` | Ejecuta AMRFinderPlus combinado con ADN, proteínas, GFF y `--plus` | FNA, FAA, GFF y `lista.txt` | `analysis/resultados_amrfinder_combinado_plus/` | Principal ARG |
| 7 | `05_diamond_bacmet_array.slurm` | Compara proteínas contra BacMet2 EXP | FAA, índice BacMet y `lista.txt` | `analysis/resultados_diamond_bacmet/` | Principal MRG |
| 8 | `06_isescan_array.slurm` | Detecta secuencias de inserción | FNA y `lista.txt` | `analysis/resultados_isescan/` | Contexto móvil |
| 9 | `07_mobsuite_array.slurm` | Clasifica contigs cromosómicos y plasmídicos | FNA, base MOB-suite y `lista.txt` | `analysis/resultados_mobsuite/` | Contexto móvil |
| 10 | `08_mobileelementfinder_array.slurm` | Detecta elementos genéticos móviles | FNA y `lista.txt` | `analysis/resultados_mobileelementfinder/` | Contexto móvil |
| 11 | `09_integronfinder_array.slurm` | Detecta integrasas, sitios `attC` e integrones | FNA y `lista.txt` | `analysis/resultados_integronfinder/` | Integrones |

### Paso 1 — Obtener el contenedor de AMRFinderPlus

Este archivo es un script Bash, no un trabajo SLURM completo. Se ejecuta directamente en el nodo autorizado para descargar imágenes:

```bash
bash "$REPO/slurm/01_amrfinder_singularity.sbash"
```

Comprobar:

```bash
test -s containers/amrfinder_4.2.7_db2026-08-07.1.sif
```

### Paso 2 — Descargar la base de AMRFinderPlus

```bash
sbatch "$REPO/slurm/02_descargar_amrfinder_db.slurm"
```

El trabajo ejecuta `amrfinder_update` y registra las versiones del programa y de la base. Debe finalizar antes de enviar el análisis combinado.

Comprobar que exista:

```text
amrfinder_db/2026-08-07.1/
```

### Paso 3 — Descargar y preparar BacMet2 EXP

Descargar las 753 referencias experimentales y su mapping:

```bash
sbatch "$REPO/slurm/descargar_bacmet.slurm"
```

El trabajo crea `checksums.sha256` y `fecha_descarga.txt`. Después de la descarga se debe construir el índice DIAMOND:

```bash
singularity exec --cleanenv \
  containers/diamond_2.1.9.sif \
  diamond makedb \
  --in bacmet_db/BacMet2_EXP_database.fasta \
  --db bacmet_db/BacMet2_EXP
```

Comprobar:

```bash
test -s bacmet_db/BacMet2_EXP.dmnd
```

La construcción del índice fue parte del análisis original, aunque su SLURM todavía no está incluido como archivo independiente.

### Paso 4 — Inicializar la base de MOB-suite

```bash
sbatch "$REPO/slurm/descargar_mobsuite_db.slurm"
```

Este trabajo ejecuta `mob_init`. Debe finalizar antes de enviar `07_mobsuite_array.slurm`.

### Paso 5 — AMRFinderPlus combinado: análisis principal de ARG

```bash
sbatch --array=1-100%25 \
  "$REPO/slurm/04_amrfinder_combinado_array.slurm" \
  "$PROJECT_ROOT/lista.txt"
```

Para cada aislado se utilizan conjuntamente:

```text
secuencias_extraidas/fna/SAMPLE.fna
secuencias_extraidas/faa/SAMPLE.faa
secuencias_extraidas/gff/SAMPLE.gff
```

La salida de cada muestra incluye:

```text
SAMPLE.amrfinder.tsv
SAMPLE.amrfinder.log
```

El resultado se publica desde el archivo temporal solamente cuando AMRFinderPlus termina correctamente. Este es el análisis utilizado para los ARG principales.

El script `03_ejecutar_amrfinder.slurm` analiza únicamente los FNA y se conserva como alternativa histórica:

```bash
sbatch "$REPO/slurm/03_ejecutar_amrfinder.slurm"
```

Sus resultados no deben sumarse a los de la corrida combinada porque ambas búsquedas se solapan.

### Paso 6 — DIAMOND contra BacMet2 EXP

Requiere que el índice `BacMet2_EXP.dmnd` exista:

```bash
sbatch "$REPO/slurm/05_diamond_bacmet_array.slurm" \
  "$PROJECT_ROOT/lista.txt"
```

El array está definido como `1-100%25`. Para cada aislado produce:

```text
SAMPLE.bacmet.raw.tsv
SAMPLE.diamond.log
SAMPLE.version.txt
```

La búsqueda usa `--very-sensitive`, `--evalue 1e-5` y `--max-target-seqs 0`. Los filtros de identidad y cobertura se aplican posteriormente en `scripts/results/analisis_arg_mrg.py`.

### Paso 7 — ISEScan

```bash
sbatch "$REPO/slurm/06_isescan_array.slurm" \
  "$PROJECT_ROOT/lista.txt"
```

El array está definido como `1-100%25`. Cada aislado se procesa en un directorio independiente:

```text
analysis/resultados_isescan/SAMPLE/
```

La tabla utilizada posteriormente para las coordenadas de las secuencias de inserción es:

```text
analysis/resultados_isescan/SAMPLE/fna/SAMPLE.fna.tsv
```

### Paso 8 — MOB-suite

```bash
sbatch "$REPO/slurm/07_mobsuite_array.slurm" \
  "$PROJECT_ROOT/lista.txt"
```

El array está definido como `1-100%25`. Las salidas principales son:

```text
analysis/resultados_mobsuite/SAMPLE/contig_report.txt
analysis/resultados_mobsuite/SAMPLE/mobtyper_results.txt
analysis/resultados_mobsuite/SAMPLE/mge.report.txt
```

`contig_report.txt` se utilizó para distinguir contexto cromosómico y plasmídico. El script se detiene si ya existe el directorio de salida de la muestra, evitando mezclar ejecuciones.

### Paso 9 — MobileElementFinder

```bash
sbatch --array=1-100%25 \
  "$REPO/slurm/08_mobileelementfinder_array.slurm" \
  "$PROJECT_ROOT/lista.txt"
```

Para cada aislado se esperan:

```text
analysis/resultados_mobileelementfinder/SAMPLE/SAMPLE.csv
analysis/resultados_mobileelementfinder/SAMPLE/SAMPLE.gff
analysis/resultados_mobileelementfinder/SAMPLE/SAMPLE_result.txt
analysis/resultados_mobileelementfinder/SAMPLE/SAMPLE_mge_sequences.fna
```

En la ejecución original se crearon 100 carpetas, pero solo 87 contienen las salidas finales. Los 13 trabajos restantes deben revisarse en `logs/mefinder_JOB_ARRAY.err` antes de interpretar ausencia de elementos móviles.

### Paso 10 — IntegronFinder

```bash
sbatch --array=1-100%25 \
  "$REPO/slurm/09_integronfinder_array.slurm" \
  "$PROJECT_ROOT/lista.txt"
```

Para cada aislado se esperan:

```text
analysis/resultados_integronfinder/SAMPLE/
└── Results_Integron_Finder_SAMPLE/
    ├── SAMPLE.integrons
    ├── SAMPLE.summary
    └── integron_finder.out
```

`SAMPLE.summary` resume integrones completos, `In0` y regiones `CALIN`. En la ejecución original hay resultados para 99 de los 100 aislados; falta `B28_A03`, que debe repetirse o documentarse como análisis no disponible.

### Paso 11 — Control de finalización

Antes de procesar los resultados se debe verificar que cada herramienta tenga una salida por aislado y revisar los archivos `.err` y `.log`. Una carpeta vacía o incompleta no se debe interpretar como ausencia biológica.

Controles mínimos:

```bash
find resultados_isescan -mindepth 1 -maxdepth 1 -type d | wc -l
find resultados_mobsuite -mindepth 1 -maxdepth 1 -type d | wc -l
find resultados_integronfinder -name '*.summary' | wc -l
find resultados_mobileelementfinder -name '*_result.txt' | wc -l
```

Los conteos esperados son 100 por herramienta, salvo que una muestra se documente explícitamente como fallida o no analizada.

## Scripts auxiliares

### Preparar FNA, FAA y GFF3

Desde archivos GBFF anotados con Bakta:

```bash
python scripts/preprocessing/extraer_fna_faa.py gbff secuencias_extraidas
python scripts/preprocessing/generar_gff_desde_gbff.py \
  gbff secuencias_extraidas/gff
```

Dependencia: Biopython.

### Procesar los resultados

| Script | Entrada principal | Salida |
|---|---|---|
| `scripts/results/analisis_arg_mrg.py` | AMRFinderPlus, DIAMOND/BacMet, FAA, FNA, GFF y metadatos | QC, tablas integradas, matriz y frecuencias |
| `scripts/results/analizar_colocalizacion_arg_mrg.py` | ARG, MRG y GFF | Pares ARG–MRG y distancias |
| `scripts/results/integrar_contexto_movil_arg_mrg.py` | Pares ARG–MRG, MOB-suite e ISEScan | Contexto móvil por aislado |
| `scripts/results/figuras_arg_mrg.R` | Matriz y frecuencias | Heatmaps y gráfico de frecuencias |
| `scripts/results/plot_atlas_gggenomes_ARG_MRG.R` | Pares, GFF, MOB-suite e ISEScan | Atlas de contexto genómico |

Algunos scripts conservan rutas del análisis original. Deben ajustarse o parametrizarse antes de usarlos en otro sistema.

## Programas principales

| Programa | Versión documentada |
|---|---:|
| AMRFinderPlus | 4.2.7 |
| Base AMRFinderPlus | 2026-08-07.1 |
| DIAMOND | 2.1.9 |
| ISEScan | 1.7.3 |
| MOB-suite | 3.1.9 |
| MobileElementFinder | 1.0.3; base 1.0.2 |
| IntegronFinder | 2.0.6 |

## Criterios principales de ARG y MRG

- ARG principales: filas `Type=AMR` de AMRFinderPlus combinado con `--plus`.
- Los hallazgos `STRESS/METAL` de AMRFinderPlus se conservaron como evidencia complementaria.
- MRG principales: candidatos BacMet con identidad ponderada ≥80 %, cobertura de query ≥80 %, cobertura de referencia ≥80 % y E-value ≤1e-5.
- Cada aislado se contó una vez por característica, sin multiplicar copias ni asociaciones multimetal.
- Las predicciones genómicas no demuestran resistencia fenotípica, movilidad ni transferencia horizontal.
