# Entradas por etapa

Los archivos biológicos completos no se incluyen en el repositorio. Esta carpeta documenta los nombres y formatos esperados.

| Etapa | Entrada mínima | Ejemplo de nombre |
|---|---|---|
| Preparación de secuencias | Archivo GenBank anotado | `gbff/B28_A03.gbff` |
| AMRFinderPlus combinado | FNA, FAA y GFF3 de la misma muestra | `B28_A03.fna`, `B28_A03.faa`, `B28_A03.gff` |
| DIAMOND/BacMet | FAA y base `BacMet2_EXP.dmnd` | `B28_A03.faa` |
| ISEScan | FNA | `B28_A03.fna` |
| MOB-suite | FNA y base local de MOB-suite | `B28_A03.fna` |
| MobileElementFinder | FNA | `B28_A03.fna` |
| IntegronFinder | FNA | `B28_A03.fna` |
| Integración ARG–MRG | Resultados AMRFinderPlus y DIAMOND, mapping BacMet, FNA, FAA, GFF y metadatos | un conjunto por aislado |
| Figuras | Matriz de presencia/ausencia y tabla de frecuencias | archivos TSV generados por el análisis |

## Convención de identificadores

El mismo identificador debe utilizarse en `lista.txt` y en todos los nombres de archivos:

```text
B28_A03
```

Los identificadores de proteínas y contigs deben incluir la muestra para evitar colisiones:

```text
>B28_A03|EGGMEA_03534
>B28_A03|contig_1
```

## Lista de muestras

`lista.txt` contiene un identificador por línea, sin encabezado. Se incluye `lista_ejemplo.txt` solamente para mostrar el formato.

Antes de enviar un array se debe verificar que la lista no contenga líneas vacías ni identificadores duplicados.
