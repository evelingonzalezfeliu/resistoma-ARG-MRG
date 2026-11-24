#!/usr/bin/env python3
"""Genera archivos GFF3 desde uno o varios GBFF anotados con Bakta."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from urllib.parse import quote

from Bio import SeqIO


DEFAULT_INPUT = Path(
    "/Users/efeliu/Documents/Mlatorre/Genomas_Cepario_Bioing_UOH"
)
DEFAULT_OUTPUT = Path(
    "/Users/efeliu/Documents/Mlatorre/analisis/gff"
)


def escape(text: object) -> str:
    """Codifica caracteres reservados dentro de atributos GFF3."""
    return quote(str(text), safe="|:._-")


def clean_sample_name(path: Path) -> str:
    """Convierte 'Copia de B28_A03.gbff' en 'B28_A03'."""
    name = re.sub(r"^Copia de\s+", "", path.stem, flags=re.IGNORECASE)
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def first_qualifier(
    qualifiers: dict[str, list[str]], key: str
) -> str | None:
    values = qualifiers.get(key)
    return values[0] if values else None


def build_attributes(sample: str, qualifiers: dict[str, list[str]]) -> str:
    locus = first_qualifier(qualifiers, "locus_tag")
    if locus is None:
        raise ValueError("CDS sin locus_tag")

    protein = f"{sample}|{locus}"
    attributes: list[tuple[str, str]] = [
        ("ID", protein),
        ("Name", protein),
        ("locus_tag", locus),
    ]

    for key in ("gene", "protein_id", "product"):
        value = first_qualifier(qualifiers, key)
        if value:
            attributes.append((key, value))

    return ";".join(f"{key}={escape(value)}" for key, value in attributes)


def convert_gbff(gbff: Path, output_dir: Path) -> dict[str, object]:
    sample = clean_sample_name(gbff)
    output_path = output_dir / f"{sample}.gff"
    temporary_path = output_dir / f".{sample}.gff.tmp"

    records = list(SeqIO.parse(gbff, "genbank"))
    if not records:
        raise ValueError("el GBFF no contiene registros GenBank")

    total_cds = 0
    written_cds = 0
    without_translation = 0
    locus_seen: set[str] = set()

    try:
        with temporary_path.open("w", encoding="utf-8", newline="\n") as out:
            out.write("##gff-version 3\n")

            for record in records:
                contig = f"{sample}|{record.id}"
                out.write(f"##sequence-region {contig} 1 {len(record.seq)}\n")

                for feature in record.features:
                    if feature.type != "CDS":
                        continue

                    total_cds += 1
                    qualifiers = feature.qualifiers

                    # Se mantiene la misma población usada para generar los FAA.
                    if "translation" not in qualifiers:
                        without_translation += 1
                        continue

                    locus = first_qualifier(qualifiers, "locus_tag")
                    if locus is None:
                        raise ValueError(
                            f"{record.id}: CDS sin locus_tag en {feature.location}"
                        )
                    if locus in locus_seen:
                        raise ValueError(f"locus_tag duplicado: {locus}")
                    locus_seen.add(locus)

                    # No se colapsan localizaciones compuestas porque se perdería
                    # la estructura real del CDS. Se exige revisión explícita.
                    if len(feature.location.parts) != 1:
                        raise ValueError(
                            f"{locus}: CDS con varias partes; requiere revisión"
                        )

                    strand_value = feature.location.strand
                    if strand_value not in (1, -1):
                        raise ValueError(f"{locus}: orientación no definida")

                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if strand_value == 1 else "-"

                    codon_start_text = first_qualifier(qualifiers, "codon_start") or "1"
                    codon_start = int(codon_start_text)
                    if codon_start not in (1, 2, 3):
                        raise ValueError(
                            f"{locus}: codon_start inválido: {codon_start_text}"
                        )
                    phase = codon_start - 1

                    attributes = build_attributes(sample, qualifiers)
                    fields = [
                        contig,
                        "Bakta",
                        "CDS",
                        str(start),
                        str(end),
                        ".",
                        strand,
                        str(phase),
                        attributes,
                    ]
                    out.write("\t".join(fields) + "\n")
                    written_cds += 1

        temporary_path.replace(output_path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    return {
        "muestra": sample,
        "archivo_gbff": str(gbff),
        "archivo_gff": str(output_path),
        "contigs": len(records),
        "cds_totales": total_cds,
        "cds_con_translation": written_cds,
        "cds_sin_translation": without_translation,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Genera un GFF3 por muestra desde un GBFF o desde todos los "
            "GBFF de una carpeta. Solo escribe CDS con /translation."
        )
    )
    parser.add_argument(
        "entrada",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help="Archivo .gbff o carpeta que contiene los GBFF.",
    )
    parser.add_argument(
        "salida",
        nargs="?",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Carpeta donde se guardarán los GFF3 y el resumen.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.entrada.expanduser().resolve()
    output_dir = args.salida.expanduser().resolve()

    if input_path.is_file():
        gbff_files = [input_path]
    elif input_path.is_dir():
        gbff_files = sorted(input_path.glob("*.gbff"))
    else:
        raise SystemExit(f"No existe la entrada: {input_path}")

    if not gbff_files:
        raise SystemExit(f"No se encontraron archivos .gbff en: {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    summary: list[dict[str, object]] = []

    for index, gbff in enumerate(gbff_files, start=1):
        print(f"[{index}/{len(gbff_files)}] {gbff.name}", flush=True)
        try:
            summary.append(convert_gbff(gbff, output_dir))
        except Exception as exc:
            raise SystemExit(f"Error en {gbff.name}: {exc}") from exc

    summary_path = output_dir / "resumen_gff.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    print("\nConversión terminada.")
    print(f"GBFF procesados: {len(summary)}")
    print(f"GFF generados: {len(summary)}")
    print(f"CDS escritos: {sum(int(row['cds_con_translation']) for row in summary)}")
    print(
        "CDS omitidos por no tener /translation: "
        f"{sum(int(row['cds_sin_translation']) for row in summary)}"
    )
    print(f"Resumen: {summary_path}")


if __name__ == "__main__":
    main()
