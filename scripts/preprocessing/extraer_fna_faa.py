#!/usr/bin/env python3
"""Extrae contigs (.fna) y proteínas anotadas (.faa) desde archivos GBFF."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


DEFAULT_INPUT = Path(
    "/Users/efeliu/Documents/Mlatorre/Genomas y Cepario Bioing UOH"
)


def clean_sample_name(path: Path) -> str:
    """Convierte 'Copia de B28_A03.gbff' en 'B28_A03'."""
    name = re.sub(r"^Copia de\s+", "", path.stem, flags=re.IGNORECASE)
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def first_qualifier(qualifiers: dict[str, list[str]], key: str, default: str) -> str:
    values = qualifiers.get(key)
    return values[0] if values else default


def unique_id(identifier: str, counts: dict[str, int]) -> str:
    """Evita encabezados duplicados dentro de un mismo FASTA."""
    counts[identifier] = counts.get(identifier, 0) + 1
    number = counts[identifier]
    return identifier if number == 1 else f"{identifier}_dup{number}"


def process_gbff(gbff: Path, fna_dir: Path, faa_dir: Path) -> dict[str, object]:
    sample = clean_sample_name(gbff)
    records = list(SeqIO.parse(gbff, "genbank"))

    # Se prefija la muestra para que los contigs sean únicos incluso al juntar
    # posteriormente los FASTA de distintos genomas.
    nucleotide_records: list[SeqRecord] = []
    for record in records:
        contig_id = f"{sample}|{record.id}"
        nucleotide_records.append(
            SeqRecord(
                record.seq,
                id=contig_id,
                name=contig_id,
                description=record.description,
            )
        )

    fna_path = fna_dir / f"{sample}.fna"
    contig_count = SeqIO.write(nucleotide_records, fna_path, "fasta")

    proteins: list[SeqRecord] = []
    cds_total = 0
    cds_without_translation = 0
    id_counts: dict[str, int] = {}

    for record in records:
        contig_id = f"{sample}|{record.id}"

        for feature in record.features:
            if feature.type != "CDS":
                continue

            cds_total += 1
            qualifiers = feature.qualifiers
            translations = qualifiers.get("translation")

            if not translations:
                cds_without_translation += 1
                continue

            locus_tag = first_qualifier(qualifiers, "locus_tag", "")
            protein_id = first_qualifier(qualifiers, "protein_id", "")
            base_id = locus_tag or protein_id or f"CDS_{cds_total:06d}"
            protein_name = unique_id(base_id, id_counts)
            fasta_id = f"{sample}|{protein_name}"

            product = first_qualifier(qualifiers, "product", "proteína no anotada")
            product = " ".join(product.split())
            protein_id_text = protein_id or "NA"

            proteins.append(
                SeqRecord(
                    Seq("".join(translations[0].split())),
                    id=fasta_id,
                    name=fasta_id,
                    description=(
                        f"protein_id={protein_id_text} "
                        f"contig={contig_id} product={product}"
                    ),
                )
            )

    faa_path = faa_dir / f"{sample}.faa"
    protein_count = SeqIO.write(proteins, faa_path, "fasta")

    return {
        "muestra": sample,
        "archivo_gbff": str(gbff),
        "contigs_extraidos": contig_count,
        "cds_totales": cds_total,
        "proteinas_extraidas": protein_count,
        "cds_sin_translation": cds_without_translation,
        "archivo_fna": str(fna_path),
        "archivo_faa": str(faa_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extrae archivos FNA y FAA desde todos los GBFF de una carpeta."
    )
    parser.add_argument(
        "input_dir",
        nargs="?",
        type=Path,
        default=DEFAULT_INPUT,
        help="Carpeta que contiene los GBFF.",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        type=Path,
        default=Path.cwd() / "secuencias_extraidas",
        help="Carpeta donde se guardarán fna/, faa/ y el resumen.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_dir: Path = args.input_dir.expanduser().resolve()
    output_dir: Path = args.output_dir.expanduser().resolve()

    if not input_dir.is_dir():
        raise SystemExit(f"No existe la carpeta de entrada: {input_dir}")

    gbff_files = sorted(
        path for path in input_dir.iterdir() if path.is_file() and path.suffix.lower() == ".gbff"
    )
    if not gbff_files:
        raise SystemExit(f"No se encontraron archivos .gbff en: {input_dir}")

    fna_dir = output_dir / "fna"
    faa_dir = output_dir / "faa"
    fna_dir.mkdir(parents=True, exist_ok=True)
    faa_dir.mkdir(parents=True, exist_ok=True)

    summary: list[dict[str, object]] = []
    for index, gbff in enumerate(gbff_files, start=1):
        print(f"[{index}/{len(gbff_files)}] {gbff.name}", flush=True)
        summary.append(process_gbff(gbff, fna_dir, faa_dir))

    summary_path = output_dir / "resumen_extraccion.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    print("\nExtracción terminada.")
    print(f"Genomas procesados: {len(summary)}")
    print(f"Contigs extraídos: {sum(int(row['contigs_extraidos']) for row in summary)}")
    print(f"Proteínas extraídas: {sum(int(row['proteinas_extraidas']) for row in summary)}")
    print(f"CDS sin /translation: {sum(int(row['cds_sin_translation']) for row in summary)}")
    print(f"Resultados: {output_dir}")


if __name__ == "__main__":
    main()
