#!/usr/bin/env python3
"""Cruza ARG de AMRFinderPlus y MRG de BacMet con GFF de Bakta."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


BASE = Path("/Users/efeliu/Documents/Mlatorre/analisis")
AMR_PATH = BASE / "analisis_arg_mrg_20261022/tablas/amrfinderplus_depurado.tsv"
BACMET_PATH = BASE / "analisis_arg_mrg_20261022/tablas/bacmet_asignaciones_por_proteina.tsv"
GFF_DIR = BASE / "gff"
OUTPUT_DIR = BASE / "analisis_arg_mrg_20261022/colocalizacion"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if separator:
            attributes[key] = unquote(value)
    return attributes


def load_gff_features() -> tuple[dict[str, dict[str, object]], set[str]]:
    features: dict[str, dict[str, object]] = {}
    samples: set[str] = set()

    for path in sorted(GFF_DIR.glob("*.gff")):
        sample = path.stem
        samples.add(sample)
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"{path.name}:{line_number}: no tiene 9 columnas")
                attributes = parse_attributes(fields[8])
                feature_id = attributes.get("ID")
                if not feature_id:
                    raise ValueError(f"{path.name}:{line_number}: CDS sin ID")
                if feature_id in features:
                    raise ValueError(f"ID GFF duplicado: {feature_id}")
                features[feature_id] = {
                    "sample_id": sample,
                    "locus_id": feature_id,
                    "contig": fields[0],
                    "start": int(fields[3]),
                    "stop": int(fields[4]),
                    "strand": fields[6],
                    "gene": attributes.get("gene", ""),
                    "product": attributes.get("product", ""),
                }
    return features, samples


def coordinate_qc(
    source: str,
    sample_id: str,
    locus_id: str,
    reported_contig: str,
    reported_start: str,
    reported_stop: str,
    gff_by_id: dict[str, object] | None,
    mapping_method: str,
    assigned_gff: dict[str, object] | None,
    overlap_bp: int | str,
) -> dict[str, object]:
    found_by_id = gff_by_id is not None
    found = assigned_gff is not None
    start_value = int(reported_start) if reported_start else None
    stop_value = int(reported_stop) if reported_stop else None
    contig_match = found and reported_contig == assigned_gff["contig"]
    start_match = found and start_value == assigned_gff["start"]
    stop_match = found and stop_value == assigned_gff["stop"]
    all_match = found and contig_match and start_match and stop_match
    return {
        "fuente": source,
        "sample_id": sample_id,
        "locus_id_reportado": locus_id,
        "encontrado_por_id": "si" if found_by_id else "no",
        "encontrado_en_gff": "si" if found else "no",
        "metodo_mapeo": mapping_method,
        "locus_id_gff_asignado": "" if not found else assigned_gff["locus_id"],
        "solapamiento_gff_bp": overlap_bp,
        "contig_reportado": reported_contig,
        "start_reportado": "" if start_value is None else start_value,
        "stop_reportado": "" if stop_value is None else stop_value,
        "contig_gff": "" if not found else assigned_gff["contig"],
        "start_gff": "" if not found else assigned_gff["start"],
        "stop_gff": "" if not found else assigned_gff["stop"],
        "coordenadas_coinciden": "si" if all_match else "no",
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    gff_features, all_samples = load_gff_features()
    features_by_contig: dict[str, list[dict[str, object]]] = defaultdict(list)
    for feature in gff_features.values():
        features_by_contig[str(feature["contig"])].append(feature)

    amr_rows = [row for row in read_tsv(AMR_PATH) if row["Type"] == "AMR"]
    bacmet_rows = [
        row for row in read_tsv(BACMET_PATH) if row["is_main_mrg"].lower() == "si"
    ]

    qc_rows: list[dict[str, object]] = []
    args_by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)
    mrgs_by_sample: dict[str, list[dict[str, object]]] = defaultdict(list)

    for row in amr_rows:
        locus_id = row["Protein id"]
        gff = None if locus_id in ("", "NA") else gff_features.get(locus_id)
        assigned_gff = gff
        mapping_method = "id_exacto" if gff is not None else "sin_mapeo"
        overlap_bp: int | str = ""

        if assigned_gff is None:
            reported_start = int(row["Start"])
            reported_stop = int(row["Stop"])
            overlap_candidates: list[tuple[int, dict[str, object]]] = []
            for feature in features_by_contig.get(row["Contig id"], []):
                overlap = max(
                    0,
                    min(reported_stop, int(feature["stop"]))
                    - max(reported_start, int(feature["start"]))
                    + 1,
                )
                if overlap:
                    overlap_candidates.append((overlap, feature))
            if overlap_candidates:
                overlap_bp, assigned_gff = max(
                    overlap_candidates,
                    key=lambda item: (item[0], -int(item[1]["start"])),
                )
                mapping_method = "solapamiento_coordenadas"

        qc_rows.append(
            coordinate_qc(
                "AMRFinderPlus",
                row["sample_id"],
                locus_id,
                row["Contig id"],
                row["Start"],
                row["Stop"],
                gff,
                mapping_method,
                assigned_gff,
                overlap_bp,
            )
        )
        if gff is None and assigned_gff is None:
            continue
        if gff is None:
            coordinate_feature = {
                "sample_id": row["sample_id"],
                "locus_id": (
                    f"AMRFinderPlus|{row['sample_id']}|{row['Contig id']}|"
                    f"{row['Start']}-{row['Stop']}"
                ),
                "contig": row["Contig id"],
                "start": int(row["Start"]),
                "stop": int(row["Stop"]),
                "strand": row["Strand"],
                "gene": assigned_gff["gene"],
                "product": assigned_gff["product"],
            }
        else:
            coordinate_feature = gff
        args_by_sample[row["sample_id"]].append(
            {
                **coordinate_feature,
                "reported_locus_id": locus_id,
                "gff_locus_id": assigned_gff["locus_id"],
                "mapping_method": mapping_method,
                "symbol": row["Element symbol"],
                "name": row["Element name"],
                "class": row["Class"],
                "subclass": row["Subclass"],
                "method": row["Method"],
                "identity": row["% Identity to reference"],
                "coverage": row["% Coverage of reference"],
            }
        )

    for row in bacmet_rows:
        locus_id = row["locus_id"]
        gff = gff_features.get(locus_id)
        qc_rows.append(
            coordinate_qc(
                "BacMet",
                row["sample_id"],
                locus_id,
                row["contig"],
                row["start"],
                row["stop"],
                gff,
                "id_exacto" if gff is not None else "sin_mapeo",
                gff,
                "",
            )
        )
        if gff is None:
            continue
        mrgs_by_sample[row["sample_id"]].append(
            {
                **gff,
                "feature_label": row["feature_label"],
                "gene_candidates": row["gene_candidates"],
                "bacmet_ids": row["bacmet_ids"],
                "metal_list": row["metal_list"],
                "metal_category": row["metal_category"],
                "ambiguity_status": row["ambiguity_status"],
            }
        )

    pair_rows: list[dict[str, object]] = []
    for sample_id in sorted(set(args_by_sample) & set(mrgs_by_sample)):
        for arg in args_by_sample[sample_id]:
            for mrg in mrgs_by_sample[sample_id]:
                same_contig = arg["contig"] == mrg["contig"]
                same_locus = arg["gff_locus_id"] == mrg["locus_id"]

                distance: int | str = ""
                center_distance: float | str = ""
                overlap: int | str = ""
                order = "contigs_distintos"
                if same_contig:
                    overlap = max(
                        0,
                        min(int(arg["stop"]), int(mrg["stop"]))
                        - max(int(arg["start"]), int(mrg["start"]))
                        + 1,
                    )
                    if int(arg["stop"]) < int(mrg["start"]):
                        distance = int(mrg["start"]) - int(arg["stop"]) - 1
                        order = "ARG_antes_MRG"
                    elif int(mrg["stop"]) < int(arg["start"]):
                        distance = int(arg["start"]) - int(mrg["stop"]) - 1
                        order = "MRG_antes_ARG"
                    else:
                        distance = 0
                        order = "solapados"
                    arg_center = (int(arg["start"]) + int(arg["stop"])) / 2
                    mrg_center = (int(mrg["start"]) + int(mrg["stop"])) / 2
                    center_distance = round(abs(arg_center - mrg_center), 1)

                pair_rows.append(
                    {
                        "sample_id": sample_id,
                        "arg_locus_id": arg["locus_id"],
                        "arg_locus_id_reportado": arg["reported_locus_id"],
                        "arg_locus_id_gff": arg["gff_locus_id"],
                        "arg_metodo_mapeo_gff": arg["mapping_method"],
                        "arg_symbol": arg["symbol"],
                        "arg_name": arg["name"],
                        "arg_class": arg["class"],
                        "arg_subclass": arg["subclass"],
                        "arg_method": arg["method"],
                        "arg_identity_pct": arg["identity"],
                        "arg_coverage_pct": arg["coverage"],
                        "arg_contig": arg["contig"],
                        "arg_start": arg["start"],
                        "arg_stop": arg["stop"],
                        "arg_strand": arg["strand"],
                        "mrg_locus_id": mrg["locus_id"],
                        "mrg_feature_label": mrg["feature_label"],
                        "mrg_gene_candidates": mrg["gene_candidates"],
                        "bacmet_ids": mrg["bacmet_ids"],
                        "metal_list": mrg["metal_list"],
                        "metal_category": mrg["metal_category"],
                        "mrg_ambiguity_status": mrg["ambiguity_status"],
                        "mrg_contig": mrg["contig"],
                        "mrg_start": mrg["start"],
                        "mrg_stop": mrg["stop"],
                        "mrg_strand": mrg["strand"],
                        "mismo_contig": "si" if same_contig else "no",
                        "mismo_locus": "si" if same_locus else "no",
                        "orden_genomico": order,
                        "distancia_bordes_bp": distance,
                        "distancia_centros_bp": center_distance,
                        "solapamiento_bp": overlap,
                    }
                )

    pair_fields = list(pair_rows[0]) if pair_rows else []
    same_contig_rows = [row for row in pair_rows if row["mismo_contig"] == "si"]
    same_contig_rows.sort(
        key=lambda row: (
            row["sample_id"],
            int(row["distancia_bordes_bp"]),
            row["arg_locus_id"],
            row["mrg_locus_id"],
        )
    )
    high_confidence_rows = [
        row
        for row in same_contig_rows
        if row["arg_metodo_mapeo_gff"] == "id_exacto"
    ]

    summary_rows: list[dict[str, object]] = []
    for sample_id in sorted(all_samples):
        sample_pairs = [row for row in pair_rows if row["sample_id"] == sample_id]
        colocated = [row for row in sample_pairs if row["mismo_contig"] == "si"]
        distinct = [row for row in colocated if row["mismo_locus"] == "no"]
        nearest = min(distinct, key=lambda row: int(row["distancia_bordes_bp"]), default=None)
        exact_distinct = [
            row for row in distinct if row["arg_metodo_mapeo_gff"] == "id_exacto"
        ]
        nearest_exact = min(
            exact_distinct,
            key=lambda row: int(row["distancia_bordes_bp"]),
            default=None,
        )
        summary_rows.append(
            {
                "sample_id": sample_id,
                "n_arg": len(args_by_sample.get(sample_id, [])),
                "n_mrg": len(mrgs_by_sample.get(sample_id, [])),
                "n_pares_arg_mrg": len(sample_pairs),
                "n_pares_mismo_contig": len(colocated),
                "n_pares_mismo_locus": sum(row["mismo_locus"] == "si" for row in colocated),
                "n_pares_distintos_loci_hasta_10kb": sum(
                    row["mismo_locus"] == "no"
                    and int(row["distancia_bordes_bp"]) <= 10_000
                    for row in colocated
                ),
                "n_pares_distintos_loci_hasta_50kb": sum(
                    row["mismo_locus"] == "no"
                    and int(row["distancia_bordes_bp"]) <= 50_000
                    for row in colocated
                ),
                "n_pares_mismo_contig_alta_confianza": sum(
                    row["arg_metodo_mapeo_gff"] == "id_exacto" for row in colocated
                ),
                "n_pares_alta_confianza_hasta_10kb": sum(
                    int(row["distancia_bordes_bp"]) <= 10_000 for row in exact_distinct
                ),
                "min_distancia_distintos_loci_bp": (
                    "" if nearest is None else nearest["distancia_bordes_bp"]
                ),
                "arg_mas_cercano": "" if nearest is None else nearest["arg_symbol"],
                "mrg_mas_cercano": "" if nearest is None else nearest["mrg_feature_label"],
                "contig_par_mas_cercano": "" if nearest is None else nearest["arg_contig"],
                "min_distancia_alta_confianza_bp": (
                    "" if nearest_exact is None else nearest_exact["distancia_bordes_bp"]
                ),
                "arg_alta_confianza_mas_cercano": (
                    "" if nearest_exact is None else nearest_exact["arg_symbol"]
                ),
                "mrg_alta_confianza_mas_cercano": (
                    "" if nearest_exact is None else nearest_exact["mrg_feature_label"]
                ),
            }
        )

    write_tsv(OUTPUT_DIR / "pares_arg_mrg.tsv", pair_rows, pair_fields)
    write_tsv(
        OUTPUT_DIR / "pares_arg_mrg_mismo_contig.tsv",
        same_contig_rows,
        pair_fields,
    )
    write_tsv(
        OUTPUT_DIR / "pares_arg_mrg_mismo_contig_alta_confianza.tsv",
        high_confidence_rows,
        pair_fields,
    )
    write_tsv(
        OUTPUT_DIR / "resumen_colocalizacion_por_muestra.tsv",
        summary_rows,
        list(summary_rows[0]),
    )
    write_tsv(
        OUTPUT_DIR / "qc_coordenadas_arg_mrg.tsv",
        qc_rows,
        list(qc_rows[0]),
    )

    qc_matches = sum(row["coordenadas_coinciden"] == "si" for row in qc_rows)
    print(f"ARG usados: {sum(len(rows) for rows in args_by_sample.values())}")
    print(f"MRG usados: {sum(len(rows) for rows in mrgs_by_sample.values())}")
    print(f"Parejas ARG-MRG en la misma muestra: {len(pair_rows)}")
    print(f"Parejas en el mismo contig: {len(same_contig_rows)}")
    print(f"Coordenadas idénticas al GFF: {qc_matches}/{len(qc_rows)}")
    print(f"Resultados: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
