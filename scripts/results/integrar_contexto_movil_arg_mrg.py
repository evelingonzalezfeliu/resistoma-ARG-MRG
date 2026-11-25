#!/usr/bin/env python3
"""Integra colocalización ARG-MRG y evidencia de elementos móviles.

La fila representativa de cada aislado es la pareja ARG-MRG de alta confianza
con menor distancia entre CDS. La salida conserva los campos de evidencia para
que las categorías del heatmap sean auditables. Se integran MOB-suite, ISEScan,
MobileElementFinder e IntegronFinder.
"""

from __future__ import annotations

import csv
from pathlib import Path


BASE = Path("/Users/efeliu/Documents/Mlatorre/analisis")
ANALYSIS = BASE / "analisis_arg_mrg_20261022"
PAIRS = ANALYSIS / "colocalizacion/pares_arg_mrg_mismo_contig_alta_confianza.tsv"
MATRIX = ANALYSIS / "tablas/matriz_presencia_ausencia.tsv"
MOB_ROOT = BASE / "resultados_mobsuite"
IS_ROOT = BASE / "resultados_isescan"
MEF_ROOT = BASE / "resultados_mobileelementfinder"
INTEGRON_ROOT = BASE / "resultados_integronfinder"
OUTPUT = ANALYSIS / "colocalizacion/contexto_movil_por_muestra.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_commented_csv(path: Path) -> list[dict[str, str]]:
    """Lee CSV con metadatos iniciales comentados mediante '#'."""
    with path.open(newline="", encoding="utf-8") as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        return list(csv.DictReader(lines))


def integer(row: dict[str, str], field: str) -> int:
    return int(float(row[field]))


def interval_distance(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    if b_end < a_start:
        return a_start - b_end - 1
    if b_start > a_end:
        return b_start - a_end - 1
    return 0


samples = [row["sample_id"] for row in read_tsv(MATRIX)]
pairs = read_tsv(PAIRS)

closest_pair: dict[str, dict[str, str]] = {}
for pair in pairs:
    sample = pair["sample_id"]
    distance = integer(pair, "distancia_bordes_bp")
    if sample not in closest_pair or distance < integer(
        closest_pair[sample], "distancia_bordes_bp"
    ):
        closest_pair[sample] = pair


def contig_mob_record(sample: str, contig: str) -> dict[str, str] | None:
    path = MOB_ROOT / sample / "contig_report.txt"
    if not path.exists():
        return None
    for row in read_tsv(path):
        reported = row.get("contig_id", "").split(" ", 1)[0]
        if reported == contig:
            return row
    return None


def plasmid_mobility(sample: str, cluster: str) -> str:
    path = MOB_ROOT / sample / "mobtyper_results.txt"
    if not path.exists() or cluster in ("", "-"):
        return "indeterminada"
    for row in read_tsv(path):
        if row.get("primary_cluster_id") == cluster:
            value = row.get("predicted_mobility", "").strip()
            return value if value and value != "-" else "indeterminada"
    return "indeterminada"


def nearest_is(
    sample: str,
    contig: str,
    arg_start: int,
    arg_end: int,
    mrg_start: int,
    mrg_end: int,
) -> tuple[dict[str, str] | None, int | None]:
    path = IS_ROOT / sample / "fna" / f"{sample}.fna.tsv"
    if not path.exists():
        return None, None
    best: dict[str, str] | None = None
    best_distance: int | None = None
    for row in read_tsv(path):
        if row.get("seqID") != contig:
            continue
        is_start = integer(row, "isBegin")
        is_end = integer(row, "isEnd")
        distance = min(
            interval_distance(arg_start, arg_end, is_start, is_end),
            interval_distance(mrg_start, mrg_end, is_start, is_end),
        )
        if best_distance is None or distance < best_distance:
            best, best_distance = row, distance
    return best, best_distance


def load_integron_regions(sample: str) -> list[dict[str, object]] | None:
    """Devuelve regiones de integrón o None si la muestra no fue analizada."""
    paths = sorted((INTEGRON_ROOT / sample).glob("**/*.integrons"))
    if not paths:
        return None
    groups: dict[tuple[str, str], dict[str, object]] = {}
    with paths[0].open(newline="", encoding="utf-8") as handle:
        lines = (line for line in handle if line.strip() and not line.startswith("#"))
        for row in csv.DictReader(lines, delimiter="\t"):
            key = (row["ID_replicon"], row["ID_integron"])
            start = int(row["pos_beg"])
            end = int(row["pos_end"])
            region = groups.setdefault(
                key,
                {
                    "contig": row["ID_replicon"],
                    "id": row["ID_integron"],
                    "type": row["type"],
                    "start": start,
                    "end": end,
                },
            )
            region["start"] = min(int(region["start"]), start)
            region["end"] = max(int(region["end"]), end)
    return list(groups.values())


def load_mef_regions(sample: str) -> list[dict[str, object]] | None:
    """Devuelve regiones MGE o None si falta la salida final de la muestra."""
    path = MEF_ROOT / sample / f"{sample}.csv"
    if not path.exists():
        return None
    regions: list[dict[str, object]] = []
    for row in read_commented_csv(path):
        start_text = row.get("start", "")
        end_text = row.get("end", row.get("stop", ""))
        if not start_text or not end_text:
            continue
        regions.append(
            {
                "contig": row.get("contig", "").split(" ", 1)[0],
                "id": row.get("mge_no", ""),
                "name": row.get("name", ""),
                "type": row.get("type", ""),
                "start": int(start_text),
                "end": int(end_text),
            }
        )
    return regions


def sample_integron_status(
    regions: list[dict[str, object]] | None,
) -> tuple[str, str, int | str]:
    if regions is None:
        return "No analizado", "NA", "NA"
    types = sorted({str(region["type"]) for region in regions})
    if "complete" in types:
        status = "Integron completo"
    elif "CALIN" in types:
        status = "CALIN"
    elif "In0" in types:
        status = "In0"
    else:
        status = "Sin integron"
    return status, ";".join(types) or "NA", len(regions)


def sample_mge_status(
    regions: list[dict[str, object]] | None,
) -> tuple[str, str, int | str]:
    if regions is None:
        return "No analizado", "NA", "NA"
    types = sorted({str(region["type"]) for region in regions})
    status = "MGE detectado" if regions else "Sin MGE"
    return status, ";".join(types) or "NA", len(regions)


def region_context(
    regions: list[dict[str, object]] | None,
    contig: str,
    arg_start: int,
    arg_end: int,
    mrg_start: int,
    mrg_end: int,
    label: str,
) -> dict[str, object]:
    """Resume solapamiento y distancia de una pareja a regiones genómicas."""
    if regions is None:
        return {
            "context": "No analizado",
            "nearest": None,
            "distance": None,
            "arg_ids": [],
            "mrg_ids": [],
            "shared_ids": [],
        }

    candidates = [region for region in regions if region["contig"] == contig]
    arg_ids: list[str] = []
    mrg_ids: list[str] = []
    best: dict[str, object] | None = None
    best_distance: int | None = None
    for region in candidates:
        start = int(region["start"])
        end = int(region["end"])
        arg_distance = interval_distance(arg_start, arg_end, start, end)
        mrg_distance = interval_distance(mrg_start, mrg_end, start, end)
        distance = min(arg_distance, mrg_distance)
        if arg_distance == 0:
            arg_ids.append(str(region["id"]))
        if mrg_distance == 0:
            mrg_ids.append(str(region["id"]))
        if best_distance is None or distance < best_distance:
            best = region
            best_distance = distance

    shared_ids = sorted(set(arg_ids) & set(mrg_ids))
    if shared_ids:
        context = f"ARG y MRG en {label}"
    elif arg_ids or mrg_ids:
        context = f"ARG o MRG en {label}"
    elif best_distance is not None and best_distance <= 10_000:
        context = f"{label} <=10 kb"
    elif best_distance is not None and best_distance <= 50_000:
        context = f"{label} 10-50 kb"
    else:
        context = f"Sin {label} <=50 kb"

    return {
        "context": context,
        "nearest": best,
        "distance": best_distance,
        "arg_ids": sorted(set(arg_ids)),
        "mrg_ids": sorted(set(mrg_ids)),
        "shared_ids": shared_ids,
    }


rows: list[dict[str, object]] = []
for sample in samples:
    sample_integrons = load_integron_regions(sample)
    sample_mges = load_mef_regions(sample)
    integron_status, integron_types, integron_count = sample_integron_status(
        sample_integrons
    )
    mge_status, mge_types, mge_count = sample_mge_status(sample_mges)
    pair = closest_pair.get(sample)
    if pair is None:
        rows.append(
            {
                "sample_id": sample,
                "arg_symbol": "NA",
                "mrg_symbol": "NA",
                "arg_mrg_contig": "NA",
                "arg_mrg_distance_bp": "NA",
                "arg_mrg_distance_class": "Sin pareja confiable",
                "molecule_type": "NA",
                "genomic_context": "Sin pareja confiable",
                "plasmid_cluster": "NA",
                "plasmid_mobility": "NA",
                "nearest_is_family": "NA",
                "nearest_is_type": "NA",
                "nearest_is_start": "NA",
                "nearest_is_end": "NA",
                "nearest_is_distance_bp": "NA",
                "is_proximity": "Sin pareja confiable",
                "integron_context": "Sin pareja confiable",
                "integron_status": integron_status,
                "integron_types": integron_types,
                "integron_count": integron_count,
                "arg_integron_ids": "NA",
                "mrg_integron_ids": "NA",
                "shared_integron_ids": "NA",
                "nearest_integron_id": "NA",
                "nearest_integron_type": "NA",
                "nearest_integron_distance_bp": "NA",
                "mef_context": "Sin pareja confiable",
                "mge_status": mge_status,
                "mge_types": mge_types,
                "mge_count": mge_count,
                "arg_mge_ids": "NA",
                "mrg_mge_ids": "NA",
                "shared_mge_ids": "NA",
                "nearest_mge_name": "NA",
                "nearest_mge_type": "NA",
                "nearest_mge_distance_bp": "NA",
                "transposon_status": "Sin pareja confiable",
            }
        )
        continue

    distance = integer(pair, "distancia_bordes_bp")
    if distance <= 10_000:
        distance_class = "<=10 kb"
    elif distance <= 50_000:
        distance_class = "10-50 kb"
    else:
        distance_class = ">50 kb"

    contig = pair["arg_contig"]
    arg_start = integer(pair, "arg_start")
    arg_end = integer(pair, "arg_stop")
    mrg_start = integer(pair, "mrg_start")
    mrg_end = integer(pair, "mrg_stop")
    mob = contig_mob_record(sample, contig)
    molecule = mob.get("molecule_type", "NA") if mob else "NA"
    cluster = mob.get("primary_cluster_id", "NA") if mob else "NA"
    mobility = (
        plasmid_mobility(sample, cluster) if molecule == "plasmid" else "NA"
    )

    is_row, is_distance = nearest_is(
        sample,
        contig,
        arg_start,
        arg_end,
        mrg_start,
        mrg_end,
    )
    if molecule == "plasmid":
        context = "Plasmido"
    elif molecule == "chromosome":
        context = "Cromosoma"
    else:
        context = "No clasificado"

    if is_distance is not None and is_distance <= 10_000:
        is_proximity = "IS <=10 kb"
    elif is_distance is not None and is_distance <= 50_000:
        is_proximity = "IS 10-50 kb"
    else:
        is_proximity = "Sin IS <=50 kb"

    integron = region_context(
        sample_integrons,
        contig,
        arg_start,
        arg_end,
        mrg_start,
        mrg_end,
        "integron",
    )
    mef = region_context(
        sample_mges,
        contig,
        arg_start,
        arg_end,
        mrg_start,
        mrg_end,
        "MGE",
    )
    nearest_integron = integron["nearest"]
    nearest_mge = mef["nearest"]

    rows.append(
        {
            "sample_id": sample,
            "arg_symbol": pair["arg_symbol"],
            "mrg_symbol": pair["mrg_feature_label"],
            "arg_mrg_contig": contig,
            "arg_mrg_distance_bp": distance,
            "arg_mrg_distance_class": distance_class,
            "molecule_type": molecule,
            "genomic_context": context,
            "plasmid_cluster": cluster,
            "plasmid_mobility": mobility,
            "nearest_is_family": is_row.get("family", "NA") if is_row else "NA",
            "nearest_is_type": is_row.get("type", "NA") if is_row else "NA",
            "nearest_is_start": is_row.get("isBegin", "NA") if is_row else "NA",
            "nearest_is_end": is_row.get("isEnd", "NA") if is_row else "NA",
            "nearest_is_distance_bp": is_distance if is_distance is not None else "NA",
            "is_proximity": is_proximity,
            "integron_context": integron["context"],
            "integron_status": integron_status,
            "integron_types": integron_types,
            "integron_count": integron_count,
            "arg_integron_ids": ";".join(integron["arg_ids"]) or "NA",
            "mrg_integron_ids": ";".join(integron["mrg_ids"]) or "NA",
            "shared_integron_ids": ";".join(integron["shared_ids"]) or "NA",
            "nearest_integron_id": nearest_integron["id"] if nearest_integron else "NA",
            "nearest_integron_type": nearest_integron["type"] if nearest_integron else "NA",
            "nearest_integron_distance_bp": integron["distance"] if integron["distance"] is not None else "NA",
            "mef_context": mef["context"],
            "mge_status": mge_status,
            "mge_types": mge_types,
            "mge_count": mge_count,
            "arg_mge_ids": ";".join(mef["arg_ids"]) or "NA",
            "mrg_mge_ids": ";".join(mef["mrg_ids"]) or "NA",
            "shared_mge_ids": ";".join(mef["shared_ids"]) or "NA",
            "nearest_mge_name": nearest_mge["name"] if nearest_mge else "NA",
            "nearest_mge_type": nearest_mge["type"] if nearest_mge else "NA",
            "nearest_mge_distance_bp": mef["distance"] if mef["distance"] is not None else "NA",
            "transposon_status": mef["context"],
        }
    )

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Muestras: {len(rows)}")
print(f"Con pareja confiable: {sum(r['arg_mrg_contig'] != 'NA' for r in rows)}")
print(f"Salida: {OUTPUT}")
