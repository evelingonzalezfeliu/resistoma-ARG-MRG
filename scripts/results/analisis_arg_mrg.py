#!/usr/bin/env python3
"""Depura e integra resultados AMRFinderPlus y BacMet/DIAMOND.

No modifica los insumos. Lee directamente los ZIP y genera tablas trazables.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
import os
import re
import sys
import zipfile
from collections import Counter, defaultdict
from datetime import date, datetime
from pathlib import Path
from urllib.parse import unquote

from openpyxl import load_workbook


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1]
TABLES = OUT / "tablas"
QC_DIR = OUT / "qc"

AMR_ZIP = ROOT / "resultados_amrfinder_combinado_plus.zip"
DIAMOND_ZIP = ROOT / "resultados_diamond_bacmet.zip"
BACMET_ZIP = ROOT / "bacmet_db.zip"
EXPECTED_FILE = ROOT / "lista.txt"
FAA_DIR = ROOT / "proteinas_faa"
GFF_DIR = ROOT / "gff"
FNA_DIR = ROOT / "secuencias_extraidas" / "fna"
METADATA_XLSX = ROOT.parent / "Genomas_Cepario_Bioing_UOH" / "CEPARIO 2025.xlsx"

DIAMOND_COLS = [
    "qseqid", "sseqid", "pident", "length", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qlen", "slen",
]

METAL_PATTERNS = {
    "Aluminio": [r"\baluminium\b", r"\baluminum\b", r"\bAl\b"],
    "Antimonio": [r"\bantimony\b", r"\bSb\b"],
    "Arsenico": [r"\barsenic\b", r"\barsenate\b", r"\barsenite\b", r"\bAs\b"],
    "Bismuto": [r"\bbismuth\b", r"\bBi\b"],
    "Cadmio": [r"\bcadmium\b", r"\bCd\b"],
    "Cobalto": [r"\bcobalt\b", r"\bCo\b"],
    "Cobre": [r"\bcopper\b", r"\bCu\b"],
    "Cromo": [r"\bchromium\b", r"\bchromate\b", r"\bCr\b"],
    "Hierro": [r"\biron\b", r"\bFe\b"],
    "Manganeso": [r"\bmanganese\b", r"\bMn\b"],
    "Mercurio": [r"\bmercury\b", r"\bmercuric\b", r"\bmercurous\b", r"\bHg\b"],
    "Molibdeno": [r"\bmolybdenum\b", r"\bmolybdate\b", r"\bMo\b"],
    "Niquel": [r"\bnickel\b", r"\bNi\b"],
    "Oro": [r"\bgold\b", r"\bAu\b"],
    "Plata": [r"\bsilver\b", r"\bAg\b"],
    "Plomo": [r"\blead\b", r"\bPb\b"],
    "Selenio": [r"\bselenium\b", r"\bselenate\b", r"\bselenite\b", r"\bSe\b"],
    "Telurio": [r"\btellurium\b", r"\btellurite\b", r"\bTe\b"],
    "Uranio": [r"\buranium\b", r"\buranyl\b", r"\bU\b"],
    "Vanadio": [r"\bvanadium\b", r"\bvanadate\b", r"\bV\b"],
    "Zinc": [r"\bzinc\b", r"\bZn\b"],
}


def ensure_dirs():
    for path in (TABLES, QC_DIR):
        path.mkdir(parents=True, exist_ok=True)


def clean(value):
    if value is None:
        return ""
    return str(value).strip()


def normalize_isolate(value):
    s = clean(value).upper().replace(" ", "").replace(".", "")
    m = re.fullmatch(r"B(\d)(\d)([A-J])(\d{2})", s)
    if m:
        return f"B{m.group(1)}{m.group(2)}_{m.group(3)}{m.group(4)}"
    return s


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def merge_intervals(intervals):
    valid = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = []
    for start, stop in valid:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, stop])
        else:
            merged[-1][1] = max(merged[-1][1], stop)
    return sum(stop - start + 1 for start, stop in merged)


def parse_faa_ids(path):
    ids = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def parse_fasta_lengths(path):
    lengths = {}
    current = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current:
                lengths[current] += len(line.strip())
    return lengths


def parse_gff(path):
    records = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            attrs = {}
            for item in parts[8].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attrs[key] = unquote(value)
            key = attrs.get("ID") or attrs.get("Name")
            if key:
                records[key] = {
                    "contig": parts[0], "start": parts[3], "stop": parts[4],
                    "strand": parts[6], "gff_gene": attrs.get("gene", ""),
                    "gff_product": attrs.get("product", ""),
                }
    return records


def metal_annotation(compound):
    metals = []
    for label, patterns in METAL_PATTERNS.items():
        if any(re.search(p, compound, flags=re.I) for p in patterns):
            metals.append(label)
    # Cada elemento separado por coma se considera un compuesto. Si alguno no
    # parece metálico, la referencia se clasifica como asociación mixta.
    entries = [x.strip() for x in compound.split(",") if x.strip()]
    nonmetal = []
    for entry in entries:
        if not any(re.search(p, entry, flags=re.I) for pats in METAL_PATTERNS.values() for p in pats):
            nonmetal.append(entry)
    if metals and nonmetal:
        ref_class = "metal_y_biocida"
    elif metals:
        ref_class = "metal"
    else:
        ref_class = "biocida_no_metal"
    metal_category = ""
    if len(metals) == 1:
        metal_category = metals[0]
    elif len(metals) > 1:
        metal_category = "Multimetal"
    return sorted(set(metals)), ref_class, metal_category, nonmetal


def informative_taxon(value):
    value = clean(value)
    if not value or value.lower() in {"na", "n/a", "-"}:
        return ""
    return value


def genus(value):
    value = informative_taxon(value)
    if not value or value.lower().startswith("unclassified"):
        return ""
    token = re.split(r"[\s/_]", value)[0]
    token = re.sub(r"[^A-Za-z]", "", token)
    aliases = {"Pseudoarthobacter": "Pseudarthrobacter"}
    return aliases.get(token, token)


def metadata():
    wb = load_workbook(METADATA_XLSX, data_only=True, read_only=False)
    ws = wb.active
    rows = []
    dates_outside_expected = []
    for row in ws.iter_rows(min_row=5, max_row=107, min_col=1, max_col=62):
        for cell in row:
            if isinstance(cell.value, (datetime, date)) and cell.column != 5:
                dates_outside_expected.append({
                    "celda": cell.coordinate,
                    "encabezado": clean(ws.cell(4, cell.column).value or ws.cell(3, cell.column).value),
                    "valor_interpretado": cell.value.isoformat(),
                    "formato_excel": cell.number_format,
                    "decision": "no_reconstruido",
                })
    for r in range(5, 108):
        source_id = clean(ws.cell(r, 1).value)
        if not source_id:
            continue
        norm = normalize_isolate(source_id)
        seq1 = informative_taxon(ws.cell(r, 12).value)
        seq2 = informative_taxon(ws.cell(r, 13).value)
        g1, g2 = genus(seq1), genus(seq2)
        discordant = bool(g1 and g2 and g1.lower() != g2.lower())
        if discordant:
            taxonomy_group = "Discordante"
            taxonomy_display = f"1a: {seq1}; 2a: {seq2}"
        else:
            taxonomy_group = g2 or g1 or "No clasificado"
            chosen = seq2 or seq1
            taxonomy_display = chosen or "No informado"
        genome_code = clean(ws.cell(r, 2).value)
        copper = "No evaluado/no informado"
        if genome_code in {"G12", "G19", "G27", "G28", "G32", "G39"}:
            copper = "Alta"
        elif genome_code in {"G17", "G53"}:
            copper = "Media"
        elif genome_code in {"G51", "G64"}:
            copper = "Baja"
        rows.append({
            "excel_row": r,
            "id_excel_original": source_id,
            "id_normalizado": norm,
            "codigo_genoma": genome_code,
            "procedencia_original": clean(ws.cell(r, 3).value),
            "procedencia_normalizada": clean(ws.cell(r, 3).value).title(),
            "taxonomia_16S": clean(ws.cell(r, 11).value),
            "taxonomia_genoma_1": seq1 or "No informado",
            "taxonomia_genoma_2": seq2 or "No informado",
            "taxonomia_discordante_genero": "si" if discordant else "no",
            "taxonomia_para_agrupacion": taxonomy_group,
            "taxonomia_mostrada": taxonomy_display,
            "fenotipo_cobre": copper,
        })
    return rows, dates_outside_expected


def read_bacmet_mapping(z):
    raw = z.read("BacMet2_EXP.753.mapping.txt").decode("utf-8-sig")
    mapping = {}
    mapping_by_id = defaultdict(list)
    row_count = 0
    for row in csv.DictReader(io.StringIO(raw), delimiter="\t"):
        row_count += 1
        row = {k: clean(v) for k, v in row.items()}
        metals, ref_class, metal_cat, nonmetal = metal_annotation(row["Compound"])
        row.update({
            "metal_list": ";".join(metals),
            "reference_class": ref_class,
            "metal_category": metal_cat,
            "nonmetal_compounds_detected": ";".join(nonmetal),
        })
        # La accession distingue las cinco reutilizaciones de BACxxxx. El nombre
        # del gen no es una clave fiable: tres encabezados FASTA difieren del mapping.
        key = (row["BacMet_ID"], row["Accession"])
        mapping[key] = row
        mapping_by_id[row["BacMet_ID"]].append(row)
    return mapping, mapping_by_id, row_count


def bacmet_subject_key(sseqid):
    parts = sseqid.split("|")
    bac_id = parts[0] if parts else ""
    gene = parts[1] if len(parts) > 1 else ""
    accession = ""
    # Encabezados típicos: BACxxxx|gene|sp|ACCESSION|... o BACxxxx|gene|ACCESSION|...
    if len(parts) > 3 and parts[2] in {"sp", "tr"}:
        accession = parts[3]
    elif len(parts) > 2:
        accession = parts[2]
    return bac_id, gene, accession


def fasta_ids_from_bytes(data):
    ids = []
    for line in data.decode("utf-8", "replace").splitlines():
        if line.startswith(">"):
            ids.append(line[1:].split()[0].split("|")[0])
    return ids


def main():
    ensure_dirs()
    expected = [x.strip() for x in EXPECTED_FILE.read_text().splitlines() if x.strip()]
    expected_set = set(expected)

    metadata_rows, date_anomalies = metadata()
    md_by_norm = defaultdict(list)
    for row in metadata_rows:
        md_by_norm[row["id_normalizado"]].append(row)
    id_map = []
    joined_md = {}
    for sample in expected:
        matches = md_by_norm.get(sample, [])
        if len(matches) == 1:
            status = "exacto_por_id_normalizado"
            selected = matches[0]
        elif len(matches) > 1:
            status = "ambiguo_id_duplicado"
            selected = None
        else:
            status = "sin_correspondencia"
            selected = None
        # B29_A02 queda sin unir: el candidato de fila 107 duplica B28_A02.
        candidate = ""
        if sample == "B29_A02":
            candidate = "Fila 107: B2.8. A.02 (posible error tipografico; no unido)"
        base = {
            "sample_id": sample,
            "mapping_status": status,
            "candidate_needing_review": candidate,
        }
        if selected:
            base.update(selected)
            joined_md[sample] = selected
        else:
            for key in metadata_rows[0]:
                base.setdefault(key, "")
            joined_md[sample] = {
                "procedencia_normalizada": "No informado",
                "taxonomia_para_agrupacion": "No clasificado",
                "taxonomia_mostrada": "No informado",
                "fenotipo_cobre": "No evaluado/no informado",
                "taxonomia_discordante_genero": "no",
            }
        id_map.append(base)
    extras = [r for r in metadata_rows if r["id_normalizado"] not in expected_set]
    for r in extras:
        row = {"sample_id": "", "mapping_status": "registro_excel_fuera_de_lista", "candidate_needing_review": ""}
        row.update(r)
        id_map.append(row)
    write_tsv(TABLES / "correspondencia_ids_metadatos.tsv", id_map)
    write_tsv(QC_DIR / "celdas_fecha_sospechosas_excel.tsv", date_anomalies)

    faa_sets, gff_records, fna_lengths = {}, {}, {}
    for sample in expected:
        faa_path = FAA_DIR / f"{sample}.faa"
        gff_path = GFF_DIR / f"{sample}.gff"
        fna_path = FNA_DIR / f"{sample}.fna"
        faa_sets[sample] = parse_faa_ids(faa_path) if faa_path.exists() else set()
        gff_records[sample] = parse_gff(gff_path) if gff_path.exists() else {}
        fna_lengths[sample] = parse_fasta_lengths(fna_path) if fna_path.exists() else {}

    qc = {s: {
        "sample_id": s,
        "en_lista": "si",
        "faa_presente": "si" if (FAA_DIR / f"{s}.faa").exists() else "no",
        "gff_presente": "si" if (GFF_DIR / f"{s}.gff").exists() else "no",
        "fna_presente": "si" if (FNA_DIR / f"{s}.fna").exists() else "no",
        "metadata_status": next(x["mapping_status"] for x in id_map if x["sample_id"] == s),
    } for s in expected}

    amr_rows = []
    amr_review = []
    with zipfile.ZipFile(AMR_ZIP) as z:
        names = set(z.namelist())
        for sample in expected:
            tsv_name = f"1057770/{sample}.amrfinder.tsv"
            log_name = f"1057770/{sample}.amrfinder.log"
            q = qc[sample]
            q["amrfinder_tsv_presente"] = "si" if tsv_name in names else "no"
            q["amrfinder_log_presente"] = "si" if log_name in names else "no"
            log_text = z.read(log_name).decode("utf-8", "replace") if log_name in names else ""
            completed = bool(re.search(r"amrfinder took \d+ seconds to complete", log_text))
            error = bool(re.search(r"(?i)\b(error|failed|killed|out of memory|oom)\b", log_text))
            q["amrfinder_log_completo"] = "si" if completed else "no"
            q["amrfinder_log_error"] = "si" if error else "no"
            rows = []
            if tsv_name in names:
                text = io.TextIOWrapper(z.open(tsv_name), encoding="utf-8-sig")
                rows = list(csv.DictReader(text, delimiter="\t"))
            q["amrfinder_hallazgos"] = len(rows)
            q["amrfinder_estado"] = (
                "fallido_o_incompleto" if not completed or error or tsv_name not in names
                else "correcto_sin_hallazgos" if not rows
                else "correcto_con_hallazgos"
            )
            missing_faa, missing_gff, missing_contig = set(), set(), set()
            locus_mismatch = set()
            invalid_nucleotide_coords = set()
            for row in rows:
                row["sample_id"] = sample
                pid = clean(row.get("Protein id"))
                if pid.lower() in {"na", "n/a", "-"}:
                    pid = ""
                if pid and pid not in faa_sets[sample]:
                    missing_faa.add(pid)
                if pid and pid not in gff_records[sample]:
                    missing_gff.add(pid)
                contig = clean(row.get("Contig id"))
                try:
                    start = int(row.get("Start") or 0); stop = int(row.get("Stop") or 0)
                except ValueError:
                    start = stop = 0
                if contig not in fna_lengths[sample]:
                    missing_contig.add(contig)
                elif not (1 <= min(start, stop) <= max(start, stop) <= fna_lengths[sample][contig]):
                    invalid_nucleotide_coords.add(f"{contig}:{start}-{stop}")
                if pid and pid in gff_records[sample]:
                    grec = gff_records[sample][pid]
                    if (contig, str(min(start, stop)), str(max(start, stop))) != (grec["contig"], grec["start"], grec["stop"]):
                        locus_mismatch.add(pid)
                ref_len = float(row.get("Reference sequence length") or 0)
                target_len = float(row.get("Target length") or 0)
                ratio = target_len / ref_len if ref_len else math.nan
                row["target_reference_length_ratio"] = "" if math.isnan(ratio) else f"{ratio:.4f}"
                flags = []
                if row.get("Method", "").startswith("PARTIAL") or row.get("Method") == "INTERNAL_STOP":
                    flags.append("coincidencia_parcial_o_stop")
                if not math.isnan(ratio) and (ratio < 0.8 or ratio > 1.2):
                    flags.append("longitud_objetivo_difiere_mas_20pct")
                if pid in locus_mismatch:
                    flags.append("coordenadas_no_coinciden_con_gff")
                if contig in missing_contig or f"{contig}:{start}-{stop}" in invalid_nucleotide_coords:
                    flags.append("contig_o_coordenadas_no_validos_en_fna")
                row["review_flags"] = ";".join(flags)
                if flags:
                    amr_review.append(row.copy())
                amr_rows.append(row)
            q["amrfinder_protein_ids_no_en_faa"] = len(missing_faa)
            q["amrfinder_protein_ids_no_en_gff"] = len(missing_gff)
            q["amrfinder_contig_ids_no_en_fna"] = len(missing_contig)
            q["amrfinder_coordenadas_fuera_fna"] = len(invalid_nucleotide_coords)
            q["amrfinder_locus_no_coincide_gff"] = len(locus_mismatch)

    if amr_rows:
        fields = list(amr_rows[0])
        write_tsv(TABLES / "amrfinderplus_depurado.tsv", amr_rows, fields)
        write_tsv(QC_DIR / "hallazgos_amrfinder_revision.tsv", amr_review, fields)

    bacmet_integrity = []
    with zipfile.ZipFile(BACMET_ZIP) as bz:
        mapping, mapping_by_id, mapping_row_count = read_bacmet_mapping(bz)
        fasta_bytes = bz.read("BacMet2_EXP_database.fasta")
        map_bytes = bz.read("BacMet2_EXP.753.mapping.txt")
        fasta_ids = fasta_ids_from_bytes(fasta_bytes)
        checksum_text = bz.read("checksums.sha256").decode()
        expected_checksums = dict(line.split(maxsplit=1)[::-1] for line in checksum_text.splitlines() if line.strip())
        bacmet_integrity = [
            {"control": "filas_mapping", "valor": mapping_row_count, "esperado": 753, "resultado": "OK" if mapping_row_count == 753 else "REVISAR"},
            {"control": "ids_bacmet_unicos_mapping", "valor": len(mapping_by_id), "esperado": "748 (5 IDs reutilizados)", "resultado": "OK" if len(mapping_by_id) == 748 else "REVISAR"},
            {"control": "secuencias_fasta", "valor": len(fasta_ids), "esperado": 753, "resultado": "OK" if len(fasta_ids) == 753 else "REVISAR"},
            {"control": "ids_mapping_no_en_fasta", "valor": len(set(mapping_by_id) - set(fasta_ids)), "esperado": 0, "resultado": "OK" if not set(mapping_by_id) - set(fasta_ids) else "REVISAR"},
            {"control": "ids_fasta_no_en_mapping", "valor": len(set(fasta_ids) - set(mapping_by_id)), "esperado": 0, "resultado": "OK" if not set(fasta_ids) - set(mapping_by_id) else "REVISAR"},
            {"control": "sha256_fasta", "valor": sha256_bytes(fasta_bytes), "esperado": expected_checksums.get("BacMet2_EXP_database.fasta", ""), "resultado": "OK" if sha256_bytes(fasta_bytes) == expected_checksums.get("BacMet2_EXP_database.fasta") else "REVISAR"},
            {"control": "sha256_mapping", "valor": sha256_bytes(map_bytes), "esperado": expected_checksums.get("BacMet2_EXP.753.mapping.txt", ""), "resultado": "OK" if sha256_bytes(map_bytes) == expected_checksums.get("BacMet2_EXP.753.mapping.txt") else "REVISAR"},
        ]
    write_tsv(QC_DIR / "integridad_bacmet_db.tsv", bacmet_integrity)

    # Agregación por pareja query-referencia, con coberturas por unión de intervalos.
    pair_data = {}
    diamond_file_ids = []
    with zipfile.ZipFile(DIAMOND_ZIP) as z:
        names = set(z.namelist())
        for sample in expected:
            name = f"{sample}.bacmet.raw.tsv"
            q = qc[sample]
            q["diamond_archivo_presente"] = "si" if name in names else "no"
            q["diamond_log_presente"] = "no"
            q["diamond_estado"] = "archivo_ausente"
            q["diamond_filas_raw"] = 0
            q["diamond_filas_invalidas"] = 0
            q["diamond_coordenadas_invalidas"] = 0
            q["diamond_pident_fuera_rango"] = 0
            q["diamond_qseqid_prefijo_incorrecto"] = 0
            q["diamond_qseqid_no_en_faa"] = 0
            if name not in names:
                continue
            diamond_file_ids.append(sample)
            invalid_qids = set()
            prefix_bad = set()
            rows_ok = 0
            with io.TextIOWrapper(z.open(name), encoding="utf-8") as handle:
                reader = csv.reader(handle, delimiter="\t")
                for line_no, vals in enumerate(reader, 1):
                    q["diamond_filas_raw"] += 1
                    if len(vals) != 12:
                        q["diamond_filas_invalidas"] += 1
                        continue
                    row = dict(zip(DIAMOND_COLS, vals))
                    try:
                        pident = float(row["pident"]); length = int(row["length"])
                        qstart = int(row["qstart"]); qend = int(row["qend"])
                        sstart = int(row["sstart"]); send = int(row["send"])
                        evalue = float(row["evalue"]); bitscore = float(row["bitscore"])
                        qlen = int(row["qlen"]); slen = int(row["slen"])
                    except ValueError:
                        q["diamond_filas_invalidas"] += 1
                        continue
                    coord_bad = not (1 <= qstart <= qlen and 1 <= qend <= qlen and 1 <= sstart <= slen and 1 <= send <= slen and length > 0)
                    if coord_bad:
                        q["diamond_coordenadas_invalidas"] += 1
                        continue
                    if not (0 <= pident <= 100):
                        q["diamond_pident_fuera_rango"] += 1
                        continue
                    if not row["qseqid"].startswith(sample + "|"):
                        prefix_bad.add(row["qseqid"])
                    if row["qseqid"] not in faa_sets[sample]:
                        invalid_qids.add(row["qseqid"])
                    bac_id, subject_gene, subject_accession = bacmet_subject_key(row["sseqid"])
                    key = (sample, row["qseqid"], row["sseqid"])
                    rec = pair_data.setdefault(key, {
                        "sample_id": sample, "qseqid": row["qseqid"], "sseqid": row["sseqid"],
                        "bacmet_id": bac_id, "subject_gene": subject_gene,
                        "subject_accession": subject_accession, "qlen": qlen, "slen": slen,
                        "q_intervals": [], "s_intervals": [], "hsps": 0,
                        "identity_num": 0.0, "identity_den": 0,
                        "pident_max": 0.0, "evalue_min": math.inf, "bitscore_max": -math.inf,
                    })
                    rec["q_intervals"].append((qstart, qend)); rec["s_intervals"].append((sstart, send))
                    rec["hsps"] += 1; rec["identity_num"] += pident * length; rec["identity_den"] += length
                    rec["pident_max"] = max(rec["pident_max"], pident)
                    rec["evalue_min"] = min(rec["evalue_min"], evalue)
                    rec["bitscore_max"] = max(rec["bitscore_max"], bitscore)
                    rows_ok += 1
            q["diamond_qseqid_prefijo_incorrecto"] = len(prefix_bad)
            q["diamond_qseqid_no_en_faa"] = len(invalid_qids)
            q["diamond_estado"] = "archivo_valido_sin_log_ejecucion" if rows_ok or q["diamond_filas_raw"] == 0 else "archivo_invalido"

    agg_rows = []
    for rec in pair_data.values():
        qcov = min(100.0, merge_intervals(rec["q_intervals"]) / rec["qlen"] * 100)
        scov = min(100.0, merge_intervals(rec["s_intervals"]) / rec["slen"] * 100)
        pident_weighted = rec["identity_num"] / rec["identity_den"]
        m = mapping.get((rec["bacmet_id"], rec["subject_accession"]), {})
        out = {
            "sample_id": rec["sample_id"], "qseqid": rec["qseqid"], "sseqid": rec["sseqid"],
            "bacmet_id": rec["bacmet_id"], "subject_gene": rec["subject_gene"],
            "subject_accession": rec["subject_accession"], "hsp_count": rec["hsps"],
            "pident_weighted": round(pident_weighted, 4), "pident_max": round(rec["pident_max"], 4),
            "query_coverage_pct": round(qcov, 4), "reference_coverage_pct": round(scov, 4),
            "evalue_min": rec["evalue_min"], "bitscore_max": rec["bitscore_max"],
            "qlen": rec["qlen"], "slen": rec["slen"],
            "Gene_name": m.get("Gene_name", ""), "Accession": m.get("Accession", ""),
            "Organism": m.get("Organism", ""), "Location": m.get("Location", ""),
            "Compound": m.get("Compound", ""), "metal_list": m.get("metal_list", ""),
            "reference_class": m.get("reference_class", "mapping_no_encontrado"),
            "metal_category": m.get("metal_category", ""),
            "mapping_found": "si" if m else "no",
        }
        for threshold in (70, 80, 90):
            out[f"pass_id{threshold}_cov80_e1e-5"] = "si" if (
                pident_weighted >= threshold and qcov >= 80 and scov >= 80 and rec["evalue_min"] <= 1e-5
            ) else "no"
        agg_rows.append(out)

    filtered = [r for r in agg_rows if r["pass_id80_cov80_e1e-5"] == "si"]
    by_query = defaultdict(list)
    for row in filtered:
        by_query[(row["sample_id"], row["qseqid"])].append(row)
    selected_groups = {}
    for key, rows in by_query.items():
        rows.sort(key=lambda r: (-float(r["bitscore_max"]), -float(r["query_coverage_pct"]), -float(r["reference_coverage_pct"]), -float(r["pident_weighted"]), float(r["evalue_min"]), r["bacmet_id"]))
        best = rows[0]
        rank_tuple = (best["bitscore_max"], best["query_coverage_pct"], best["reference_coverage_pct"], best["pident_weighted"], best["evalue_min"])
        tied = [r for r in rows if (r["bitscore_max"], r["query_coverage_pct"], r["reference_coverage_pct"], r["pident_weighted"], r["evalue_min"]) == rank_tuple]
        compatible = len({(r["Gene_name"], r["reference_class"], r["metal_list"]) for r in tied}) == 1
        status = "mejor_unico" if len(tied) == 1 else "empate_compatible" if compatible else "empate_incompatible"
        for idx, row in enumerate(rows, 1):
            row["candidate_rank"] = idx
            row["selected_top"] = "si" if row in tied else "no"
            row["top_tie_count"] = len(tied)
            row["ambiguity_status"] = status
        selected_groups[key] = tied
    filtered.sort(key=lambda r: (r["sample_id"], r["qseqid"], int(r["candidate_rank"])))
    filtered_fields = list(filtered[0]) if filtered else list(agg_rows[0]) + ["candidate_rank", "selected_top", "top_tie_count", "ambiguity_status"]
    write_tsv(TABLES / "bacmet_filtrado_anotado.tsv", filtered, filtered_fields)

    sensitivity = []
    for threshold in (70, 80, 90):
        pass_rows = [r for r in agg_rows if r[f"pass_id{threshold}_cov80_e1e-5"] == "si"]
        sensitivity.append({
            "identity_threshold_pct": threshold,
            "query_coverage_threshold_pct": 80,
            "reference_coverage_threshold_pct": 80,
            "evalue_max": "1e-5",
            "candidate_pairs": len(pass_rows),
            "proteins_with_candidates": len({(r["sample_id"], r["qseqid"]) for r in pass_rows}),
            "isolates_with_candidates": len({r["sample_id"] for r in pass_rows}),
            "metal_or_mixed_candidate_pairs": sum(r["reference_class"] in {"metal", "metal_y_biocida"} for r in pass_rows),
        })
    write_tsv(TABLES / "comparacion_filtros_bacmet.tsv", sensitivity)

    # Asignación final de una etiqueta por proteína; los empates incompatibles
    # permanecen como una sola etiqueta compuesta para no multiplicar genes.
    protein_assignments = []
    for (sample, qid), tied in sorted(selected_groups.items()):
        names = sorted({r["Gene_name"] or r["bacmet_id"] for r in tied})
        metals = sorted({m for r in tied for m in r["metal_list"].split(";") if m})
        classes = sorted({r["reference_class"] for r in tied})
        is_mrg = any(c in {"metal", "metal_y_biocida"} for c in classes)
        ambiguity = tied[0]["ambiguity_status"]
        label = names[0] if len(names) == 1 else "AMBIG:[" + "|".join(names) + "]"
        metal_category = metals[0] if len(metals) == 1 else "Multimetal" if len(metals) > 1 else ""
        gff = gff_records[sample].get(qid, {})
        protein_assignments.append({
            "sample_id": sample, "locus_id": qid, "feature_label": label,
            "gene_candidates": ";".join(names), "bacmet_ids": ";".join(sorted({r["bacmet_id"] for r in tied})),
            "ambiguity_status": ambiguity, "reference_classes": ";".join(classes),
            "metal_list": ";".join(metals), "metal_category": metal_category,
            "is_main_mrg": "si" if is_mrg else "no",
            "contig": gff.get("contig", ""), "start": gff.get("start", ""), "stop": gff.get("stop", ""),
            "gff_gene": gff.get("gff_gene", ""), "gff_product": gff.get("gff_product", ""),
        })
    write_tsv(TABLES / "bacmet_asignaciones_por_proteina.tsv", protein_assignments)

    # Integración por locus/proteína. ARG principales: AMRFinder Type=AMR.
    integrated = []
    amr_by_locus = defaultdict(list)
    for row in amr_rows:
        protein_id = clean(row.get("Protein id"))
        if protein_id.lower() in {"na", "n/a", "-"}:
            protein_id = ""
        locus = protein_id or f"{row.get('Contig id','')}:{row.get('Start','')}-{row.get('Stop','')}"
        amr_by_locus[(row["sample_id"], locus)].append(row)
    bac_by_locus = {(r["sample_id"], r["locus_id"]): r for r in protein_assignments}
    all_loci = sorted(set(amr_by_locus) | set(bac_by_locus))
    for sample, locus in all_loci:
        arows = amr_by_locus.get((sample, locus), [])
        brow = bac_by_locus.get((sample, locus))
        arg_rows = [r for r in arows if r.get("Type") == "AMR"]
        stress_rows = [r for r in arows if r.get("Type") == "STRESS" and r.get("Subtype") == "METAL"]
        if arg_rows:
            integrated.append({
                "sample_id": sample, "locus_id": locus, "analysis_category": "ARG_principal",
                "feature_label": ";".join(sorted({r["Element symbol"] for r in arg_rows})),
                "evidence_sources": "AMRFinderPlus", "amrfinder_type": "AMR",
                "amrfinder_methods": ";".join(sorted({r["Method"] for r in arg_rows})),
                "amrfinder_class": ";".join(sorted({r["Class"] for r in arg_rows})),
                "bacmet_ids": "", "metal_list": "", "ambiguity_status": "",
            })
        if brow and brow["is_main_mrg"] == "si":
            sources = ["BacMet_DIAMOND"]
            if stress_rows:
                sources.append("AMRFinderPlus_STRESS_METAL")
            integrated.append({
                "sample_id": sample, "locus_id": locus, "analysis_category": "MRG_principal",
                "feature_label": brow["feature_label"], "evidence_sources": ";".join(sources),
                "amrfinder_type": "STRESS/METAL" if stress_rows else "",
                "amrfinder_methods": ";".join(sorted({r["Method"] for r in stress_rows})),
                "amrfinder_class": ";".join(sorted({r["Class"] for r in stress_rows})),
                "bacmet_ids": brow["bacmet_ids"], "metal_list": brow["metal_list"],
                "metal_category": brow["metal_category"], "ambiguity_status": brow["ambiguity_status"],
            })
        elif stress_rows:
            integrated.append({
                "sample_id": sample, "locus_id": locus, "analysis_category": "MRG_complementario_AMRFinder",
                "feature_label": ";".join(sorted({r["Element symbol"] for r in stress_rows})),
                "evidence_sources": "AMRFinderPlus_STRESS_METAL", "amrfinder_type": "STRESS/METAL",
                "amrfinder_methods": ";".join(sorted({r["Method"] for r in stress_rows})),
                "amrfinder_class": ";".join(sorted({r["Class"] for r in stress_rows})),
                "bacmet_ids": "", "metal_list": "", "metal_category": "",
                "ambiguity_status": "sin_candidato_bacmet_principal_en_mismo_locus",
            })
    write_tsv(TABLES / "resultados_integrados.tsv", integrated)

    # Presencia/ausencia principal por aislado, no por copias.
    arg_presence = defaultdict(set)
    arg_meta = {}
    for row in amr_rows:
        if row.get("Type") == "AMR":
            label = row.get("Element symbol") or row.get("Element name")
            arg_presence[row["sample_id"]].add(label)
            arg_meta[label] = {"block": "ARG", "feature": label, "category": row.get("Class", "ARG")}
    mrg_presence = defaultdict(set)
    mrg_meta = {}
    for row in protein_assignments:
        if row["is_main_mrg"] == "si":
            label = row["feature_label"]
            mrg_presence[row["sample_id"]].add(label)
            old = mrg_meta.get(label)
            cat = row["metal_category"] or "Metal no resuelto"
            if old and old["category"] != cat:
                cat = "Multimetal"
            mrg_meta[label] = {"block": "MRG", "feature": label, "category": cat}
    features = list(arg_meta.values()) + list(mrg_meta.values())
    feature_order = [(r["block"], r["feature"]) for r in features]
    matrix_rows = []
    for sample in expected:
        md = joined_md[sample]
        row = {
            "sample_id": sample,
            "analysis_amrfinder": qc[sample]["amrfinder_estado"],
            "analysis_diamond": qc[sample]["diamond_estado"],
            "taxonomia": md["taxonomia_mostrada"],
            "taxonomia_grupo": md["taxonomia_para_agrupacion"],
            "procedencia": md["procedencia_normalizada"],
            "fenotipo_cobre": md["fenotipo_cobre"],
        }
        for block, feature in feature_order:
            present = feature in (arg_presence[sample] if block == "ARG" else mrg_presence[sample])
            available = qc[sample]["amrfinder_estado"].startswith("correcto") if block == "ARG" else qc[sample]["diamond_estado"].startswith("archivo_valido")
            row[f"{block}:{feature}"] = 1 if present else 0 if available else "NA"
        matrix_rows.append(row)
    matrix_fields = list(matrix_rows[0])
    write_tsv(TABLES / "matriz_presencia_ausencia.tsv", matrix_rows, matrix_fields)

    denominators = {
        "ARG": sum(qc[s]["amrfinder_estado"].startswith("correcto") for s in expected),
        "MRG": sum(qc[s]["diamond_estado"].startswith("archivo_valido") for s in expected),
    }
    frequencies = []
    for meta in features:
        block, feature = meta["block"], meta["feature"]
        numerator = sum(feature in (arg_presence[s] if block == "ARG" else mrg_presence[s]) for s in expected)
        denom = denominators[block]
        frequencies.append({
            "block": block, "feature": feature, "category": meta["category"],
            "positive_isolates_n": numerator, "analyzed_isolates_N": denom,
            "frequency_pct": round(numerator / denom * 100, 2) if denom else "",
            "definition": "aislados positivos / aislados analizados correctamente",
        })
    frequencies.sort(key=lambda r: (r["block"], -r["positive_isolates_n"], r["feature"]))
    write_tsv(TABLES / "frecuencias_arg_mrg.tsv", frequencies)

    # Resumen de concordancia de MRG por locus exacto (ID de proteína).
    concordance = []
    for sample in expected:
        bac_loci = {r["locus_id"] for r in protein_assignments if r["sample_id"] == sample and r["is_main_mrg"] == "si"}
        amr_loci = {locus for (s, locus), rows in amr_by_locus.items() if s == sample and any(r.get("Type") == "STRESS" and r.get("Subtype") == "METAL" for r in rows)}
        concordance.append({
            "sample_id": sample, "bacmet_mrg_loci": len(bac_loci),
            "amrfinder_stress_metal_loci": len(amr_loci),
            "shared_exact_protein_ids": len(bac_loci & amr_loci),
            "bacmet_only_loci": len(bac_loci - amr_loci),
            "amrfinder_only_loci": len(amr_loci - bac_loci),
        })
    write_tsv(TABLES / "concordancia_mrg_por_muestra.tsv", concordance)

    for sample in expected:
        q = qc[sample]
        q["bacmet_candidatos_filtrados"] = sum(r["sample_id"] == sample for r in filtered)
        q["bacmet_proteinas_asignadas"] = sum(r["sample_id"] == sample for r in protein_assignments)
        q["mrg_principales"] = len(mrg_presence[sample])
        q["arg_principales"] = len(arg_presence[sample])
        q["qc_global"] = "apto_con_limitacion_sin_log_diamond" if (
            q["amrfinder_estado"].startswith("correcto") and q["diamond_estado"].startswith("archivo_valido")
        ) else "revisar"
    write_tsv(TABLES / "qc_por_muestra.tsv", [qc[s] for s in expected])

    summary = {
        "expected_samples": len(expected),
        "amr_complete": sum(qc[s]["amrfinder_estado"].startswith("correcto") for s in expected),
        "amr_no_hits": sum(qc[s]["amrfinder_estado"] == "correcto_sin_hallazgos" for s in expected),
        "diamond_valid": sum(qc[s]["diamond_estado"].startswith("archivo_valido") for s in expected),
        "diamond_logs": sum(qc[s]["diamond_log_presente"] == "si" for s in expected),
        "amr_rows": len(amr_rows),
        "arg_rows": sum(r.get("Type") == "AMR" for r in amr_rows),
        "amr_stress_metal_rows": sum(r.get("Type") == "STRESS" and r.get("Subtype") == "METAL" for r in amr_rows),
        "bacmet_pairs_80": len(filtered),
        "bacmet_proteins_80": len(protein_assignments),
        "mrg_features": len(mrg_meta),
        "arg_features": len(arg_meta),
        "amr_review_rows": len(amr_review),
        "bacmet_incompatible_ties": sum(r["ambiguity_status"] == "empate_incompatible" for r in protein_assignments),
        "metadata_exact": sum(r["mapping_status"] == "exacto_por_id_normalizado" for r in id_map if r["sample_id"]),
        "taxonomy_discordant": sum(joined_md[s]["taxonomia_discordante_genero"] == "si" for s in expected),
        "excel_date_anomalies": len(date_anomalies),
    }
    write_tsv(QC_DIR / "resumen_ejecucion.tsv", [summary])
    print("\n".join(f"{k}\t{v}" for k, v in summary.items()))


if __name__ == "__main__":
    main()
