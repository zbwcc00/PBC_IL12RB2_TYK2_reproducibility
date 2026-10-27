from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import requests


REPO = Path(__file__).resolve().parents[2]
PROCESSED_DIR = REPO / "work" / "genetics" / "regulatory_annotation" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

VARIANTS = [
    {
        "rsid": "rs6679356",
        "component": "L1",
        "pbc_pip": 0.940257624753562,
        "grch37_pos1": 67820194,
        "gwas_effect_allele": "C",
        "gwas_other_allele": "T",
        "gwas_beta": 0.439362,
        "gwas_p": 6.61382e-63,
        "grch38_variant_id": "1_67354511_C_T",
    },
    {
        "rsid": "rs6702599",
        "component": "L1",
        "pbc_pip": 0.059742375246438,
        "grch37_pos1": 67825399,
        "gwas_effect_allele": "A",
        "gwas_other_allele": "C",
        "gwas_beta": 0.433251,
        "gwas_p": 4.74331e-62,
        "grch38_variant_id": "1_67359716_A_C",
    },
    {
        "rsid": "rs3790569",
        "component": "L2",
        "pbc_pip": 0.260844507774078,
        "grch37_pos1": 67836860,
        "gwas_effect_allele": "C",
        "gwas_other_allele": "T",
        "gwas_beta": -0.246527,
        "gwas_p": 2.06337e-16,
        "grch38_variant_id": "1_67371177_T_C",
    },
    {
        "rsid": "rs3828069",
        "component": "L2",
        "pbc_pip": 0.575115492225937,
        "grch37_pos1": 67839573,
        "gwas_effect_allele": "C",
        "gwas_other_allele": "T",
        "gwas_beta": -0.249879,
        "gwas_p": 6.51302e-17,
        "grch38_variant_id": "1_67373890_T_C",
    },
    {
        "rsid": "rs11209064",
        "component": "L2",
        "pbc_pip": 0.161479991291523,
        "grch37_pos1": 67865853,
        "gwas_effect_allele": "T",
        "gwas_other_allele": "G",
        "gwas_beta": -0.254964,
        "gwas_p": 4.91551e-16,
        "grch38_variant_id": "1_67400170_G_T",
    },
    {
        "rsid": "rs72932868",
        "component": "L3",
        "pbc_pip": 0.137106139801697,
        "grch37_pos1": 67868626,
        "gwas_effect_allele": "A",
        "gwas_other_allele": "C",
        "gwas_beta": 0.0665103,
        "gwas_p": 0.107487552570303,
        "grch38_variant_id": "1_67402943_C_A",
    },
    {
        "rsid": "rs112588573",
        "component": "L3",
        "pbc_pip": 0.817829860198303,
        "grch37_pos1": 67888931,
        "gwas_effect_allele": "A",
        "gwas_other_allele": "G",
        "gwas_beta": 0.0553225,
        "gwas_p": 0.196306495822297,
        "grch38_variant_id": "1_67423248_G_A",
    },
]

OPEN_TARGETS_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
OPEN_TARGETS_QUERY = """
query v($id: String!) {
  variant(variantId: $id) {
    id
    rsIds
    mostSevereConsequence { id label }
    variantEffect {
      method
      score
      normalisedScore
      assessment
      target { id approvedSymbol }
    }
    enhancerToGenes(page: {index: 0, size: 200}) {
      count
      rows {
        datasourceId
        biosampleName
        biosampleFromSourceId
        distanceToTss
        score
        intervalType
        chromosome
        start
        end
        studyId
        pmid
        qualityControls
        resourceScore { name value }
        target { id approvedSymbol }
        biosample { biosampleId biosampleName description }
      }
    }
  }
}
"""


def request_json(method: str, url: str, **kwargs):
    response = requests.request(method, url, timeout=40, **kwargs)
    response.raise_for_status()
    return response.json()


def ncbi_mapping(rsid: str):
    refsnp_id = rsid.removeprefix("rs")
    data = request_json(
        "GET", f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{refsnp_id}"
    )
    mappings = {}
    for placement in data["primary_snapshot_data"]["placements_with_allele"]:
        assemblies = [
            row.get("assembly_name")
            for row in placement["placement_annot"].get(
                "seq_id_traits_by_assembly", []
            )
        ]
        for assembly in assemblies:
            if assembly not in {"GRCh37.p13", "GRCh38.p14"}:
                continue
            alleles = [row["allele"]["spdi"] for row in placement["alleles"]]
            mappings[assembly] = {
                "seq_id": alleles[0]["seq_id"],
                "pos1": int(alleles[0]["position"]) + 1,
                "reference": alleles[0]["deleted_sequence"],
                "alleles": "/".join(row["inserted_sequence"] for row in alleles),
            }
    return mappings


def ucsc_track(track: str, pos1: int):
    start0 = pos1 - 1
    params = {
        "genome": "hg38",
        "track": track,
        "chrom": "chr1",
        "start": start0,
        "end": pos1,
    }
    data = request_json(
        "GET", "https://api.genome.ucsc.edu/getData/track", params=params
    )
    return data.get(track, [])


def open_targets_variant(variant_id: str):
    payload = {"query": OPEN_TARGETS_QUERY, "variables": {"id": variant_id}}
    data = request_json("POST", OPEN_TARGETS_ENDPOINT, json=payload)
    if data.get("errors"):
        raise RuntimeError(json.dumps(data["errors"]))
    return data["data"]["variant"]


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    variant_rows = []
    e2g_rows = []
    for variant in VARIANTS:
        mapping = ncbi_mapping(variant["rsid"])
        grch38_pos1 = mapping["GRCh38.p14"]["pos1"]
        expected_pos1 = int(variant["grch38_variant_id"].split("_")[1])
        if grch38_pos1 != expected_pos1:
            raise RuntimeError(
                f"GRCh38 mapping mismatch for {variant['rsid']}: "
                f"{grch38_pos1} != {expected_pos1}"
            )

        ccres = ucsc_track("encodeCcreCombined", grch38_pos1)
        tf_clusters = ucsc_track("encRegTfbsClustered", grch38_pos1)
        ot_variant = open_targets_variant(variant["grch38_variant_id"])
        enhancer_rows = ot_variant["enhancerToGenes"]["rows"]
        il12rb2_rows = [
            row
            for row in enhancer_rows
            if row.get("target", {}).get("approvedSymbol") == "IL12RB2"
        ]
        t_cell_rows = [
            row for row in il12rb2_rows if "T cell" in row.get("biosampleName", "")
        ]

        for row in il12rb2_rows:
            e2g_rows.append(
                {
                    "rsid": variant["rsid"],
                    "component": variant["component"],
                    "grch38_variant_id": variant["grch38_variant_id"],
                    "datasource_id": row.get("datasourceId"),
                    "biosample_name": row.get("biosampleName"),
                    "distance_to_tss": row.get("distanceToTss"),
                    "score": row.get("score"),
                    "interval_type": row.get("intervalType"),
                    "chromosome": row.get("chromosome"),
                    "start": row.get("start"),
                    "end": row.get("end"),
                    "study_id": row.get("studyId"),
                    "pmid": row.get("pmid"),
                    "target": row.get("target", {}).get("approvedSymbol"),
                }
            )

        variant_rows.append(
            {
                **variant,
                "ncbi_grch37_pos1": mapping["GRCh37.p13"]["pos1"],
                "ncbi_grch38_pos1": grch38_pos1,
                "ncbi_grch38_reference": mapping["GRCh38.p14"]["reference"],
                "ncbi_grch38_alleles": mapping["GRCh38.p14"]["alleles"],
                "most_severe_consequence": ot_variant["mostSevereConsequence"][
                    "label"
                ],
                "encode_ccre_count": len(ccres),
                "encode_ccre_ids": ";".join(row.get("name", "") for row in ccres),
                "encode_ccre_classes": ";".join(
                    row.get("ccre", "") for row in ccres
                ),
                "encode_tf_cluster_count": len(tf_clusters),
                "encode_tf_names": ";".join(
                    row.get("name", "") for row in tf_clusters
                ),
                "e2g_total_count": len(enhancer_rows),
                "e2g_il12rb2_count": len(il12rb2_rows),
                "e2g_il12rb2_t_cell_count": len(t_cell_rows),
                "e2g_il12rb2_biosamples": ";".join(
                    sorted({row.get("biosampleName", "") for row in il12rb2_rows})
                ),
                "e2g_il12rb2_score_min": min(
                    (row["score"] for row in il12rb2_rows), default=""
                ),
                "e2g_il12rb2_score_max": max(
                    (row["score"] for row in il12rb2_rows), default=""
                ),
            }
        )

    variant_fields = list(variant_rows[0])
    e2g_fields = list(e2g_rows[0]) if e2g_rows else []
    write_tsv(PROCESSED_DIR / "01_credible_set_variant_regulatory_annotation.tsv", variant_rows, variant_fields)
    if e2g_rows:
        write_tsv(PROCESSED_DIR / "02_IL12RB2_E2G_links.tsv", e2g_rows, e2g_fields)

    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_conventions": {
            "NCBI_SPDI": "0-based position converted to 1-based POS in outputs",
            "UCSC_tracks": "0-based half-open query interval [POS-1, POS)",
            "Open_Targets_variant_id": "GRCh38 1-based chromosome_position_ref_alt",
        },
        "sources": {
            "NCBI_Variation": "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/",
            "UCSC_ENCODE_cCRE": "https://api.genome.ucsc.edu/getData/track?track=encodeCcreCombined",
            "UCSC_ENCODE_TF_cluster": "https://api.genome.ucsc.edu/getData/track?track=encRegTfbsClustered",
            "Open_Targets_Platform": OPEN_TARGETS_ENDPOINT,
        },
        "limitations": [
            "ENCODE cCRE and TF-cluster tracks are aggregated annotations and are not NK/T-cell specific.",
            "Open Targets E2G links are predictive enhancer-gene maps and do not establish allele-specific function.",
            "Absence of an annotation is not evidence that a variant lacks regulatory activity.",
        ],
    }
    (PROCESSED_DIR / "03_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
