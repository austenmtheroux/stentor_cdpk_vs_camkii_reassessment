import sys
import csv
from pathlib import Path

if len(sys.argv) != 4:
    sys.exit(
        "Usage: python make_cdpk_vs_camkii_score_table.py"
        "<Stentor_vs_CDPK.tbl> <Stentor_vs_CaMKII.tbl> <output.tsv>"
    )

cdpk_tbl = Path(sys.argv[1])
camkii_tbl = Path(sys.argv[2])
out_tsv = Path(sys.argv[3])

def parse_tblout(path):
    """
    Parse HMMER --tblout output.

    Returns:
        dict mapping clean Stentor ID -> score/evalue information.

    Uses the full-sequence HMMER bit score and E-value because the target sequences are already PF00069 kinase-domain slices.
    """
    hits = {}

    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.split()

            raw_target = fields[0]
            full_seq_evalue = float(fields[4])
            full_seq_bitscore = float(fields[5])
            full_seq_bias = float(fields[6])

            clean_id = raw_target.split("|")[0]

            if clean_id not in hits or full_seq_bitscore > hits[clean_id]["bitscore"]:
                hits[clean_id] = {
                    "raw_target": raw_target,
                    "evalue": full_seq_evalue,
                    "bitscore": full_seq_bitscore,
                    "bias": full_seq_bias,
                }

    return hits

cdpk = parse_tblout(cdpk_tbl)
camkii = parse_tblout(camkii_tbl)

all_ids = sorted(set(cdpk) | set(camkii))

with open(out_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "Stentor_ID",
        "CDPK_bitscore",
        "CDPK_Evalue",
        "CaMKII_bitscore",
        "CaMKII_Evalue",
        "Delta_bitscore_CDPK_minus_CaMKII",
        "Preferred_model"
    ])

    for stentor_id in all_ids:
        cdpk_score = cdpk.get(stentor_id, {}).get("bitscore")
        cdpk_evalue = cdpk.get(stentor_id, {}).get("evalue")
        camkii_score = camkii.get(stentor_id, {}).get("bitscore")
        camkii_evalue = camkii.get(stentor_id, {}).get("evalue")

        if cdpk_score is None or camkii_score is None:
            delta = "NA"
            preferred = "missing_score"
        else:
            delta_value = cdpk_score - camkii_score
            delta = f"{delta_value:.2f}"
            preferred = "CDPK-like" if delta_value > 0 else "CaMKII-like"

        writer.writerow([
            stentor_id,
            "NA" if cdpk_score is None else f"{cdpk_score:.2f}",
            "NA" if cdpk_evalue is None else f"{cdpk_evalue:.2e}",
            "NA" if camkii_score is None else f"{camkii_score:.2f}",
            "NA" if camkii_evalue is None else f"{camkii_evalue:.2e}",
            delta,
            preferred,
        ])