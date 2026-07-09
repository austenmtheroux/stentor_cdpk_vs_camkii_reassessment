"""
Generate the CDPK-vs-CaMKII reannotation SUMMARY TABLE for the Stentor coeruleus
proteins previously described as "CaMKII homologs".

This is the single analysis script for the paper. It parses the candidate FASTA and
the per-protein InterProScan TSVs and writes ONE table: the reannotation summary.
It deliberately does NOT produce the raw-domain-calls table, the architecture-summary
table, a report, or the domain-architecture schematic (the schematic is made by hand).

The emitted table is the AUTO (pre-curation) version. Downstream, it is curated by
hand to add:
  * `new_gene_id_if_mapped` — the mapping to the manuscript genome annotation, obtained
    by BLASTp of the SteCoe_ sequences against the newly annotated proteome;
  * a free-text `Notes` column.
Keep the AUTO output and the curated table as separate files (see README).

Inputs  (edit the paths below to match your layout):
  data/stentor_CaMKII_homolog_candidates_full_length.fasta   # 14 candidate sequences
  data/interpro_files/SteCoe_<ID>_interpro.tsv               # one InterProScan TSV each
Output:
  results/reannotation_summary_AUTO.tsv

Design notes:
  * Nothing is invented. Coordinates, counts, and family calls come only from the
    provided annotations. Where a signature is absent it is reported as
    "not detected in provided annotations", never as an absolute biological absence.
  * The reannotation call is driven by architecture, not by any single DB name:
    kinase catalytic domain + own C-terminal EF-hand/CaM-like lobe + CDPK PANTHER/
    InterPro family + NO CaMKII association/hub domain  ->  CDPK-like, not CaMKII.
  * Confidence and per-protein Notes are computed internally (Revised_annotation
    depends on the confidence tier) but are NOT written to the table; they are printed
    to the terminal as a QC aid for manual curation.
"""

import os
import re
import glob
import pandas as pd

# paths
DATA_DIR     = "data"
FASTA        = os.path.join(DATA_DIR, "stentor_CaMKII_homolog_candidates_full_length.fasta")
INTERPRO_DIR = os.path.join(DATA_DIR, "interpro_files")
RESULTS_DIR  = "results"
OUT_TABLE    = os.path.join(RESULTS_DIR, "reannotation_summary_AUTO.tsv")

os.makedirs(RESULTS_DIR, exist_ok=True)
DASH = "\u2013"  # en dash for coordinate ranges

# load FASTAs
def load_fasta(path):
    recs = {}
    pid, ann, seq = None, "", []
    def flush():
        if pid is not None:
            recs[pid] = {"annotation": ann.strip(), "seq": "".join(seq)}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                parts = re.split(r"\s+", header, maxsplit=1)
                pid = parts[0]
                ann = parts[1] if len(parts) > 1 else ""
                seq = []
            else:
                seq.append(line.strip())
        flush()
    return recs

def x_report(seq):
    # locate ambiguous residues (X); map any suspected duplications around them
    positions = [i + 1 for i, c in enumerate(seq) if c.upper() == "X"]  
    dup = False
    for i0 in [p - 1 for p in positions]:
        pre = seq[max(0, i0 - 10):i0]
        post = seq[i0 + 1:i0 + 40]
        if len(pre) >= 8 and pre[-8:] in post:
            dup = True
    return positions, dup

# load interpro TSVs
CANON = [
    "signature_accession", "signature_name", "short_name", "source_db",
    "feature_type", "integrated_into", "integrated_signatures", "go_terms",
    "protein_uniprot_acc", "protein_length", "matches",
]
HEADER_TOKENS = {"accession", "source database", "matches", "protein length"}

def normalise_columns(df):
    lower = [str(c).strip().lower() for c in df.columns]
    looks_like_header = sum(tok in lower for tok in HEADER_TOKENS) >= 2
    if looks_like_header and len(df.columns) == len(CANON):
        df = df.copy(); df.columns = CANON
        inferred = "positional map onto standard InterPro-match columns (header present)"
    elif looks_like_header:
        name_map = {}
        for c, lc in zip(df.columns, lower):
            if lc.startswith("accession"):     name_map[c] = "signature_accession"
            elif lc == "name":                 name_map[c] = "signature_name"
            elif "short" in lc:                name_map[c] = "short_name"
            elif "source" in lc:               name_map[c] = "source_db"
            elif lc == "type":                 name_map[c] = "feature_type"
            elif "integrated into" in lc:      name_map[c] = "integrated_into"
            elif "integrated signat" in lc:    name_map[c] = "integrated_signatures"
            elif "go" in lc:                   name_map[c] = "go_terms"
            elif "protein accession" in lc:    name_map[c] = "protein_uniprot_acc"
            elif "protein length" in lc:       name_map[c] = "protein_length"
            elif "match" in lc:                name_map[c] = "matches"
        df = df.rename(columns=name_map)
        inferred = "fuzzy header map (non-standard column count)"
    else:
        df = df.copy(); df.columns = CANON[:len(df.columns)]
        inferred = "no header detected; positional map onto standard InterPro-match columns"
    for c in CANON:
        if c not in df.columns:
            df[c] = ""
    return df[CANON], inferred

def load_interpro(path):
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, header=0)
    df, inferred = normalise_columns(raw)
    df = df[~(df["signature_accession"].str.strip() == "")].reset_index(drop=True)
    return df, inferred

# coordinate and classification helpers
def parse_segments(matches):
    segs = []
    for chunk in str(matches).split(","):
        chunk = chunk.strip()
        m = re.match(r"^(\d+)\.\.(\d+)$", chunk)
        if m:
            segs.append((int(m.group(1)), int(m.group(2))))
        elif re.match(r"^\d+$", chunk):
            segs.append((int(chunk), int(chunk)))
    return segs

KIN_RE  = re.compile(r"kinase|stkc|s_tkc|pkinase|phosphotransferase|phosphorylase kinase|kinase-like", re.I)
CDPK_RE = re.compile(r"calcium-dependent .*protein kinase|cdpk", re.I)
HUB_RE  = re.compile(r"association domain|hub domain|camkii_ad|camk2.*assoc|calcium/calmodulin-dependent protein kinase ii assoc", re.I)
EF_RE   = re.compile(r"ef-?hand|efh|ef_hand|calmodulin|calcium.binding", re.I)

def is_hub_row(row):
    # CaMKII association/hub-domain signature (Pfam PF08332 / IPR013543)
    blob = f"{row['signature_accession']} {row['signature_name']} {row['short_name']}"
    return row["signature_accession"] in {"PF08332", "IPR013543"} or bool(HUB_RE.search(blob))

def span(rows):
    # outer (min start, max end) across a list of (start,end)
    if not rows:
        return None
    return (min(s for s, _ in rows), max(e for _, e in rows))

def fmt(a, b):
    return f"{a}{DASH}{b}"

# parse candidates
fasta = load_fasta(FASTA)
tsv_paths = sorted(glob.glob(os.path.join(INTERPRO_DIR, "SteCoe_*_interpro.tsv")))
if not tsv_paths:
    raise SystemExit(f"No SteCoe_*_interpro.tsv files found in {INTERPRO_DIR!r}")

# kinase-domain source priority (domain-level Pfam/InterPro before broad superfamily/CDD)
KIN_PRIORITY = [("pfam", "PF00069"), ("interpro", "IPR000719"),
                ("smart", "SM00220"), ("profile", "PS50011"), ("cdd", "cd05117")]

summary_rows = []  
flags        = []   
inferred_notes = {}

for path in tsv_paths:
    pid = os.path.basename(path).replace("_interpro.tsv", "")
    df, inferred = load_interpro(path)
    inferred_notes[pid] = inferred

    fa = fasta.get(pid, {})
    seq = fa.get("seq", "")
    fasta_len = len(seq)
    ip_len_vals = sorted({int(v) for v in df["protein_length"] if str(v).strip().isdigit()})
    ip_len = ip_len_vals[0] if ip_len_vals else None
    length_used = fasta_len if fasta_len else ip_len
    length_mismatch = (ip_len is not None and fasta_len and ip_len != fasta_len)
    xpos, xdup = x_report(seq)

    # collect coordinate segments per (source, accession)
    per_source, hub_hits = {}, []
    for _, r in df.iterrows():
        segs = parse_segments(r["matches"])
        per_source.setdefault((r["source_db"].lower(), r["signature_accession"]), []).extend(segs)
        if is_hub_row(r):
            s, e = span(segs) if segs else (None, None)
            hub_hits.append((r["source_db"], r["signature_accession"], s, e))

    def segs_of(src, acc):
        return per_source.get((src, acc), [])

    # kinase evidence
    kin_source, kin_span = None, None
    for src, acc in KIN_PRIORITY:
        if segs_of(src, acc):
            kin_source, kin_span = src, span(segs_of(src, acc)); break
    cdd_stkc  = span(segs_of("cdd", "cd05117"))
    serthr_as = segs_of("prosite", "PS00108")
    tyr_as    = segs_of("prosite", "PS00109")
    smart_kin = segs_of("smart", "SM00220")
    has_serthr_as, has_tyr_as = bool(serthr_as), bool(tyr_as)
    pseudokinase = not has_serthr_as

    kin_bits = []
    if kin_source:
        pretty = {"pfam": "Pfam Pkinase (PF00069)", "interpro": "InterPro Prot_kinase_dom (IPR000719)",
                  "smart": "SMART S_TKc (SM00220)", "profile": "PROSITE profile (PS50011)",
                  "cdd": "CDD STKc_CAMK (cd05117)"}[kin_source]
        kin_bits.append(f"{pretty}: {fmt(*kin_span)}")
    if cdd_stkc and kin_source != "cdd":
        kin_bits.append(f"CDD STKc_CAMK (cd05117): {fmt(*cdd_stkc)} [CAMK kinase group]")
    if has_serthr_as:
        kin_bits.append(f"Ser/Thr active site (PS00108): {fmt(*span(serthr_as))}")
    elif has_tyr_as:
        kin_bits.append(f"NO Ser/Thr active-site signature; a Tyr active-site pattern "
                        f"(PS00109) is present at {fmt(*span(tyr_as))} \u2014 divergent/degenerate "
                        f"catalytic loop, not evidence of Tyr-kinase identity")
    else:
        kin_bits.append("no Ser/Thr or Tyr active-site signature detected (degenerate catalytic loop)")
    if not smart_kin:
        kin_bits.append("SMART S_TKc not detected")
    kin_evidence = "; ".join(kin_bits)

    # EF-hand / CaM-like evidence
    pfam_pair_rows, pfam_single_rows = [], []
    for (src, acc), segs in per_source.items():
        if src != "pfam":
            continue
        nm = df[(df["source_db"].str.lower() == "pfam") &
                (df["signature_accession"] == acc)]["signature_name"].iloc[0]
        if EF_RE.search(nm):
            (pfam_pair_rows if re.search("pair", nm, re.I) else pfam_single_rows).append((acc, nm, segs))
    smart_ef   = segs_of("smart", "SM00054")
    prosite_ef = segs_of("prosite", "PS00018")
    profile_ef = segs_of("profile", "PS50222")
    n_pfam_pair, n_smart_ef  = sum(len(s) for _, _, s in pfam_pair_rows), len(smart_ef)
    n_prosite_ef, n_profile_ef = len(prosite_ef), len(profile_ef)

    ef_bits = []
    if pfam_pair_rows:
        pr = [f"{acc}: " + " & ".join(fmt(a, b) for a, b in segs) for acc, _, segs in pfam_pair_rows]
        ef_bits.append(f"Pfam EF-hand-pair [{'; '.join(pr)}] = {n_pfam_pair} pair-domain(s)")
    else:
        ef_bits.append("Pfam EF-hand-pair: not detected")
    if pfam_single_rows:
        ps = [f"{acc} ({nm}): " + " / ".join(fmt(a, b) for a, b in segs) for acc, nm, segs in pfam_single_rows]
        ef_bits.append("Pfam single EF-hand [" + "; ".join(ps) + "]")
    if smart_ef:
        ef_bits.append("SMART EFh (SM00054): " + " / ".join(fmt(a, b) for a, b in smart_ef) +
                       f" = {n_smart_ef} motif(s)")
    else:
        ef_bits.append("SMART EFh (SM00054): no individual motifs detected")
    ef_evidence = "; ".join(ef_bits)

    # CDPK family
    has_cdpk = df["signature_accession"].isin(["IPR050205", "PTHR24349"]).any()
    fam_bits = []
    if segs_of("panther", "PTHR24349"):
        fam_bits.append(f"PANTHER PTHR24349 CDPK_Ser/Thr_kinases: {fmt(*span(segs_of('panther','PTHR24349')))}")
    if (df["signature_accession"] == "IPR050205").any():
        fam_bits.append("InterPro IPR050205 (Calcium-dependent Ser/Thr protein kinases)")
    family_str = "; ".join(fam_bits) if fam_bits else "no CDPK family assignment detected"

    # CaMKII association/hub domain
    has_hub = len(hub_hits) > 0
    if has_hub:
        s_db, acc, s, e = hub_hits[0]
        hub_str = f"{s_db} {acc} {fmt(s, e)}"
    else:
        hub_str = "Not detected in provided annotations"
    has_ef = bool(pfam_pair_rows or pfam_single_rows or smart_ef or profile_ef)

    suspected_truncation = bool(xpos and length_used and min(length_used - p for p in xpos) <= 5) \
                           and (n_pfam_pair < 2 and n_smart_ef == 0)
    if not has_cdpk or kin_source is None:
        revised = "manual review required"
    elif suspected_truncation:
        revised = "probable truncated CDPK-like calcium-dependent protein kinase"
    elif pseudokinase:
        revised = "CDPK-like calcium-dependent protein kinase (predicted pseudokinase; atypical catalytic loop)"
    else:
        revised = "CDPK-like calcium-dependent protein kinase"

    # collect
    # NOTE: `new_gene_id_if_mapped` and a `Notes` column are added during manual
    summary_rows.append({
        "Protein_ID": pid,
        "Length_aa": length_used,
        "StentorDB_annotation": fa.get("annotation", ""),
        "Kinase_domain_evidence": kin_evidence,
        "EF_hand_CaM_like_evidence": ef_evidence,
        "InterPro_PANTHER_family": family_str,
        "CaMKII_association_hub_domain": hub_str,
        "Revised_annotation": revised,
    })
    flags.append({"has_kinase": kin_span is not None, "has_ef": has_ef,
                  "has_cdpk": bool(has_cdpk), "has_hub": has_hub})

# write table
t1 = pd.DataFrame(summary_rows).sort_values("Protein_ID").reset_index(drop=True)
t1.to_csv(OUT_TABLE, sep="\t", index=False)

# summary output to terminal
flags_df = pd.DataFrame(flags)
n = len(t1)
print("CDPK-vs-CaMKII reannotation \u2014 summary table")
print(f"Proteins processed:                              {n}")
print(f"  with kinase-domain evidence:                   {int(flags_df['has_kinase'].sum())}")
print(f"  with C-terminal EF-hand/CaM-like evidence:     {int(flags_df['has_ef'].sum())}")
print(f"  with CDPK PANTHER/InterPro family evidence:    {int(flags_df['has_cdpk'].sum())}")
print(f"  with CaMKII association/hub-domain evidence:   {int(flags_df['has_hub'].sum())}")
print("Revised-annotation tally:")
for k, v in t1["Revised_annotation"].value_counts().items():
    print(f"   {v:>2}  {k}")
print(f"Wrote: {OUT_TABLE}")
