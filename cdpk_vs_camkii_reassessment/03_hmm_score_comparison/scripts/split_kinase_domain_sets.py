import sys
from pathlib import Path

if len(sys.argv) != 7:
    sys.exit(
        "Usage: python split_kinase_domain_sets.py"
        "<all_kinase_domains.fasta> "
        "<cdpk_accessions.txt> "
        "<camkii_accessions.txt> "
        "<cdpk_out.fasta> "
        "<camkii_out.fasta> "
        "<stentor_out.fasta>"
    )

fasta_file = Path(sys.argv[1])
cdpk_list = Path(sys.argv[2])
camkii_list = Path(sys.argv[3])
cdpk_out = Path(sys.argv[4])
camkii_out = Path(sys.argv[5])
stentor_out = Path(sys.argv[6])

cdpk_accs = {x.strip() for x in cdpk_list.read_text().splitlines() if x.strip()}
camkii_accs = {x.strip() for x in camkii_list.read_text().splitlines() if x.strip()}

def read_fasta(path):
    records = []
    header = None
    seq_chunks = []

    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

    if header is not None:
        records.append((header, "".join(seq_chunks)))

    return records

def write_fasta(records, path):
    with open(path, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")

def extract_uniprot_accession(header):
    fields = header.split("|")
    if len(fields) >= 3 and fields[0] in {"sp", "tr"}:
        return fields[1]
    return None


records = read_fasta(fasta_file)

cdpk_records = []
camkii_records = []
stentor_records = []
unassigned_records = []

for header, seq in records:
    first_token = header.split()[0]
    accession = extract_uniprot_accession(first_token)

    if first_token.startswith("SteCoe_"):
        stentor_records.append((header, seq))
    elif accession in cdpk_accs:
        cdpk_records.append((header, seq))
    elif accession in camkii_accs:
        camkii_records.append((header, seq))
    else:
        unassigned_records.append((header, seq))

write_fasta(cdpk_records, cdpk_out)
write_fasta(camkii_records, camkii_out)
write_fasta(stentor_records, stentor_out)

print(f"CDPK/CDPK-like training sequences: {len(cdpk_records)}")
print(f"CaMKII training sequences: {len(camkii_records)}")
print(f"Stentor test sequences: {len(stentor_records)}")
print(f"Unassigned kinase-domain sequences: {len(unassigned_records)}")

if len(stentor_records) == 0:
    sys.exit("ERROR: no Stentor sequences found. Make sure headers begin with 'SteCoe_'.")

missing_cdpk = sorted(cdpk_accs - {extract_uniprot_accession(h.split()[0]) for h, s in cdpk_records})
missing_camkii = sorted(camkii_accs - {extract_uniprot_accession(h.split()[0]) for h, s in camkii_records})

if missing_cdpk:
    print("WARNING: CDPK accessions not found in kinase-domain FASTA:")
    for acc in missing_cdpk:
        print(f"  {acc}")

if missing_camkii:
    print("WARNING: CaMKII accessions not found in kinase-domain FASTA:")
    for acc in missing_camkii:
        print(f"  {acc}")