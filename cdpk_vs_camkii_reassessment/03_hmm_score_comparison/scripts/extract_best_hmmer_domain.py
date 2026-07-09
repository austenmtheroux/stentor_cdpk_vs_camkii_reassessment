import sys
from pathlib import Path

if len(sys.argv) != 4:
    sys.exit(
        "Usage: python extract_best_hmmer_domain.py"
        "<hmmer.domtblout> <input_full_length.fasta> <output_domain_slices.fasta>"
    )

domtblout = Path(sys.argv[1])
fasta_in = Path(sys.argv[2])
fasta_out = Path(sys.argv[3])

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

records = read_fasta(fasta_in)
seqs = {header.split()[0]: seq for header, seq in records}
header_by_id = {header.split()[0]: header for header, seq in records}
order = [header.split()[0] for header, seq in records]

best = {}

with open(domtblout) as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.split()

        target = fields[0]
        domain_score = float(fields[13])
        hmm_from = int(fields[15])
        hmm_to = int(fields[16])
        ali_from = int(fields[17])
        ali_to = int(fields[18])
        env_from = int(fields[19])
        env_to = int(fields[20])
        mean_post_prob = float(fields[21])

        if target not in best or domain_score > best[target]["domain_score"]:
            best[target] = {
                "domain_score": domain_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "env_from": env_from,
                "env_to": env_to,
                "mean_post_prob": mean_post_prob,
            }

with open(fasta_out, "w") as out:
    for protein_id in order:
        if protein_id not in best:
            print(f"WARNING: no PF00069 domain hit for {protein_id}", file=sys.stderr)
            continue

        hit = best[protein_id]
        seq = seqs[protein_id]

        start = hit["ali_from"]
        end = hit["ali_to"]
        subseq = seq[start - 1:end]

        out_header = (
            f"{protein_id}"
            f"|PF00069_ali_{start}-{end}"
            f"|hmm_{hit['hmm_from']}-{hit['hmm_to']}"
            f"|score_{hit['domain_score']}"
        )

        out.write(f">{out_header}\n")
        for i in range(0, len(subseq), 60):
            out.write(subseq[i:i + 60] + "\n")