set -euo pipefail

mkdir -p results logs generated

# count input sequences 
grep -c "^>" data/stentor_CaMKII_homolog_candidates_full_length.fasta \
    > logs/stentor_input_sequence_count.txt

grep -c "^>" data/camk_reference_uniprot_full_length.fasta \
    > logs/reference_input_sequence_count.txt

# combine full-length reference and Stentor FASTAs
cat data/camk_reference_uniprot_full_length.fasta \
    data/stentor_CaMKII_homolog_candidates_full_length.fasta \
    > generated/all_full_length_proteins.fasta

# search all full-length proteins with PF00069"
hmmsearch \
    --cut_ga \
    --domtblout results/PF00069_vs_all_full_length.domtblout \
    data/PF00069.hmm \
    generated/all_full_length_proteins.fasta \
    > logs/PF00069_vs_all_full_length.hmmsearch.txt

# extract best PF00069 kinase-domain slice per protein
python scripts/extract_best_hmmer_domain.py \
    results/PF00069_vs_all_full_length.domtblout \
    generated/all_full_length_proteins.fasta \
    results/all_PF00069_kinase_domains.fasta \
    > logs/extract_best_hmmer_domain.log \
    2> logs/extract_best_hmmer_domain.err

# split kinase-domain sequences into CDPK training, CaMKII training, and Stentor test sets
python scripts/split_kinase_domain_sets.py \
    results/all_PF00069_kinase_domains.fasta \
    data/cdpk_training_accessions.txt \
    data/camkii_training_accessions.txt \
    results/CDPK_training_KD.fasta \
    results/CaMKII_training_KD.fasta \
    results/stentor_CaMKII_homolog_candidates_test_KD.fasta \
    > logs/split_kinase_domain_sets.log

# check stentor kinase-domain test set count
STENTOR_COUNT=$(grep -c "^>" results/stentor_CaMKII_homolog_candidates_test_KD.fasta)
if [[ "$STENTOR_COUNT" -ne 14 ]]; then
    echo "ERROR: Expected 14 Stentor kinase-domain sequences, found $STENTOR_COUNT" >&2
    exit 1
fi

# align CDPK/CDPK-like training kinase domains
mafft --localpair --maxiterate 1000 \
    results/CDPK_training_KD.fasta \
    > results/CDPK_training_KD.mafft.fasta \
    2> logs/CDPK_training_KD.mafft.log

# align CaMKII training kinase domains
mafft --localpair --maxiterate 1000 \
    results/CaMKII_training_KD.fasta \
    > results/CaMKII_training_KD.mafft.fasta \
    2> logs/CaMKII_training_KD.mafft.log

# build custom CDPK/CDPK-like kinase-domain HMM
hmmbuild \
    results/CDPK_KD.hmm \
    results/CDPK_training_KD.mafft.fasta \
    > logs/CDPK_KD.hmmbuild.log

# build custom CaMKII kinase-domain HMM
hmmbuild \
    results/CaMKII_KD.hmm \
    results/CaMKII_training_KD.mafft.fasta \
    > logs/CaMKII_KD.hmmbuild.log

# score stentor kinase domains against CDPK/CDPK-like model
hmmsearch \
    --max \
    -E 1000 \
    --domE 1000 \
    --tblout results/stentor_CaMKII_homolog_candidates_vs_CDPK_KD.tbl \
    results/CDPK_KD.hmm \
    results/stentor_CaMKII_homolog_candidates_test_KD.fasta \
    > logs/stentor_CaMKII_homolog_candidates_vs_CDPK_KD.hmmsearch.txt

# score stentor kinase domains against CaMKII model
hmmsearch \
    --max \
    -E 1000 \
    --domE 1000 \
    --tblout results/stentor_CaMKII_homolog_candidates_vs_CaMKII_KD.tbl \
    results/CaMKII_KD.hmm \
    results/stentor_CaMKII_homolog_candidates_test_KD.fasta \
    > logs/stentor_CaMKII_homolog_candidates_vs_CaMKII_KD.hmmsearch.txt

# make final CDPK-vs-CaMKII score table
python scripts/make_cdpk_vs_camkii_score_table.py \
    results/stentor_CaMKII_homolog_candidates_vs_CDPK_KD.tbl \
    results/stentor_CaMKII_homolog_candidates_vs_CaMKII_KD.tbl \
    results/stentor_CaMKII_homolog_candidates_CDPK_vs_CaMKII_KD_scores.tsv

echo "Done."