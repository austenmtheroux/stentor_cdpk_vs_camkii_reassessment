set -euo pipefail

# Search the CaMKII association-domain HMM PF08332 against the 2017 Slabodnick et al. proteome (S_coeruleus_Nov2017_proteins.fasta) 
# and present-study proteome (SteCoeF_v1.0.3.final.CDS.pep.fasta)

mkdir -p results

hmmsearch --max \
  -E 1e200 \
  --domE 1e200 \
  --tblout results/PF08332_vs_S_coeruleus_Nov2017_proteins.tbl \
  --domtblout results/PF08332_vs_S_coeruleus_Nov2017_proteins.domtbl \
  data/PF08332.hmm \
  data/S_coeruleus_Nov2017_proteins.fasta \
  > results/PF08332_vs_S_coeruleus_Nov2017_proteins.out

hmmsearch --max \
  -E 1e200 \
  --domE 1e200 \
  --tblout results/PF08332_vs_SteCoeF_v1.0.3.final.CDS.pep.tbl \
  --domtblout results/PF08332_vs_SteCoeF_v1.0.3.final.CDS.pep.domtbl \
  data/PF08332.hmm \
  data/SteCoeF_v1.0.3.final.CDS.pep.fasta \
  > results/PF08332_vs_SteCoeF_v1.0.3.final.CDS.pep.out

echo "Done"