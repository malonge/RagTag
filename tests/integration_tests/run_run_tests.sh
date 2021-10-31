#!/usr/bin/env bash

bash run_tests.sh ~/Projects/Reference_Genomes/Arabidopsis/TAIR10_chr_all.fasta \
    ~/Projects/ragtag_workspace/ara_test_env/Arabidopsis.MaSuRCA.3.4.0.fasta \
    ~/Projects/ragtag_workspace/ara_test_env/Arabidopsis.MaSuRCA.3.4.0.genes.gff \
    ~/Projects/ragtag_workspace/ecoli_tests/GCF_000008865.2_ASM886v2_genomic.fna \
    ~/Projects/ragtag_workspace/ecoli_tests/ecoli.contigs.fasta \
    ~/Projects/ragtag_workspace/ecoli_tests/reads.fofn \
    ~/Projects/ragtag_workspace/agp_reconcile_tests/ara_panel/ragtag_output_Ler/ragtag.scaffolds.agp \
    ~/Projects/ragtag_workspace/agp_reconcile_tests/ara_panel/ragtag_output_Ler/ragtag.scaffolds.cp.agp \
    ~/Projects/ragtag_workspace/agp_reconcile_tests/ara_panel/ragtag_output_Col_0/ragtag.scaffolds.agp
