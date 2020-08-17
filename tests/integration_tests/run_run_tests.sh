#!/usr/bin/env bash

bash run_tests.sh ~/Projects/Reference_Genomes/Tomato/SL4.0.genome.fasta \
    ~/Projects/ragtag_workspace/tomato_test_env/Brandywine.assembly.polished.fasta \
    ~/Projects/ragtag_workspace/tomato_test_env/brandywine_ctg_genes.gff \
    ~/Projects/Reference_Genomes/Arabidopsis/GCF_000001735.3_TAIR10_genomic.fna \
    ~/Projects/ragtag_workspace/ara_test_env/Arabidopsis.MaSuRCA.3.4.0.fasta \
    ~/Projects/ragtag_workspace/ara_test_env/Arabidopsis.MaSuRCA.3.4.0.genes.gff \
    ~/Projects/ragtag_workspace/ecoli_tests/GCF_000008865.2_ASM886v2_genomic.fna \
    ~/Projects/ragtag_workspace/ecoli_tests/ecoli.contigs.fasta \
    ~/Projects/ragtag_workspace/ecoli_tests/reads.fofn
