#!/usr/bin/env bash

# position args:
## 1. A. thaliana reference
## 2. A. thaliana query
## 3. A. thaliana gff
## 4. E. coli reference
## 5. E. coli query
## 6. E.coli reads fofn
## 7. A. thaliana solution against Ler
## 8. Copy of A. thaliana solution against Ler
## 9. A. thaliana solution against Col-0

# Assumes fasta suffixes are ".fasta"
# Assumes gff suffixes are ".gff"

mecho() {
    NAME=`basename $0`
    echo "$NAME:" $1
}

if [ $# -lt 6 ] ; then
    Usage
    exit 1
fi

A_REF=$1
A_QUERY=$2
A_QUERY_PREF=`basename $A_QUERY .fasta`
A_GFF=$3
A_GFF_PREF=`basename $A_GFF .gff`

E_REF=$4
E_QUERY=$5
E_QUERY_PREF=`basename $E_QUERY .fasta`
E_VAL_FOFN=$6

M_ASM=$A_QUERY
M_AGP_1=${7}
M_AGP_2=${8}
M_AGP_3=${9}


# Run RagTag
# Settings: default
# Data: A. thaliana
echo ""
echo "******************************************************************"
echo "*** Running RagTag on A. thaliana data with default parameters ***"
echo "******************************************************************"
echo ""

OUTDIR=ragtag_output_Ara_default
bash scripts/run_default.sh $A_REF \
    $A_QUERY \
    $A_GFF \
    $OUTDIR

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh $A_QUERY \
    ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/ragtag.correct.agp

bash scripts/validate_agp.sh ${OUTDIR}/ragtag.scaffold.fasta \
    ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/ragtag.scaffold.agp

# Validate the gff files
echo ""
mecho "Validating GFF files:"
echo ""

bash scripts/validate_gff.sh $A_QUERY \
    $A_GFF \
    ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/$A_GFF_PREF.corr.gff

bash scripts/validate_gff.sh ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/$A_GFF_PREF.corr.gff \
    ${OUTDIR}/ragtag.scaffold.fasta \
    ${OUTDIR}/$A_GFF_PREF.scaf.gff


# Validate the results
echo ""
mecho "Validating results:"
echo ""

# TODO replace with new results
bash scripts/validate_results.sh ${OUTDIR}/ragtag.correct.agp \
    ~/Projects/ragtag_workspace/static_results/tomato/ragtag.correct.agp

bash scripts/validate_results.sh ragtag_output_tomato_default/ragtag.scaffold.agp \
    ~/Projects/ragtag_workspace/static_results/tomato/ragtag.scaffold.agp


# Run RagTag
# Settings: with nucmer
# Data: Arabidopsis
echo ""
echo "**************************************************************"
echo "***     Running RagTag on Arabidopsis data with Nucmer     ***"
echo "**************************************************************"
echo ""
OUTDIR=ragtag_output_Ara_nucmer

bash scripts/run_nucmer.sh $A_REF \
    $A_QUERY \
    $A_GFF \
    $OUTDIR

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh $A_QUERY \
    ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/ragtag.correct.agp

bash scripts/validate_agp.sh ${OUTDIR}/ragtag.scaffold.fasta \
    ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/ragtag.scaffold.agp

# Validate the gff files
echo ""
mecho "Validating GFF files:"
echo ""

bash scripts/validate_gff.sh $A_QUERY \
    $A_GFF \
    ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/$A_GFF_PREF.corr.gff

bash scripts/validate_gff.sh ${OUTDIR}/ragtag.correct.fasta \
    ${OUTDIR}/$A_GFF_PREF.corr.gff \
    ${OUTDIR}/ragtag.scaffold.fasta \
    ${OUTDIR}/$A_GFF_PREF.scaf.gff

# Validate the unique anchor filtering
echo ""
mecho "Validating alignment filtering:"
echo ""

bash scripts/validate_uaf.sh ${OUTDIR}/ragtag.correct.asm.delta \
    ${OUTDIR}/ragtag.correct.debug.filtered.paf \
    1000

bash scripts/validate_uaf.sh ${OUTDIR}/ragtag.scaffold.asm.delta \
    ${OUTDIR}/ragtag.scaffold.debug.filtered.paf \
    1000

# Validate the results
echo ""
mecho "Validating results:"
echo ""

# TODO update results
bash scripts/validate_results.sh ragtag_output_Ara_nucmer/ragtag.correct.agp \
    ~/Projects/ragtag_workspace/static_results/ara/ragtag.correct.agp

bash scripts/validate_results.sh ragtag_output_Ara_nucmer/ragtag.scaffold.agp \
    ~/Projects/ragtag_workspace/static_results/ara/ragtag.scaffold.agp


# Run RagTag
# Settings: with error-corrected long-read validation
# Data: E. coli
echo ""
echo "**************************************************************"
echo "***     Running RagTag on E. coli data with validation     ***"
echo "**************************************************************"
echo ""

bash scripts/run_val.sh $E_REF \
    $E_QUERY \
    $E_VAL_FOFN \
    ragtag_output_ecoli_val

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh $E_QUERY \
    ragtag_output_ecoli_val/ragtag.correct.fasta \
    ragtag_output_ecoli_val/ragtag.correct.agp


# Run ragtag merge
# Settings: default
# Data: Arabidopsis, same AGP twice
echo ""
echo "*******************************************************************************"
echo "***     Running ragtag merge on identical AGPs with default parameters     ***"
echo "*******************************************************************************"
echo ""

bash scripts/run_merge_default.sh ragtag_output_merge_same \
    $M_ASM \
    $M_AGP_1 \
    $M_AGP_2

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh ragtag_output_merge_same/ragtag.merge.fasta \
    $M_ASM \
    ragtag_output_merge_same/ragtag.merge.agp

# Since the input files define the same objects, check that the input and output AGPs define the same objects
bash scripts/validate_agp_eq.sh $M_ASM \
    $M_AGP_1 \
    $M_AGP_2

bash scripts/validate_agp_eq.sh $M_ASM \
    $M_AGP_1 \
    ragtag_output_merge_same/ragtag.merge.agp

bash scripts/validate_agp_eq.sh $M_ASM \
    $M_AGP_2 \
    ragtag_output_merge_same/ragtag.merge.agp


# Run ragtag merge
# Settings: default
# Data: Arabidopsis, two different references
echo ""
echo "*****************************************************************************************************"
echo "***     Running ragtag merge on combing scaffolds from two references with default parameters     ***"
echo "*****************************************************************************************************"
echo ""

bash scripts/run_merge_default.sh ragtag_output_merge_diff \
    $M_ASM \
    $M_AGP_1 \
    $M_AGP_3

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh ragtag_output_merge_diff/ragtag.merge.fasta \
    $M_ASM \
    ragtag_output_merge_diff/ragtag.merge.agp