#!/usr/bin/env bash

# position args:
## 1. tomato reference
## 2. tomato query
## 3. tomato gff
## 4. ara reference
## 5. ara query
## 6. ara gff
## 7. E. coli reference
## 8. E. coli query
## 9. E.coli reads fofn

# Assumes fasta suffixes are ".fasta"
# Assumes gff suffixes are ".gff"

Usage() {
    echo "Usage: $0 tomato.ref.fa tomato.query.fa tomato.genes.gff ara.ref.fa ara.query.fa ara.genes.gff"
}

mecho() {
    NAME=`basename $0`
    echo "$NAME:" $1
}

if [ $# -lt 6 ] ; then
    Usage
    exit 1
fi

T_REF=$1
T_QUERY=$2
T_QUERY_PREF=`basename $T_QUERY .fasta`
T_GFF=$3
T_GFF_PREF=`basename $T_GFF .gff`

A_REF=$4
A_QUERY=$5
A_QUERY_PREF=`basename $A_QUERY .fasta`
A_GFF=$6
A_GFF_PREF=`basename $A_GFF .gff`

E_REF=$7
E_QUERY=$8
E_QUERY_PREF=`basename $E_QUERY .fasta`
E_VAL_FOFN=$9

M_ASM=$A_QUERY
M_AGP_1=${10}
M_AGP_2=${11}
M_AGP_3=${12}


# Run RagTag
# Settings: default
# Data: Tomato
echo ""
echo "*************************************************************"
echo "*** Running RagTag on tomato data with default parameters ***"
echo "*************************************************************"
echo ""
bash scripts/run_default.sh $T_REF \
    $T_QUERY \
    $T_GFF \
    ragtag_output_tomato_default

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh $T_QUERY \
    ragtag_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragtag_output_tomato_default/ragtag.correction.agp

bash scripts/validate_agp.sh ragtag_output_tomato_default/ragtag.scaffolds.fasta \
    ragtag_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragtag_output_tomato_default/ragtag.scaffolds.agp

# Validate the gff files
echo ""
mecho "Validating GFF files:"
echo ""

bash scripts/validate_gff.sh $T_QUERY \
    $T_GFF \
    ragtag_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragtag_output_tomato_default/$T_GFF_PREF.corr.gff

bash scripts/validate_gff.sh ragtag_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragtag_output_tomato_default/$T_GFF_PREF.corr.gff \
    ragtag_output_tomato_default/ragtag.scaffolds.fasta \
    ragtag_output_tomato_default/$T_GFF_PREF.scaf.gff


# Validate the results
echo ""
mecho "Validating results:"
echo ""

bash scripts/validate_results.sh ragtag_output_tomato_default/ragtag.correction.agp \
    ~/Projects/ragtag_workspace/static_results/tomato/ragtag.correction.agp

bash scripts/validate_results.sh ragtag_output_tomato_default/ragtag.scaffolds.agp \
    ~/Projects/ragtag_workspace/static_results/tomato/ragtag.scaffolds.agp


# Run RagTag
# Settings: with nucmer
# Data: Arabidopsis
echo ""
echo "**************************************************************"
echo "***     Running RagTag on Arabidopsis data with Nucmer     ***"
echo "**************************************************************"
echo ""

bash scripts/run_nucmer.sh $A_REF \
    $A_QUERY \
    $A_GFF \
    ragtag_output_Ara_nucmer

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh $A_QUERY \
    ragtag_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragtag_output_Ara_nucmer/ragtag.correction.agp

bash scripts/validate_agp.sh ragtag_output_Ara_nucmer/ragtag.scaffolds.fasta \
    ragtag_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragtag_output_Ara_nucmer/ragtag.scaffolds.agp

# Validate the gff files
echo ""
mecho "Validating GFF files:"
echo ""

bash scripts/validate_gff.sh $A_QUERY \
    $A_GFF \
    ragtag_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragtag_output_Ara_nucmer/$A_GFF_PREF.corr.gff

bash scripts/validate_gff.sh ragtag_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragtag_output_Ara_nucmer/$A_GFF_PREF.corr.gff \
    ragtag_output_Ara_nucmer/ragtag.scaffolds.fasta \
    ragtag_output_Ara_nucmer/$A_GFF_PREF.scaf.gff

# Validate the unique anchor filtering
echo ""
mecho "Validating alignment filtering:"
echo ""

bash scripts/validate_uaf.sh ragtag_output_Ara_nucmer/c_query_against_ref.delta \
    ragtag_output_Ara_nucmer/ragtag.correction.debug.filtered.paf \
    1000

bash scripts/validate_uaf.sh ragtag_output_Ara_nucmer/query_against_ref.delta \
    ragtag_output_Ara_nucmer/ragtag.scaffolds.debug.filtered.paf \
    1000

# Validate the results
echo ""
mecho "Validating results:"
echo ""

bash scripts/validate_results.sh ragtag_output_Ara_nucmer/ragtag.correction.agp \
    ~/Projects/ragtag_workspace/static_results/ara/ragtag.correction.agp

bash scripts/validate_results.sh ragtag_output_Ara_nucmer/ragtag.scaffolds.agp \
    ~/Projects/ragtag_workspace/static_results/ara/ragtag.scaffolds.agp


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
    ragtag_output_ecoli_val/$E_QUERY_PREF.corrected.fasta \
    ragtag_output_ecoli_val/ragtag.correction.agp


# Run ragtag merge
# Settings: default
# Data: Arabidopsis, same AGP twice
echo ""
echo "*******************************************************************************"
echo "***     Running ragtag merge on identical AGPs with default parameters     ***"
echo "*******************************************************************************"
echo ""

bash scripts/run_merge_default.sh ragtag_output_merge_same \
    $M_ASM
    $M_AGP_1
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
    $M_ASM
    $M_AGP_1
    $M_AGP_3

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh ragtag_output_merge_diff/ragtag.merge.fasta \
    $M_ASM \
    ragtag_output_merge_same/ragtag.merge.agp