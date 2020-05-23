#!/usr/bin/env bash

# position args:
## 1. tomato reference
## 2. tomato query
## 3. tomato gff
## 4. ara reference
## 5. ara query
## 6. ara gff

# Assumes fasta suffixes are ".fasta"
# Assumes gff suffixes are ".gff"

Usage() {
    echo "Usage: $0 tomato.ref.fa tomato.query.fa tomato.genes.gff ara.ref.fa ara.query.fa ara.genes.gff"
}

mecho() {
    echo "Master:" $1
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


# Run RaGOO2
# Settings: default
# Data: Tomato
echo ""
echo "*************************************************************"
echo "*** Running RaGOO2 on tomato data with default parameters ***"
echo "*************************************************************"
echo ""
bash scripts/run_default.sh $T_REF \
    $T_QUERY \
    $T_GFF \
    ragoo2_output_tomato_default

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh $T_QUERY \
    ragoo2_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragoo2_output_tomato_default/ragoo2.correction.agp

bash scripts/validate_agp.sh ragoo2_output_tomato_default/ragoo2.scaffolds.fasta \
    ragoo2_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragoo2_output_tomato_default/ragoo2.scaffolds.agp

# Validate the gff files
echo ""
mecho "Validating GFF files:"
echo ""

bash scripts/validate_gff.sh $T_QUERY \
    $T_GFF \
    ragoo2_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragoo2_output_tomato_default/$T_GFF_PREF.corr.gff

bash scripts/validate_gff.sh ragoo2_output_tomato_default/$T_QUERY_PREF.corrected.fasta \
    ragoo2_output_tomato_default/$T_GFF_PREF.corr.gff \
    ragoo2_output_tomato_default/ragoo2.scaffolds.fasta \
    ragoo2_output_tomato_default/$T_GFF_PREF.scaf.gff





# Run RaGOO2
# Settings: with nucmer
# Data: Arabidopsis
echo ""
echo "*************************************************************"
echo "***       Running RaGOO2 on tomato data with Nucmer       ***"
echo "*************************************************************"
echo ""

bash scripts/run_nucmer.sh $A_REF \
    $A_QUERY \
    $A_GFF \
    ragoo2_output_Ara_nucmer

# Validate the agp files
echo ""
mecho "Validating AGP files and associated fasta files:"
echo ""

bash scripts/validate_agp.sh $A_QUERY \
    ragoo2_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragoo2_output_Ara_nucmer/ragoo2.correction.agp

bash scripts/validate_agp.sh ragoo2_output_Ara_nucmer/ragoo2.scaffolds.fasta \
    ragoo2_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragoo2_output_Ara_nucmer/ragoo2.scaffolds.agp

# Validate the gff files
echo ""
mecho "Validating GFF files:"
echo ""

bash scripts/validate_gff.sh $A_QUERY \
    $A_GFF \
    ragoo2_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragoo2_output_Ara_nucmer/$A_GFF_PREF.corr.gff

bash scripts/validate_gff.sh ragoo2_output_Ara_nucmer/$A_QUERY_PREF.corrected.fasta \
    ragoo2_output_Ara_nucmer/$A_GFF_PREsF.corr.gff \
    ragoo2_output_Ara_nucmer/ragoo2.scaffolds.fasta \
    ragoo2_output_Ara_nucmer/$A_GFF_PREF.scaf.gff

# Validate the unique anchor filtering
echo ""
mecho "Validating alignment filtering:"
echo ""

bash scripts/validate_uaf.sh ragoo2_output_Ara_nucmer/c_query_against_ref.delta \
    ragoo2_output_Ara_nucmer/ragoo2.correction.debug.filtered.paf \
    1000

bash scripts/validate_uaf.sh ragoo2_output_Ara_nucmer/query_against_ref.delta \
    ragoo2_output_Ara_nucmer/ragoo2.scaffolds.debug.filtered.paf \
    1000