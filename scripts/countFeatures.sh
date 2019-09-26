#!/bin/bash

set -eo pipefail

declare -a DEPENDENCIES=(featureCounts samtools bedtools bc mktemp parallel datamash Rscript dstat)
for DEP in ${DEPENDENCIES[@]}; do $(command -v ${DEP} > /dev/null 2> /dev/null) || (echo "Cannot find ${DEP}" >&2; exit 1); done

#   Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -l|--sample-list -a|--annotation -o|--outdir -p|--project [-s|--structural] [-e|--expression] \n\
Where:  -l|--sample-list is a list of proecessed BAM files to analyze \n\
        -a|--annotation is a reference annotation file in GTF/GFF3 format \n\
        -o|--outdir is an output directory to use \n\
            will be modified to \${OUTDIR}/Quantify_Summarize \n\
        -p|--project is a project name for output naming \n\
        [-s|--structural] is an optional annotation (GTF/GFF3/BED) file with structural RNA annotations \n\
            if not passed, percent structural RNA will not be included in output file \n\
        [-e|--expression] is an optional reference expression file for dstat calculation \n\
            the first column should be the gene ID and all subsequent columns should be median/average \n\
            expression values for the reference population (eg. median expression in GTEx tissues) \n\
            a header row is expected for this file
" >&2
    exit 1
}

export -f Usage

#   Calculate percent duplicated
function PercDup() {
    local bamfile=$1 # What BAM file are we working with
    local bamname=$(basename ${bamfile} .bam)
    local tempfile=$(mktemp)
    (set -x; samtools rmdup -S ${bamfile} ${tempfile})
    local ndedup=$(samtools view -F 0x0100 -c ${tempfile})
    local ntotal=$(samtools view -F 0x0100 -c ${bamfile})
    local perc=$(echo "scale=4; (${ndedup} / ${ntotal})" | bc -l)
    (set -x; rm -f ${tempfile})
    echo -e "${bamname},${perc}"
}

export -f PercDup

#   Calculate percent structural
function PercStruct() {
    local bamfile=$1
    local structural=$2
    local bamname=$(basename ${bamfile} .bam)
    local nstruct=$(set -x; bedtools intersect -a ${bamfile} -b ${structural} | samtools view -F 0x0100 -c -)
    local ntotal=$(set -x; samtools view -F 0x0100 -c ${bamfile})
    local perc=$(echo "scale=4; (${nstruct} / ${ntotal})" | bc -l)
    echo -e "${bamname},${perc}"
}

export -f PercStruct

#   Calculate expression diversity
function ExprDiversity() {
    local bamfile=$1
    local counts=$2
    local bamname=$(basename ${bamfile} .bam)
    local length=$(grep -v '#' ${counts} | head -1 | tr '\t' '\n' | grep -n 'Length' | cut -f 1 -d ':')
    local -a lengths=($(grep -v '#' ${counts} | tail -n +2 | cut -f ${length}))
    local index=$(grep -v '#' ${counts} | head -1 | tr '\t' '\n' | grep -n ${bamfile} | cut -f 1 -d ':')
    local -a sample_counts=($(grep -v '#' ${counts} | tail -n +2 | cut -f ${index}))
    local depth=$(echo ${sample_counts[@]} | tr ' ' '\n' | datamash sum 1)
    local tempfile=$(mktemp)
    paste <(echo ${lengths[@]} | tr ' ' '\n') <(echo ${sample_counts[@]} | tr ' ' '\n') > ${tempfile}
    local expr_count=$(set -x; Rscript -e "cat(sum(apply(X = read.table('${tempfile}', header = FALSE), MARGIN = 1, FUN = function(x) x[2] / (x[1] / 1000) / (${depth} / 1e6)) > 0.1))")
    (set -x; rm -f ${tempfile})
    echo -e "${bamname},${expr_count}"
}

export -f ExprDiversity

while [[ $# -gt 1 ]]; do
    KEY=$1
    case $KEY in
        -l|--sample-list) # Sample list
            SAMPLE_LIST=$2
            shift
            ;;
        -a|--annotation) # Reference annotation in GTF/GFF format
            ANNOTATION=$2
            shift
            ;;
        -s|--structural) # Annotation of structural RNA
            STRUCTURAL=$2
            shift
            ;;
        -e|--expression) # Reference expression for dstat
            MEDIAN_EXPR=$2
            shift
            ;;
        -o|--outdir) # Output directory
            OUTDIR="${2}/Quantify_Summarize"
            shift
            ;;
        -p|--project) # Project name
            PROJECT=$2
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

[[ -z "${SAMPLE_LIST}" || -z "${ANNOTATION}" || -z "${OUTDIR}" || -z "${PROJECT}" ]] && Usage
[[ -f "${ANNOTATION}" ]] || (echo "Cannot find annotation file ${ANNOTATION}" >&2; exit 1)
[[ -f "${SAMPLE_LIST}" ]] || (echo "Cannot find sample list ${SAMPLE_LIST}" >&2; exit 1)

THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)
(set -x; mkdir -p "${OUTDIR}")

declare -a BAMS=($(<${SAMPLE_LIST}))
[[ "${#BAMS[@]}" -lt 1 ]] && (echo "No samples found in sample list ${SAMPLE_LIST}" >&2; exit 1)
for sample in ${BAMS[@]}; do
    [[ -f "${sample}" ]] || (echo "Cannot find sample ${sample}" >&2; exit 1)
done

#   Generate feature counts
echo "Counting features in ${#BAMS[@]} BAM files" >&2
(set -x; featureCounts --primary --verbose -T "${THREADS}" -a "${ANNOTATION}" -o "${OUTDIR}/${PROJECT}.counts" ${BAMS[@]})

#   Generate percent duplicated
echo "Calculating percent duplicated" >&2
declare -A DUPLICATED
declare -a HOLDING=($(parallel --verbose --keep-order "PercDup {}" ::: ${BAMS[@]}))
for i in ${HOLDING[@]}; do DUPLICATED[$(echo $i | cut -f 1 -d ',')]=$(echo $i | cut -f 2 -d ','); done
unset HOLDING

#   Calculate expression diversity
echo "Calculating expression diversity" >&2
declare -A EXPR_DIV
declare -a HOLDING=($(parallel --verbose --keep-order "ExprDiversity {} ${OUTDIR}/${PROJECT}.counts" ::: ${BAMS[@]}))
for i in ${HOLDING[@]}; do EXPR_DIV[$(echo $i | cut -f 1 -d ',')]=$(echo $i | cut -f 2 -d ','); done
unset HOLDING

#   Calculate dstat scores
if [[ -z "${MEDIAN_EXPR}" ]]; then
    echo "Not calculating dstat scores" >&2
else
    [[ -f "${MEDIAN_EXPR}" ]] || (echo "Cannot find median expression reference file ${MEDIAN_EXPR}" >&2; exit 1)
    (set -x; dstat -c "${OUTDIR}/${PROJECT}.counts" -t "${MEDIAN_EXPR}" -o "${OUTDIR}/${PROJECT}_dstat.tsv")
fi

#   Get percent structural RNA
if [[ -z "${STRUCTURAL}" ]]; then
    echo "Not calculating percent structrual" >&2
else
    [[ -f "${STRUCTURAL}" ]] || (echo "Cannot find structural RNA BED file ${STRUCTURAL}" >&2; exit 1)
    declare -A PSTRUCT
    declare -a HOLDING=($(parallel --verbose --keep-order "PercStruct {} ${STRUCTURAL}" ::: ${BAMS[@]}))
    for i in ${HOLDING[@]}; do PSTRUCT[$(echo $i | cut -f 1 -d ',')]=$(echo $i | cut -f 2 -d ','); done
    unset HOLDING
fi

#   Assemble table
OFILE="${OUTDIR}/${PROJECT}_stats.tsv"
declare -a COLNAMES=(p_dup expr_div)

## Figure out if we're adding dstat scores
if [[ ! -z ${MEDIAN_EXPR} ]]; then
    declare -a TISSUES=($(cut -f 1 "${OUTDIR}/${PROJECT}_dstat.tsv" | grep -vE '^$' | tr -d ' '))
    for t in ${TISSUES[@]}; do
        COLNAMES+=($(echo -n "dstat_${t}"))
    done
fi

## Figure our if we're adding percent structural RNA
[[ -z "${STRUCTURAL}" ]] || COLNAMES+=(p_struct)

## Write the table
echo -e sample$(for i in ${COLNAMES[@]}; do echo -n "\t${i}"; done) > ${OFILE}
for sample in ${!DUPLICATED[@]}; do
    declare -a ROW=(${sample} ${DUPLICATED[${sample}]} ${EXPR_DIV[${sample}]})
    if [[ ! -z "${MEDIAN_EXPR}" ]]; then
        INDEX=$(head -1 "${OUTDIR}/${PROJECT}_dstat.tsv" | tr '\t' '\n' | grep -n ${sample} | cut -f 1 -d ':')
        ROW+=($(cut -f ${INDEX} "${OUTDIR}/${PROJECT}_dstat.tsv" | tail -n +2))
        unset INDEX
    fi
    [[ -z "${STRUCTURAL}" ]] || ROW+=(${PSTRUCT[${sample}]})
    echo ${ROW[@]} | tr ' ' '\t' >> ${OFILE}
    unset ROW
done
