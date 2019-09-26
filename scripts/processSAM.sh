#!/bin/bash

set -eo pipefail

#   Check for dependencies
declare -a DEPENDENCIES=(samtools java)
for dep in ${DEPENDENCIES[@]}; do $(command -v ${dep} > /dev/null 2> /dev/null) || (echo "Cannot find ${dep}" >&2; exit 1); done

#   Defaults
PICARD_DEFAULT='/gpfs/commons/groups/lappalainen_lab/scastel/projects/tl_ase_pipeline/software/picard-tools-1.124/picard.jar'
REF_DEFAULT='/gpfs/commons/home/scastel/references/fasta/GRCh37.fa'
MEM_DEFAULT='7000m'
OUT_DEFAULT="$(pwd -P)/SAM_Processing"
PROJECT_DEFAULT="SAM_Processing"

#   Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -s|--sam-file [-r|--reference] [-o|--outdir] [-p|--picard] [-m|--max-mem] [--csi] [-p|--project] \n\
Where:  -s|--sam is the SAM file to process \n\
        -r|--reference is an optional reference genome \n\
            defaults to '${REF_DEFAULT}' \n\
        -o|--outdir is an optional output directory \n\
            will be modified to \${OUTDIR}/SAM_Processing n\
            defaults to '$(dirname ${OUT_DEFAULT})' \n\
        -p|--picard is an optional path to the Picard JAR file \n\
            defaults to '${PICARD_DEFAULT}' \n\
        -m|--max-mem is an optional maximal amount of memory for Picard \n\
            defaults to '${MEM_DEFAULT}' \n\
        --csi is an optional flag to specify CSI indices \n\
            if not passed, makes BAI indices \n\
        -p|--project is an optional project name for the final sample list \n\
            defaults to '${PROJECT_DEFAULT}' \n\
" >&2
    exit 1
}

export -f Usage

#   Check that we have enough arguments
[[ "$#" -lt 1 ]] && Usage

#   Parse the arguments
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -s|--sam-file) # SAM file
            SAM="$2"
            shift
            ;;
        -r|--reference) # Referene genome
            REF_GEN="$2"
            shift
            ;;
        -o|--outdir) # Output directory
            OUTDIR="${2}/SAM_Processing"
            shift
            ;;
        --csi) # Use CSI index
            INDEX_TYPE='-c'
            ;;
        -p|--picard) # Picard JAR
            PICARD="$2"
            shift
            ;;
        -m|--max-mem) # Max memory
            MAX_MEM="$2"
            shift
            ;;
        --project) # Project name
            PROJECT="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

#   Check the arguments
[[ -z "${REF_GEN}" ]] && REF_GEN="${REF_DEFAULT}"
[[ -z "${SAM}" || -z "${REF_GEN}" ]] && Usage
[[ -f "${SAM}" ]] || (echo "Cannot find SAM file ${SAM}" >&2; exit 1)
[[ -f "${REF_GEN}" ]] || (echo "Cannot find reference genome ${REF_GEN}" >&2; exit 1)
[[ -z "${PICARD}" ]] && PICARD="${PICARD_DEFAULT}"
[[ -f "${PICARD}" ]] || (echo "Cannot find Picard JAR file ${PICARD}" >&2; exit 1)
[[ -z "${OUTDIR}" ]] && OUTDIR="${OUT_DEFAULT}"
[[ -z "${MAX_MEM}" ]] && MAX_MEM="${MEM_DEFAULT}"
[[ -z "${INDEX_TYPE}" ]] && INDEX_TYPE='-b'
[[ -z "${PROJECT}" ]] && PROJECT="${PROJECT_DEFAULT}"

#   Check for faidx
[[ -f "${REF_GEN}.fai" ]] || (echo "Generating faidx for reference genome" >&2; set -x; samtools faidx "${REF_GEN}")

#   Create some output names
[[ "${OUTDIR:$((${#OUTDIR} - 1))}" == '/' ]] && OUTDIR="${OUTDIR:0:$((${#OUTDIR} - 1))}"
RAW_DIR="${OUTDIR}/raw"
SORT_DIR="${OUTDIR}/sorted"
RG_DIR="${OUTDIR}/addrg"
FIN_DIR="${OUTDIR}/finished"
SAM_BASE="$(basename ${SAM} .sam)"
(set -x; mkdir -p "${RAW_DIR}" "${SORT_DIR}" "${RG_DIR}" "${FIN_DIR}/metrics")

#   Convert from SAM to BAM
(set -x; samtools view -bhT "${REF_GEN}" "${SAM}" > "${RAW_DIR}/${SAM_BASE}_raw.bam")

#   Sort the BAM file
(set -x; samtools sort "${RAW_DIR}/${SAM_BASE}_raw.bam" > "${SORT_DIR}/${SAM_BASE}_sorted.bam")

#   Add read groups
(
    set -x
    java -Xmx"${MAX_MEM}" -jar "${PICARD}" \
        AddOrReplaceReadGroups \
        I="${SORT_DIR}/${SAM_BASE}_sorted.bam" \
        O="${RG_DIR}/${SAM_BASE}_addrg.bam" \
        RGID=id \
        RGLB=library \
        RGPL=platform \
        RGPU=machine \
        RGSM=sample
)

#   Mark duplicates
(
    set -x
    java -Xmx"${MAX_MEM}" -jar "${PICARD}" \
        MarkDuplicates \
        I="${RG_DIR}/${SAM_BASE}_addrg.bam" \
        O="${FIN_DIR}/${SAM_BASE}_finished.bam" \
        VALIDATION_STRINGENCY=SILENT \
        M="${FIN_DIR}/metrics/${SAM_BASE}.metrics"
)

#   Create index
case "${INDEX_TYPE}" in
    -b)
        echo "Using BAI indices" >&2
        ;;
    -c)
        echo "Using CSI indices" >&2
        ;;
    *)
        echo "Something happened" >&2
        exit 1
        ;;
esac

(set -x; samtools index "${INDEX_TYPE}" "${FIN_DIR}/${SAM_BASE}_finished.bam")
echo "${FIN_DIR}/${SAM_BASE}_finished.bam" >> "${OUTDIR}/${PROJECT}_bams.txt"
