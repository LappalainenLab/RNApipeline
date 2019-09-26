#!/bin/bash

#   Run STAR 2map

set -eo pipefail

#   Check dependencies
declare -a DEPENDENCIES=(STAR zcat bzcat)
for dep in ${DEPENDENCIES[@]}; do $(command -v ${dep} > /dev/null 2> /dev/null) || (echo "Cannot find ${dep}" >&2; exit 1); done

#   Defaults
THREADS_DEFAULT='1'
GENOME_DEFAULT='/gpfs/commons/home/scastel/references/fasta/GRCh37.fa'
INDEX_DEFAULT='/gpfs/commons/home/scastel/references/star_noj_GRCh37'
PROJECT_DEFAULT='starMap'

#   Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -f|--forward [-r|--reverse] [-i|--index] [-g|--genome] [-t|--threads] [-n|--sample-name] [-p|--project] [-o|--outdir] \n\
Where:  -f|--forward is a forward or single-ended FASTQ file \n\
        [-r|--reverse] is an optional reverse FASTQ file \n\
        [-g|--genome] is an optional genomic FASTA file \n\
            defaults to '${GENOME_DEFAULT}' \n\
        [-i|--index] is an optional STAR index directory \n\
            defaults to '${INDEX_DEFAULT}' \n\
        [-t|--threads] is an optional specifier for the number of threads \n\
            defaults to '${THREADS_DEFAULT}' \n\
        [-n|--sample-name] is an optional sample name \n\
            defaults to the basename of the forward/single-ended FASTQ file \n\
            use to customize output naming of SAM files \n\
        [-p|--project] is an optional project name for final sample list \n\
            defaults to '${PROEJCT_DEFAULT}' \n\
        [-o|--outdir] is an optional output directory \n\
            will be modified to \${OUTDIR}/Read_Mapping \n\
            defaults to '$(pwd -P)' \n\
" >&2
    exit 1
}

export -f Usage

#   Function to make a genome directory
function MakeGenomeDirectory() {
    local refGen="$1"; shift # What reference genome are we using?
    local genomeDir="$1"; shift # Where do we make the genome directory?
    local threads="$1"; shift # How many threads are we using?
    local -a extra=("${@}") # Anything else
    (set -x; mkdir -p "${genomeDir}")
    [[ "${genomeDir:$((${#genomeDir} - 1))}" != '/' ]] && local genomeDir="${genomeDir}/"
    (
        set -x
        STAR --runMode genomeGenerate \
            --genomeDir "${genomeDir}" \
            --genomeFastaFiles "${refGen}" \
            --runThreadN "${threads}" \
            --outTmpDir "${genomeDir}tmp" \
            ${extra[@]}
    )
}

export -f MakeGenomeDirectory

#   Function to check compression levels
function CheckCompression() {
    local file="$1" # What file are we working with?
    local extension="$(basename ${file} | rev | cut -f 1 -d '.' | rev)"
    case "${extension}" in
        gz) # gzipped
            echo "Working with gzipped FASTQs" >&2
            local cmd='--readFilesCommand zcat'
            ;;
        bz2) # bzipped
            echo "Working with bzipped FASTQs" >&2
            local cmd='--readFilesCommand bzcat'
            ;;
        *) # Assume uncompressed
            echo "Neither gzipped nor bzipped, assuming uncompressed" >&2
            local cmd=''
            ;;
    esac
    echo "${cmd}"
}

export -f CheckCompression

#   Ensure we have the proper number of arguments
[[ "$#" -lt 4 ]] && Usage

#   Parse arguments
while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -f|--forward) # Forward read
            FWD="$2"
            shift
            ;;
        -r|--reverse) # Reverse read
            REV="$2"
            shift
            ;;
        -i|--index) # Genome directory
            GEN_DIR="$2"
            shift
            ;;
        -g|--genome) # Reference fasta
            REF_GEN="$2"
            shift
            ;;
        -t|--threads) # Number of threads
            THREADS="$2"
            shift
            ;;
        -n|--sample-name) # Basename of sample
            SAMPLE_NAME="$2"
            shift
            ;;
        -p|--project) # Project name
            PROJECT="$2"
            shift
            ;;
        -o|--outdir) # Set output directory
            OUTDIR="${2}/Read_Mapping"
            shift
            ;;
        *) # Anything else
            Usage
            ;;
    esac
    shift
done

#   Check arguments
[[ -z "${REF_GEN}" ]] && REF_GEN="${GENOME_DEFAULT}"
[[ -z "${GEN_DIR}" ]] && GEN_DIR="${INDEX_DEFAULT}"
[[ -z "${FWD}" || -z "${REF_GEN}" ]] && Usage
[[ -z "${OUTDIR}" ]] && OUTDIR="$(pwd -P)/Read_Mapping"
[[ -f "${FWD}" ]] || (echo "Cannot find input file ${FWD}" >&2; exit 1)
[[ ! -z "${REV}" && ! -f "${REV}" ]] && (echo "Cannot find reverse file ${REV}" >&2; exit 1)
[[ -z "${THREADS}" ]] && THREADS="${THREADS_DEFAULT}"
[[ -z "${SAMPLE_NAME}" ]] && SAMPLE_NAME="$(basename ${FWD} .$(basename ${FWD} | cut -f 2- -d '.'))"
[[ -z "${PROJECT}" ]] && PROJECT="${PROJECT_DEFAULT}"

#   Check genome directory
if [[ -z "${GEN_DIR}" || ! -d "${GEN_DIR}" ]]; then
    [[ -d "${GEN_DIR}" ]] || echo "Cannot find genome directory ${GEN_DIR}" >&2
    REF_BASE="$(basename ${REF_GEN} .$(echo ${REF_GEN} | rev | cut -f 1 -d '.' | rev))"
    GEN_DIR="$(pwd -P)/${REF_BASE}_Index"
    MakeGenomeDirectory "${REF_GEN}" "${GEN_DIR}" "${THREADS}"
else
    echo "Assuming complete genome directory" >&2
fi

#   Check single- vs. paired-end mode
READ_FILES="${FWD}"
OUT_BASE="${OUTDIR}/${SAMPLE_NAME}/"
PASS_ONE="${OUT_BASE}Round1/"
TEMP_ONE="${PASS_ONE}tmp/"
(set -x; mkdir -p "${OUT_BASE}" "${PASS_ONE}")
[[ -d "${TEMP_ONE}" ]] && (set -x; rm -rf "${TEMP_ONE}")
READ_COMPRESSION="$(CheckCompression ${FWD})"

if [[ ! -z "${REV}" ]]; then
    READ_FILES="${READ_FILES} ${REV}"
    REV_COMPRESSION="$(CheckCompression ${REV})"
    [[ "${READ_COMPRESSION}" != "${REV_COMPRESSION}" ]] && (echo "The two files are not equally compressed, please fix" >&2, exit 1)
fi

#   First pass
(
    set -x
    STAR --genomeDir "${GEN_DIR}" \
        --readFilesIn ${READ_FILES} \
        --runThreadN "${THREADS}" \
        --outFileNamePrefix "${PASS_ONE}" \
        --outTmpDir "${TEMP_ONE}" \
        ${READ_COMPRESSION}
)

#   New genome
SECOND_GEN="${OUT_BASE}GenDir"
MakeGenomeDirectory "${REF_GEN}" "${SECOND_GEN}" "${THREADS}" --sjdbFileChrStartEnd "${PASS_ONE}SJ.out.tab" --sjdbOverhang 75

#   Second pass
echo "Second pass" >&2
PASS_TWO="${OUT_BASE}Round2/"
TEMP_TWO="${PASS_TWO}tmp/"
(set -x; mkdir -p ${PASS_TWO})
[[ -d ${TEMP_TWO} ]] && (set -x; rm -rf ${TEMP_TWO})
(
    set -x
    STAR --genomeDir ${SECOND_GEN} \
        --readFilesIn ${READ_FILES} \
        --runThreadN ${THREADS} \
        --outFileNamePrefix ${PASS_TWO} \
        --outTmpDir ${TEMP_TWO} \
        ${READ_COMPRESSION}
)

(set -x; rm -rf ${TEMP_ONE} ${TEMP_TWO})

(set -x; mv "${PASS_TWO}/Aligned.out.sam" "${PASS_TWO}/${SAMPLE_NAME}.out.sam")
echo "${PASS_TWO}/${SAMPLE_NAME}.out.sam" >> "${OUTDIR}/${PROJECT}_mapped.txt"
