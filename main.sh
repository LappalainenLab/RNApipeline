#!/bin/bash

set -eo pipefail

$(command -v qsub > /dev/null 2> /dev/null) || (echo "This tool submits jobs to a queue, please run on a cluster of some sort"; exit 1)

#   Usage message
function Usage() {
    echo -e "\
Usage:  $(basename $0) <handler> proj.conf \n\
Where:  <handler> is one of:
            1 or Quality_Assessment \n\
            2 or Sequence_Trimming \n\
            3 or Read_Mapping \n\
            4 or SAM_Processing \n\
            5 or Gene_Counting \n\
And:    proj.conf is the full file path to the configuration file
" >&2
    exit 1
}

#   Export the function
export -f Usage

#   Function to find maximum memory
function MaxMem() {
    local qsub="$1"
    local memraw=$(echo "${qsub}" | grep -oE 'mem=[[:alnum:]]+' | cut -f 2 -d '=')
    local digits=$(echo "${memraw}" | -oE '[[:digit:]]+')
    if $(echo "${memraw}" | grep -i 'g' > /dev/null 2> /dev/null); then
        local maxmem="${digits}G"
    elif $(echo "${memraw}" | grep -i 'm' > /dev/null 2> /dev/null); then
        local maxmem="${digits}M"
    elif $(echo "${memraw}" | grep -i 'k' > /dev/null 2> /dev/null); then
        local maxmem="${digits}K"
    else
        local maxmem="${digits}"
    fi
    echo "${maxmem}"
}

#   Export the function
export -f MaxMem

#   Function to find number of threads requested
function Nthreads() {
    local qsub="$1"
    local nthreads=$(echo "${qsub}" | grep -oE 'pe smp [[:digit:]]+' | cut -f 3 -d ' ')
    echo "${nthreads}"
}

#   Export the function
export -f Nthreads

#   Check to make sure our samples exist
function checkSamples() {
    local sample_list="$1" # Sample ist
    [[ -f "${sample_list}" ]] || (echo "Cannot find sample list ${sample_list}" >&2; exit 1)
    for sample in $(<${sample_list}); do [[ -f "${sample}" ]] || (echo "Cannot find sample ${sample}" >&2; exit 1); done
}

#   Export the function to be used elsewhere
export -f checkSamples

#   Where is 'sequence_handling' located?
SEQUENCE_HANDLING=$(pwd -P)

#   For clusters that don't use login shells
$(type module > /dev/null 2> /dev/null) || function module() { eval $(/usr/bin/modulecmd bash $*); }
export -f module

#   If we have less than two arguments
if [[ "$#" -lt 2 ]]; then Usage; fi # Display the usage message and exit

ROUTINE="$1" # What routine are we running?
CONFIG="$2" # Where is our config file?

#   Find full path to the config file
#   Because SGE is stupid...
CONFIG="$(dirname ${CONFIG})/$(basename ${CONFIG})"
#   For the home directory
CONFIG="${CONFIG/'~/'/${HOME}/}"
#   For a directory above us
CONFIG="${CONFIG/'../'/$(pwd -P | rev | cut -f 2- -d '/' | rev)/}"
#   For this directory
CONFIG="${CONFIG/'./'/$(pwd -P)/}"

#   Check the existence of the config file
[[ -f "${CONFIG}" ]] || (echo "Cannot find configuration file" >&2; exit 1)
source "${CONFIG}" # Source it, providing parameters and software
bash "${CONFIG}" > /dev/null 2> /dev/null # Load any modules

#   Make our output directory
[[ -z "${OUT_DIR}" ]] && OUT_DIR="${SEQUENCE_HANDLING}/${PROJECT}"
mkdir -p "${OUT_DIR}"

[[ -z "${EMAIL}" ]] && QSUB_EMAIL='' || QSUB_EMAIL="-m abe -M ${EMAIL}"

case "${ROUTINE}" in
    1|Quality_Assessment)
        echo "$(basename $0): Assessing quality..." >&2
        checkSamples ${RAW_SAMPLES}
        echo "${SEQUENCE_HANDLING}/scripts/assessQuality.sh --sample-list ${RAW_SAMPLES} --outdir ${OUT_DIR} --project ${PROJECT}" | qsub ${QA_QSUB} ${QSUB_EMAIL} -N Quality_Assessment
        ;;
    2|Sequence_Trimming)
        echo "$(basename $0): Trimming read sequences..." >&2
        checkSamples ${RAW_SAMPLES}
        #   Partition samples into forward, reverse, and single-ended samples
        if [[ -z "${FORWARD_NAMING}" || -z "${REVERSE_NAMING}" ]]; then
            declare -a FORWARD=()
            declare -a REVERSE=()
            declare -a SINGLES=($(<${RAW_SAMPLES}))
        else
            declare -a FORWARD=($(grep -E "${FORWARD_NAMING}" "${RAW_SAMPLES}" | sort))
            declare -a REVERSE=($(grep -E "${REVERSE_NAMING}" "${RAW_SAMPLES}" | sort))
            declare -a SINGLES=($(grep -vEf <(echo -e "${FORWARD_NAMING}\n${REVERSE_NAMING}") "${RAW_SAMPLES}" | sort))
        fi
        $(${PHRED64}) && QUAL='-64' || QUAL=''
        #   Submit paired samples
        [[ "${#FORWARD[@]}" -ne "${#REVERSE[@]}" ]] && (echo "Unequal numbers of forward and reverse samples" >&2; exit 1)
        for i in $(seq ${#FORWARD[@]}); do
            NAME=$(basename ${FORWARDS[$(($i - 1))]} ${FORWARD_NAMING})
            echo "${SEQUENCE_HANDLING}/scripts/trimFASTQ.sh --forward ${FORWARD[$(($i - 1))]} --reverse ${REVERSE[$(($i - 1))]} --outdir ${OUT_DIR} --adapters ${ADAPTERS} --project ${PROJECT} ${QUAL}" | qsub ${ST_QSUB} ${QSUB_EMAIL} -N "Trim_${NAME}_paired"
        done
        #   Submit single samples
        for i in $(seq ${#SINGLES[@]}); do
            NAME=$(basename ${SINGLES[$((i - 1))]})
            echo "${SEQUENCE_HANDLING}/scripts/trimFASTQ.sh --forward ${SINGLES[$(($i - 1))]} --outdir ${OUT_DIR} --adapters ${ADAPTERS} ${QUAL}" | qsub ${ST_QSUB} ${QSUB_EMAIL} -N "Trim_${NAME}"
        done
        ;;
    3|Read_Mapping)
        echo "$(basename $0): Mapping reads..." >&2
        checkSamples "${TRIMMED_LIST}"
        NTHREADS=$(Nthreads "${RM_QSUB}")
        #   Partition samples into forward, reverse, and single-ended reads
        declare -a FORWARD=($(grep -E "${FORWARD_TRIMMED}" ${TRIMMED_LIST} | sort))
        declare -a REVERSE=($(grep -E "${REVERSE_TRIMMED}" ${TRIMMED_LIST} | sort))
        declare -a SINGLES=($(grep -E "${SINGLES_TRIMMED}" ${TRIMMED_LIST} | sort))
        #   Submit paired samples
        [[ "${#FORWARD[@]}" -ne "${#REVERSE[@]}" ]] && (echo "Unequal numbers of forward and reverse samples" >&2; exit 1)
        for i in $(seq ${#FORWARD[@]}); do
            NAME=$(basename ${FORWARD[$(($i - 1))]} ${FORWARD_TRIMMED})
            echo "${SEQUENCE_HANDLING}/scripts/starMap.sh --forward ${FORWARD[$(($i - 1))]} --reverse ${REVERSE[$(($i - 1))]} --index ${REF_IND} --genome ${REF_GEN} --threads ${NTHREADS} --sample-name ${NAME} --project ${PROJECT} --outdir ${OUT_DIR}" | qsub ${RM_QSUB} ${QSUB_EMAIL} -N "Map_${NAME}_paired"
        done
        #   Submit single samples
        for i in $(seq ${#SINGLES[@]}); do
            NAME=$(basename ${SINGLES[$(($i - 1))]} ${FORWARD_TRIMMED})
            echo "${SEQUENCE_HANDLING}/scripts/starMap.sh --forward ${SINGLES[$(($i - 1))]} --index ${REF_IND} --genome ${REF_GEN} --threads ${NTHREADS} --sample-name ${NAME} --project ${PROJECT} --outdir ${OUT_DIR}" | qsub ${RM_QSUB} ${QSUB_EMAIL} -N "Map_${NAME}_single"
        done
        ;;
    4|SAM_Processing)
        echo "$(basename $0): Processing SAM files" >&2
        checkSamples "${MAPPED_LIST}"
        MAXMEM=$(MaxMem "${SP_QSUB}")
        case ${INDEX_TYPE} in
            BAI)
                INDEXCMD=''
            CSI)
                INDEXCMD='--csi'
            *)
                echo "Unknown index type ${INDEX_TYPE}, please choose from 'BAI' or 'CSI'"
                exit 1
                ;;
        esac
        for sample in $(<${MAPPED_LIST}); do
            echo "${SEQUENCE_HANDLING}/scripts/processSAM.sh --sam-file ${sample} --reference ${REF_GEN} --outdir ${OUT_DIR} --max-mem ${MAXMEM} --project ${PROJECT} ${INDEXCMD}" | qsub ${SP_QSUB} ${QSUB_EMAIL} -N "Process_$(basename ${sample} .sam)"
        done
        ;;
    5|Gene_Counting)
        echo "$(basename $0): Counting genes..." >&2
        checkSamples "${BAM_LIST}"
        echo "${SEQUENCE_HANDLING}/scripts/countFeatures.sh --sample-list ${BAM_LIST} --annotation ${REF_ANN} --output ${OUT_DIR}/Gene_Counting/${PROJECT}.counts" | qsub ${GC_QSUB} ${QSUB_EMAIL} -N "Gene_Counting"
        ;;
    Variant_Calling)
        echo "$(basename $0): Calling variants..." >&2
        echo "Variant_Calling is not yet implemented, exiting..." >&2; exit 1
        ;;
    *)
        Usage
        ;;
esac
