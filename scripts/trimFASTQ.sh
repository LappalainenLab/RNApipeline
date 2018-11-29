#!/bin/bash

set -eo pipefail

#   Check dependencies
$(command -v java > /dev/null 2> /dev/null) || (echo "Cannot find Java" >&2; exit 1)

#   Set some defaults
TRIMMOMATIC='/nfs/sw/trimmomatic/trimmomatic-0.36/trimmomatic-0.36.jar'
ADAPTERS_DEFAULT='/nfs/sw/trimmomatic/trimmomatic-0.36/adapters/NexteraPE-PE.fa'
OUTDIR_DEFAULT="$(pwd -P)/Sequence_Trimming"
PROJECT_DEFAULT='Trimmomatic'

#   Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -f|--forward [-r|--reverse] [-o|--outdirectory] [-a|--adapters] [-s|--split] [-64] \n\
" >&2
    exit 1
}

#   Export the function
export -f Usage

#   Ensure we have the proper number of arguments
[[ "$#" -lt 1 ]] && Usage

#   Parse arguments
while [[ "$#" -ge 1 ]]; do
    case $1 in
        -f|--forward)
            FORWARD="$2"
            shift
            ;;
        -r|--reverse)
            REVERSE="$2"
            shift
            ;;
        -o|--outdirectory)
            OUTDIR="${2}/Sequence_Trimming"
            shift
            ;;
        -a|--adapters)
            ADAPTERS="$2"
            shift
            ;;
        -s|--split)
            SPLIT="$2"
            shift
            ;;
        -p|--project)
            PROJECT="$2"
            shift
            ;;
        -64)
            QUALITY='phred64'
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

#   Check arguments
[[ -z "${FORWARD}" ]] && Usage
[[ -z "${REVERSE}" ]] && MODE='SE' || MODE='PE'
[[ -z "${OUTDIR}" ]] && OUTDIR="${OUTDIR_DEFAULT}"
[[ -z "${ADAPTERS}" ]] && ADAPTERS="${ADAPTERS_DEFAULT}"
[[ -z "${QUALITY}" ]] && QUALITY='phred33'

[[ -f "${FORWARD}" ]] || (echo "Cannot find forward FASTQ ${FORWARD}" >&2; exit 1)
[[ "${MODE}" == 'PE' && ! -f "${REVERSE}" ]] && (echo "Cannot find reverse FASTQ ${REVERSE}" >&2; exit 1)
[[ -f "${ADAPTERS}" ]] || (echo "Cannot find adapters file ${ADAPTERS}" >&2; exit 1)

#   Make an output directory
(set -x; mkdir -p "${OUTDIR}")

#   Figure out our outbase
if [[ -z "${SPLIT}" ]]; then
    EXTENSION="$(echo ${FORWARD} | rev | cut -f 1 -d '.' | rev)"
    [[ "${EXTENSION}" == 'gz' ]] && EXTENSION=".$(echo ${FORWARD} | rev | cut -f 1,2 -d '.' | rev)"
    OUT_BASE=$(basename ${FORWARD} ${EXTENSION})
else
    OUT_BASE="$(basename ${FORWARD})"
    OUT_BASE="${OUT_BASE/${SPLIT}*/}"
fi

#   Make an array of output file names
declare -a OUTPUTS=()
if [[ "${MODE}" == 'PE' ]]; then
    for direction in forward reverse; do
        for pair in paired unpaired; do
            OUTPUTS+=("${OUTDIR}/${OUT_BASE}_${direction}_${pair}.fastq.gz")
        done
    done
else
    OUTPUTS+=("${OUTDIR}/${OUT_BASE}_trimmed.fastq.gz")
fi

#   Run Trimmomatic
(
    set -x;
    java -jar "${TRIMMOMATIC}" \
        "${MODE}" \
        -"${QUALITY}" \
        "${FORWARD}" \
        "${REVERSE}" \
        ${OUTPUTS[@]} \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
)

echo "${OUTPUTS[@]}" | tr ' ' '\n' | grep -E 'paired\.fastq\.gz|trimmed\.fastq\.gz' >> "${OUTDIR}/${PROJECT}_trimmed.txt"
