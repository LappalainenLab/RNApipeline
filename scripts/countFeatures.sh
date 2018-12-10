#!/bin/bash

set -eo pipefail

$(command -v featureCounts > /dev/null 2> /dev/null) || (echo "Please install featureCounts..." >&2; exit 1)

function Usage() {
    echo -e "\
Usage: $(basename $0) -l|--sample-list -a|--annotation -o|--output
" >&2
    exit 1
}

export -f Usage


while [[ $# -gt 1 ]]; do
    KEY=$1
    case $KEY in
        -l|--sample-list)
            SAMPLE_LIST=$2
            shift
            ;;
        -a|--annotation)
            ANNOTATION=$2
            shift
            ;;
        -o|--output)
            OUTPUT=$2
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

[[ -z "${SAMPLE_LIST}" || -z "${ANNOTATION}" || -z "${OUTPUT}" ]] && Usage
[[ -f "${ANNOTATION}" ]] || (echo "Cannot find annotation file ${ANNOTATION}" >&2; exit 1)
[[ -f "${SAMPLE_LIST}" ]] || (echo "Cannot find sample list ${SAMPLE_LIST}" >&2; exit 1)

(set -x; mkdir -p "$(dirname ${OUTPUT})")

declare -a BAMS=($(<${SAMPLE_LIST}))
[[ "${#BAMS[@]}" -lt 1 ]] && (echo "No samples found in sample list ${SAMPLE_LIST}" >&2; exit 1)
for sample in ${BAMS[@]}; do
    [[ -f "${sample}" ]] || (echo "Cannot find sample ${sample}" >&2; exit 1)
done

echo "Counting features in ${#BAMS[@]} BAM files" >&2

THREADS=$(cat /proc/cpuinfo | grep processor | wc -l)

(set -x; featureCounts --primary --verbose -T "${THREADS}" -a "${ANNOTATION}" -o "${OUTPUT}" ${BAMS[@]})
