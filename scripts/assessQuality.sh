#!/bin/sh

set -eo pipefail

#   Check dependencies
declare -a DEPENDENCIES=(fastqc parallel unzip)
for dep in ${DEPENDENCIES[@]}; do $(command -v ${dep} > /dev/null 2> /dev/null) || echo ("Cannot find ${dep}" >&2; exit 1); done

#   Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) ... \n\
Where:      ... \n\
" >&2
    exit 1
}

export -f Usage

# A function to make specific outdirectories
function makeOutDirs() {
    local sample="$1"
    local outBase="$2"
    local sampleName="$(basename ${sample} .$(echo ${sample} | cut -f 2- -d '.'))"
    local outDir="${outBase}/${sampleName}"
    (set -x; mkdir -p "${outDir}")
    echo "${outDir}"
}

#   Export the function
export -f makeOutDirs

#   A function to summarize the results
function summarizeQC() {
    local zipFile="$1"
    local outDir="$2"
    #   Get some information
    local zipDir="${zipFile/.zip}"
    local sampleName="$(basename ${zipFile} _fastqc.zip)"
    #   Unzip our results
    (set -x; unzip -q "${zipFile}" -d "$(dirname ${zipFile})")
    #   Create our category files
    local -a notPass=($(grep -v 'PASS' "${zipDir}/summary.txt" | cut -f 1 | sort -u))
    for category in "${notPass[@]}"
    do
        (set -x; grep "${category}" "${zipDir}/summary.txt" > "${outDir}/${sampleName}_${category}.txt")
    done
    #   Clean up
    (set -x; rm -rf "${zipDir}")
}

#   Export the function
export -f summarizeQC

#   Make full table
function qcTable {
    local zipFile="$1"
    local tableName="$2"
    #   Get some information
    local zipDir="${zipFile/.zip}"
    local sampleName="$(basename ${zipFile} _fastqc.zip)"
    local dataFile="${zipDir}/fastqc_data.txt"
    #   Unzip results
    (set -x; unzip -q "${zipFile}" -d "$(dirname ${zipFile})")
    #   Make our table
    local -A results=()
    results=[total]=$(grep 'Total Sequences' "${dataFile}" | cut -f 2)
    results=[poor]=$(grep 'poor quality' "${dataFile}" | cut -f 2)
    results[length]=$(grep 'Sequence length' "${dataFile}" | cut -f 2)
    results[gc]=$(grep '%GC' "${dataFile}" | cut -f 2)
    (set -x; echo -e "${sampleName}\t${results[total]}\t${results[poor]}\t${results[length]}\t${results[gc]}" >> "${tableName}")
    #   Clean up
    (set -x; rm -rf "${zipDir}")
}

#   Export the function
export -f qcTable


#   Ensure we have the proper number of arguments
[[ "$#" -lt 3 ]] && Usage

#   Parse arguments
while [[ "$#" -gt 1 ]]; do
    case $1 in
        -s|--sample-list) # What is our list of samples?
            SAMPLE_LIST="$2"
            shift
            ;;
        -o|--outdir) # Where are we storing our results?
            OUT_DIR="${2}/Quality_Assesment"
            shift
            ;;
        -p|--project) # What do we call our results?
            PROJECT="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

#   Check arguments
[[ -z "${SAMPLE_LIST}" || -z "${OUT_DIR}" || -z "${PROJECT}" ]] && Usage
[[ -f "${SAMPLE_LIST}" ]] || (echo "Cannot find sample list ${SAMPLE_LIST}" >&2; exit 1)
for sample in $(<${SAMPLE_LIST}); do [[ -f "${sample}" ]] || (echo "Cannot find sample ${sample}" >&2; exit 1); done

#   Make our output directories
declare -a OUTDIRS=($(xargs -a "${SAMPLE_LIST}" -d '\n' -I {} sh -c "makeOutDirs {} ${OUT_DIR}"))
#   Run FastQC in parallel
(set -x; parallel --verbose --xapply "fastqc --outdir {2} {1}" :::: ${SAMPLE_LIST} ::: ${OUTDIRS[@]})
#   Collect our outputs
declare -a ZIP_FILES=($(find "${OUT_DIR}" -name "*.zip" | sort))
declare -a SORTED=($(for dir in ${OUTDIRS[@]}; do echo ${dir}; done | sort))
#   Summarize the results of FastQC
(set -x; parallel --verbose --xapply "summarizeQC {1} {2}" ::: ${ZIP_FILES[@]} ::: ${SORTED[@]})
#   Make project table
TABLE="${OUT_DIR}/${PROEJCT}_qc.txt"
echo -e "#File\tTotal Sequences\tPoor Quality\tSequence Length\t% GC" > "${TABLE}"
(set -x; parallel --verbose --keep-order "qcTable {} ${TABLE}" ::: ${ZIP_FILES[@]})
