#!/bin/bash

set -eo pipefail

#   Usage message
function Usage() {
    echo -e "\
Usage: $(basename $0) -d|--directory directory [-p|--project project] \n\
Where:  -d|--directory is the STAR output directory \n\
        -p|--project is an optional project name \n\
            Defaults to the basename of the STAR output directory \n\
            eg. if -d|--directory is /home/user/STARoutput \n\
            -p|--project defaults to 'STARoutput' \n\
" >&2
    exit 1
}

#   Check if we have enough arguments
[[ "$#" -lt 1 ]] && Usage

#   Parse the arguments
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        -d|--directory)
            DIRECTORY="$2"
            shift
            ;;
        -p|--project)
            PROJECT="$2"
            shift
            ;;
        --no-wait)
            WAIT=false
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

#   Check the arguments
[[ -z "${DIRECTORY}" ]] && Usage
[[ -d "${DIRECTORY}" ]] || (echo "Cannot find STAR directory ${DIRECTORY}" >&2; exit 1)
[[ -z "${PROJECT}" ]] && PROJECT="$(basename ${DIRECTORY})"
[[ -z "${WAIT}" ]] && WAIT=true


declare -a SAMS=($(find "${DIRECTORY}" -name 'Aligned.out.sam'))
[[ "${#SAMS[@]}" -lt 1 ]] && (echo "No SAMs found" >&2; exit 1)

echo "Renaming ${#SAMS[@]} SAM files from 'Aligned.out.sam' to '${PROJECT}.out.sam'" >&2
if $(${WAIT}); then
    echo "Press ^c to cancel" >&2
    for i in {7..1}; do echo -en "\rContinuing in $i"; sleep 1; done
    echo
fi

for sam in ${SAMS[@]}; do (set -x; mv ${sam} ${sam/Aligned/${PROJECT}}); done
