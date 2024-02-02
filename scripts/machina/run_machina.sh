#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh

if [[ $# -eq 0 ]] ; then
    echo "Usage: run_machina.sh --edges <edges_file> --labels <labels_file> --colors <colors_file> --primary-tissue <tissue> --outdir <dir_path>"
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--edges) EDGES="$2"; shift ;;
        -l|--labels) LABELS="$2"; shift ;;
        -c|--colors) COLORS="$2"; shift ;;
        -p|--primary-tissue) PTISSUE="$2"; shift ;;
        -o|--outdir) OUTDIR=$2; shift ;;

    *) echo "Usage: run_machina.sh --edges <edges_file> --labels <labels_file> --colors <colors_file> --primary-tissue <tissue> --outdir <dir_path>" ; exit 1 ;;
    esac
    shift
done

if ! command -v pmh_tr &> /dev/null
then
    echo "pmh_tr from the MACHINA installation could not be found. Exiting!"
    exit
fi

if [ -z "${THREADS}" ]
then
    THREADS=1
fi

if [ -z "${EDGES}" ]
then
    echo "Missing --edges parameter. Exiting!"
    exit
fi
if [ -z "${LABELS}" ]
then
     echo "Missing --labels parameter. Exiting!"
    exit
fi
if [ -z "${COLORS}" ]
then
     echo "Missing --colors parameter. Exiting!"
    exit
fi
if [ -z "${PTISSUE}" ]
then
    echo "Missing --primary-tissue parameter. Exiting!"
    exit
fi
if [ -z "${OUTDIR}" ]
then
    echo "Missing --outdir parameter. Exiting!"
    exit
fi

MACHINA="pmh_sankoff"

if [[ -d $OUTDIR ]]; then
    :
else
    mkdir ${OUTDIR}
fi

${MACHINA} -p ${PTISSUE} -c ${COLORS} -o ${OUTDIR} ${EDGES} ${LABELS} &> ${OUTDIR}/machina_results.txt