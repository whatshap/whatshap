#!/bin/bash

fPrintUsage()
{
	echo "USAGE:"
	echo "$0 cns_dir config_file"
}

if [ $# -ne 2 ]; then
	fPrintUsage;
	exit 1;
fi

CNS_DIR=$1
CONFIG_FILE=$2

### parse arguments
while read LINE ; do
    eval "${LINE}"
done < ${CONFIG_FILE}

echo ${THREADS}
echo ${ONT_READ_LIST}
echo ${MIN_READ_LENGTH}
echo ${OVLP_FAST_OPTIONS}
echo ${OVLP_SENSITIVE_OPTIONS}
echo ${CNS_FAST_OPTIONS}
echo ${CNS_SENSITIVE_OPTIONS}
echo ${NUM_ITER}
echo ${CLEANUP}

READ_LIST="-n ${ONT_READ_LIST}"

if [ ! -d ${CNS_DIR} ]; then
    mkdir -p ${CNS_DIR}
fi

PPRR_DIR="${CNS_DIR}/raw_reads"
if [ ! -d ${PPRR_DIR} ]; then
	mkdir -p ${PPRR_DIR}
fi
TMP_READ_LIST="${PPRR_DIR}/raw_read_list.txt"
PPRR_FINISHED="${PPRR_DIR}/pprr.finished"
if [ -f ${PPRR_FINISHED} ]; then
	echo "Job PPRR has been finished, skip it."
else
	PPRR_CMD="oc2pprr ${ONT_READ_LIST} ${PPRR_DIR} ${TMP_READ_LIST} ${MIN_READ_LENGTH}"
	oc2cmd.sh ${PPRR_CMD}
	if [ $? -ne 0 ]; then
    		exit 1;
	fi
	touch ${PPRR_FINISHED}
fi
READ_LIST="-n ${TMP_READ_LIST}"


for ((i=1;i<=${NUM_ITER};i++))
do
    DIR_NAME="cns_iter${i}"
    WRK_DIR="${CNS_DIR}/${DIR_NAME}"
    CNS_READS="${WRK_DIR}/cns.fasta"
    UNCNS_READS="${WRK_DIR}/raw.fasta"

    if [ ! -d ${WRK_DIR} ]; then
        mkdir -p ${WRK_DIR}
    fi

    if [ $i -eq 1 ]; then
        ITER_CAN_OPTIONS="-t ${THREADS} ${OVLP_SENSITIVE_OPTIONS}"
        ITER_CNS_OPTIONS="-t ${THREADS} ${CNS_SENSITIVE_OPTIONS} -r 0"
    else
        ITER_CAN_OPTIONS="-t ${THREADS} ${OVLP_FAST_OPTIONS}"
        ITER_CNS_OPTIONS="-t ${THREADS} ${CNS_FAST_OPTIONS} -r 1"
    fi

    if [ $i -eq ${NUM_ITER} ]; then
        ITER_CNS_OPTIONS="${ITER_CNS_OPTIONS} -f 0"
    else
        ITER_CNS_OPTIONS="${ITER_CNS_OPTIONS} -f 1"
    fi
    export ONTCNS_CAN_OPTIONS=${ITER_CAN_OPTIONS}
    export ONTCNS_CNS_OPTIONS=${ITER_CNS_OPTIONS}

    CMD="oc2cns_iter.sh -t ${THREADS} -w ${WRK_DIR} -r ${CNS_READS} -u ${UNCNS_READS} ${READ_LIST} -c ${CLEANUP}"
    oc2cmd.sh ${CMD}
    if [ $? -ne 0 ]; then
        exit 1;
    fi

    if [ $i -ne ${NUM_ITER} ]; then
        NEXT_READ_LIST="${WRK_DIR}/ReadList.txt"
        echo "${CNS_READS}" > "${NEXT_READ_LIST}"
        echo "${UNCNS_READS}" >> "${NEXT_READ_LIST}"
        READ_LIST="-n ${NEXT_READ_LIST}"
    fi
done
