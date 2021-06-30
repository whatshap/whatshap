#!/bin/bash

fPrintUsage()
{
    echo "USAGE:"
    echo "$0 WRK_DIR NUM_THREADS CONTIGS OUTPUT"
}

if [ $# -ne 4 ]; then 
    fPrintUsage;
    exit 1;
fi

WRK_DIR=$1
NTHREADS=$2
CONTIGS=$3
OUTPUT=$4

if [ ! -d ${WRK_DIR} ]; then
    mkdir -p ${WRK_DIR}
fi

CTG_READS="${WRK_DIR}/ctg_reads.fasta"

SPLIT_CTG_CMD="oc2SplitCtgs ${CONTIGS} ${CTG_READS}"
oc2cmd.sh ${SPLIT_CTG_CMD}
if [ $? -ne 0 ]; then
    exit 1;
fi

CTG_READ_LIST="${WRK_DIR}/ctg_read_list.txt"
echo ${CTG_READS} > ${CTG_READ_LIST}

PAC_CMD="oc2mkdb ${WRK_DIR} ${CTG_READ_LIST}"
oc2cmd.sh ${PAC_CMD}
if [ $? -ne 0 ]; then
    exit 1;
fi

TMP_CANDIDATES="${WRK_DIR}/tmp_candidates.txt"
CAN_CMD="oc2pm -t ${NTHREADS} -j 0 ${WRK_DIR} ${TMP_CANDIDATES}"
oc2cmd.sh ${CAN_CMD}
if [ $? -ne 0 ]; then
    exit 1;
fi

CANDIDATES="${WRK_DIR}/candidates.txt"
FIX_CAN_CMD="oc2FixCanInfo ${CTG_READS} ${TMP_CANDIDATES} ${CANDIDATES}"
oc2cmd.sh ${FIX_CAN_CMD}
if [ $? -ne 0 ]; then
    exit 1;
fi

PM_CMD="oc2ctgpm -t ${NTHREADS} ${CONTIGS} ${CANDIDATES} ${OUTPUT}"
oc2cmd.sh ${PM_CMD}
if [ $? -ne 0 ]; then
    exit 1;
fi
