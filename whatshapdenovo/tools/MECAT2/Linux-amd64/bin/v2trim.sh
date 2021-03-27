#!/bin/bash

fPrintUsage()
{
	echo "USAGE:"
	echo "$0 wrk_dir reads cpu_threads"
}

if [ $# -ne 3 ]; then
	fPrintUsage;
	exit 1;
fi

WRK_DIR=$1
INPUT=$2
NTHREADS=$3

if [ ! -d ${WRK_DIR} ]; then
	mkdir -p ${WRK_DIR}
fi

ALL_FINISHED=${WRK_DIR}/all.finished

if [ -f ${ALL_FINISHED} ]; then
	echo "All Jobs Have Been Finished. Exit Normally."
	exit 0;
fi

TRIM_PM_DIR="${WRK_DIR}/trim_pm_dir"
if [ ! -d ${TRIM_PM_DIR} ]; then
	mkdir -p ${TRIM_PM_DIR}
fi

TRIM_PM_RESULT="${WRK_DIR}/trim_pm.m4"
TRIM_PM_CMD="v2asmpm.sh ${TRIM_PM_DIR}  ${INPUT} ${NTHREADS} ${TRIM_PM_RESULT} -B"
${TRIM_PM_CMD}
if [ $? -ne 0 ]; then
	echo "Failed at running ${TRIM_PM_CMD}"
	exit 1;
fi

PM4_CMD="v2pm4 ${TRIM_PM_DIR} ${TRIM_PM_RESULT} 0.09 ${NTHREADS}"
${PM4_CMD}
if [ $? -ne 0 ]; then
	echo "Failed at running ${PM4_CMD}"
	exit 1;
fi

LCR_RESULT="${WRK_DIR}/lcr.txt"
LCR_CMD="v2lcr ${TRIM_PM_RESULT} ${TRIM_PM_DIR} 0.09 1 1 500 ${LCR_RESULT} ${NTHREADS}"
${LCR_CMD}
if [ $? -ne 0 ]; then
	echo "Failed at running ${LCR_CMD}"
	exit 1;
fi

SR_RESULT="${WRK_DIR}/sr.txt"
SR_CMD="v2sr ${TRIM_PM_RESULT} ${TRIM_PM_DIR} ${LCR_RESULT} 500 ${SR_RESULT} ${NTHREADS}"
${SR_CMD}
if [ $? -ne 0 ]; then
	echo "Failed at running ${SR_CMD}"
	exit 1;
fi

TRIM_READS="${WRK_DIR}/trimReads.fasta"
TB_CMD="v2tb ${INPUT} ${SR_RESULT} ${TRIM_READS}"
${TB_CMD}
if [ $? -ne 0 ]; then
	echo "Failed at running ${TB_CMD}"
	exit 1;
fi

ASM_PM_DIR="${WRK_DIR}/asm_pm_dir"
ASM_PM_RESULT="${WRK_DIR}/asm_pm.m4"
ASM_PM_CMD="v2asmpm.sh ${ASM_PM_DIR} ${TRIM_READS} ${NTHREADS} ${ASM_PM_RESULT}"
${ASM_PM_CMD}
if [ $? -ne 0 ]; then
	echo "Failed at running ${ASM_PM_CMD}"
	exit 1;
fi

touch ${ALL_FINISHED}
exit 0
