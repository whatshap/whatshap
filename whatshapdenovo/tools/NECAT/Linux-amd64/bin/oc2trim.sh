#!/bin/bash

function usage() {
	echo "USAGE:"
	echo "$0 wrk_dir reads num_threads error_cutoff min_ovlp min_cov min_size output_reads output_ovlps"
}

if [ $# -ne 9 ]; then
	usage;
	exit 1;
fi

WRK_DIR=$1
INPUT_READS=$2
NTHREADS=$3
ERR_CUTOFF=$4
MIN_OVLP=$5
MIN_COV=$6
MIN_SIZE=$7
OUTPUT_READS=$8
OUTPUT_OVLPS=$9

TMP_TRIMMED_READS="${WRK_DIR}/tmp_trimReads.fasta"
TMP_TRIMMED_READS_OVLPS="${WRK_DIR}/tmp_trimReads.fasta.pm.m4"

if [ ! -d ${WRK_DIR} ]; then
	mkdir -p ${WRK_DIR}
fi

ALL_FINISHED="${WRK_DIR}/all.finished"
if [ -f ${ALL_FINISHED} ]; then
	echo "[$0 $@] is done. Exit normally..."
	exit 0
fi

### renumber reads
RENUM_READS_FINISHED="${WRK_DIR}/renum_reads.finished"
RENUM_READS="${WRK_DIR}/renum_reads.fasta"
RENUM_READS_CMD="oc2renumberSeqs ${INPUT_READS} ${RENUM_READS}"
if [ -f ${RENUM_READS_FINISHED} ]; then
	echo "[${RENUM_READS_CMD}] is done. Exir normally..."
else
	echo "[${RENUM_READS_CMD}] STARTS"
	${RENUM_READS_CMD}
	if [ $? -ne 0 ]; then
		echo "[${RENUM_READS_CMD}] FAIL"
		exit 1;
	fi
	echo "[${RENUM_READS_CMD}] ENDS"
	touch ${RENUM_READS_FINISHED}
fi

READ_LIST="${WRK_DIR}/read_list.txt"
echo ${RENUM_READS} > ${READ_LIST}

### all reads pairwise mapping

ALL_READS_PM_DIR="${WRK_DIR}/all_reads_ovlps"
ALL_READS_PM_RESULTS="${WRK_DIR}/all_reads_ovlps.m4.bin"
if [ ! -d ${ALL_READS_PM_DIR} ]; then
	mkdir -p ${ALL_READS_PM_DIR}
fi

JOB_FINISHED="${WRK_DIR}/all_reads_ovlps.finished"
ALL_READS_PM_CMD="oc2asmpm.sh ${ALL_READS_PM_DIR} ${READ_LIST} ${ALL_READS_PM_RESULTS} -t ${NTHREADS} -u 1 -k 13"
if [ -f ${JOB_FINISHED} ]; then
	echo "[${ALL_READS_PM_CMD}] is done. Exit normally..."
else
	echo "[${ALL_READS_PM_CMD}] STARTS"
	${ALL_READS_PM_CMD}
	if [ $? -ne 0 ]; then
		echo "[${ALL_READS_PM_CMD}] FAIL"
		exit 1
	fi
	echo "[${ALL_READS_PM_CMD}] ENDS"
	touch ${JOB_FINISHED}
fi

### partition m4
JOB_FINISHED="${WRK_DIR}/pm4.finished"
PM4_CMD="oc2pm4 ${ALL_READS_PM_DIR} ${ALL_READS_PM_RESULTS} ${ERR_CUTOFF} ${NTHREADS}"
if [ -f ${JOB_FINISHED} ]; then
	echo "[${PM4_CMD}] is done. Exit normally..."
else
	echo "[${PM4_CMD}] STARTS"
	${PM4_CMD}
	if [ $? -ne 0 ]; then
		echo "[${PM4_CMD}] FAIL"
		exit 1
	fi
	echo "[${PM4_CMD}] ENDS"
	touch ${JOB_FINISHED}
fi

### largest cover range
JOB_FINISHED="${WRK_DIR}/lcr.finished"
LCR_RESULTS="${WRK_DIR}/clipped_ranges.txt"
LCR_CMD="oc2lcr ${ALL_READS_PM_RESULTS} ${ALL_READS_PM_DIR} ${ERR_CUTOFF} ${MIN_OVLP} ${MIN_COV} ${MIN_SIZE} ${NTHREADS} ${LCR_RESULTS}"
if [ -f ${JOB_FINISHED} ]; then
	echo "[${LCR_CMD}] is done. Exit normally..."
else
	echo "[${LCR_CMD}] STARTS"
	${LCR_CMD}
	if [ $? -ne 0 ]; then
		echo "[${LCR_CMD}] FAIL"
		exit 1;
	fi
	echo "[${LCR_CMD}] ENDS"
	touch ${JOB_FINISHED}
fi

### extract trimmed reads
COMPLETE_READS="${WRK_DIR}/complete_reads.fasta"
UNCOMPLETE_READS="${WRK_DIR}/uncomplete_reads.fasta"
JOB_FINISHED="${WRK_DIR}/etr.finished"
ETR_CMD="oc2etr ${LCR_RESULTS} ${RENUM_READS} ${ALL_READS_PM_RESULTS} ${COMPLETE_READS} ${UNCOMPLETE_READS} ${TMP_TRIMMED_READS_OVLPS}"
if [ -f ${JOB_FINISHED} ]; then
	echo "[${ETR_CMD}] is done. Exit normally..."
else
	echo "[${ETR_CMD}] STARTS"
	${ETR_CMD}
	if [ $? -ne 0 ]; then
		echo $"[${ETR_CMD}] FAIL"
		exit 1;
	fi
	echo "[${ETR_CMD}] ENDS"
	touch ${JOB_FINISHED}
fi

### complete v.s. uncomplete
JOB_FINISHED="${WRK_DIR}/cu_ovlps.finished"
CU_OVLP_RESULTS="${WRK_DIR}/cu_ovlps.rm.m4"
CU_OVLP_CMD="oc2rm -k 13 -t ${NTHREADS} ${COMPLETE_READS} ${UNCOMPLETE_READS} ${CU_OVLP_RESULTS}"
if [ -f ${JOB_FINISHED} ]; then
	echo "[${CU_OVLP_CMD}] is done. Exit normally..."
else
	echo "[${CU_OVLP_CMD}] STARTS"
    ${CU_OVLP_CMD}
	if [ $? -ne 0 ]; then
		echo "[${CU_OVLP_CMD}] FAIL"
		exit 1;
	fi
	echo "[${CU_OVLP_CMD}] ENDS"
	touch ${JOB_FINISHED}
fi

### uncomplete v.s. uncomplte
JOB_FINISHED="${WRK_DIR}/uu_ovlps.finished"
UU_OVLP_RESULTS="${WRK_DIR}/uu_ovlps.pm.m4"
UU_OVLP_DIR="${WRK_DIR}/uu_ovlps"
UNCOMPLETE_READ_LIST="${WRK_DIR}/uncomplete_read_list.txt"
echo ${UNCOMPLETE_READS} > ${UNCOMPLETE_READ_LIST}
UU_OVLP_CMD="oc2asmpm.sh ${UU_OVLP_DIR} ${UNCOMPLETE_READ_LIST} ${UU_OVLP_RESULTS} -t ${NTHREADS} -k 13"
if [ -f ${JOB_FINISHED} ]; then
	echo "[${UU_OVLP_CMD}] is done. Exit normally..."
else
	echo "[${UU_OVLP_CMD}] STARTS"
	${UU_OVLP_CMD}
	if [ $? -ne 0 ]; then
		echo "[${UU_OVLP_CMD}] FAIL"
		exit 1;
	fi
	echo "[${UU_OVLP_CMD}] ENDS"
	touch ${JOB_FINISHED}
fi

cat ${COMPLETE_READS} > ${TMP_TRIMMED_READS}
cat ${UNCOMPLETE_READS} >> ${TMP_TRIMMED_READS}

cat ${CU_OVLP_RESULTS} >> ${TMP_TRIMMED_READS_OVLPS}
cat ${UU_OVLP_RESULTS} >> ${TMP_TRIMMED_READS_OVLPS}

### order results
JOB_FINISHED="${WRK_DIR}/order_results.finished"
ODR_RSLTS_CMD="oc2orderResults ${TMP_TRIMMED_READS} ${TMP_TRIMMED_READS_OVLPS} ${OUTPUT_READS} ${OUTPUT_OVLPS}"
if [ -f ${JOB_FINISHED} ]; then
	echo "[${ODR_RSLTS_CMD}] is done. Exit normally..."
else
	echo "[${ODR_RSLTS_CMD}] STARTS"
	${ODR_RSLTS_CMD}
	if [ $? -ne 0 ]; then
		echo "[${ODR_RSLTS_CMD}] FAIL"
		exit 1;
	fi
	echo "[${ODR_RSLTS_CMD}] ENDS"
	touch ${JOB_FINISHED}
fi

#if [ -f ${TMP_TRIMMED_READS} ]; then
#	rm ${TMP_TRIMMED_READS}
#fi

#if [ -f ${TMP_TRIMMED_READS_OVLPS} ]; then
#	rm ${TMP_TRIMMED_READS_OVLPS}
#fi

touch ${ALL_FINISHED}
