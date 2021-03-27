#!/bin/bash

fPrintUsage()
{
	echo "USAGE:"
	echo "$0 wrk_dir reads num_threads error_cutoff min_ovlp min_cov min_read_size trim_reads"
}

if [ $# -ne 8 ]; then
	fPrintUsage
	exit 1
fi

WRK_DIR=$1
READS_IN=$2
NTHREADS=$3
ERROR_CUTOFF=$4
MIN_OVLP=$5
MIN_COV=$6
MIN_SIZE=$7
READS_OUT=$8
TRIM_OVLP_OPTIONS=${OSA_TRIM_OVLP_OPTIONS}

if [ ! -d ${WRK_DIR} ]; then
	mkdir -p ${WRK_DIR}
fi

ALL_FINISHED="${WRK_DIR}/all.finished"
if [ -f ${ALL_FINISHED} ]; then
    echo "triming bases has been finished. exit"
    exit 0;
fi

### renumber reads
RENUM_READS_FINISHED="${WRK_DIR}/renum_reads.finished"
RENUM_READS="${WRK_DIR}/renum_reads.fasta"
if [ ! -f ${RENUM_READS_FINISHED} ]; then
    RENUM_SEQ_CMD="oc2renumberSeqs ${READS_IN} ${RENUM_READS}"
    oc2cmd.sh ${RENUM_SEQ_CMD}
    if [ $? -ne 0 ]; then
	    eixt 1;
    fi
    touch ${RENUM_READS_FINISHED}
fi

READ_LIST="${WRK_DIR}/read_list.txt"
echo ${RENUM_READS} > ${READ_LIST}

### pairwise mapping

PM_WRK_DIR="${WRK_DIR}/pac_in"
PM_RESULT="${WRK_DIR}/pm.m4"
if [ ! -d ${PM_WRK_DIR} ]; then
	mkdir -p ${PM_WRK_DIR}
fi

PM_FINISHED="${WRK_DIR}/pm.finished"
if [ ! -f ${PM_FINISHED} ]; then
    PM_CMD="oc2asmpm.sh ${PM_WRK_DIR} ${READ_LIST} ${PM_RESULT} -t ${NTHREADS} ${TRIM_OVLP_OPTIONS}"
    oc2cmd.sh ${PM_CMD}
    if [ $? -ne 0 ]; then
	    exit 1
    fi
    touch ${PM_FINISHED}
fi

### partition m4
PM4_FINISHED="${WRK_DIR}/pm4.finised"
if [ ! -f ${PM4_FINISHED} ]; then
    PM4_CMD="oc2pm4 ${PM_WRK_DIR} ${PM_RESULT} ${ERROR_CUTOFF} ${NTHREADS}"
    oc2cmd.sh ${PM4_CMD}
    if [ $? -ne 0 ]; then
	    exit 1;
    fi
    touch ${PM4_FINISHED}
fi

### compute max coverage range
CLIPPED_FINISHED="${WRK_DIR}/clipped_ranges.finisehd"
CLIPPED_RANGES="${WRK_DIR}/clipped_ranges.txt"
if [ ! -f ${CLIPPED_FINISHED} ]; then
#    LCR_CMD="oc2lcr ${PM_RESULT} ${PM_WRK_DIR} ${ERROR_CUTOFF} ${MIN_OVLP} ${MIN_COV} ${MIN_SIZE} ${CLIPPED_RANGES} ${NTHREADS}"
    LCR_CMD="oc2lcr ${PM_RESULT} ${PM_WRK_DIR} ${ERROR_CUTOFF} ${MIN_OVLP} ${MIN_COV} ${MIN_SIZE} ${READS_OUT} ${NTHREADS}"
    oc2cmd.sh ${LCR_CMD}
    if [ $? -ne 0 ]; then
	    exit 1;
    fi
    touch ${CLIPPED_FINISHED}
fi

### trim reads
#TB_FINISHED="${WRK_DIR}/tb.finished"
#if [ ! -f ${TB_FINISHED} ]; then
#    TB_CMD="oc2trimBases ${RENUM_READS} ${CLIPPED_RANGES} ${READS_OUT}"
#    oc2cmd.sh ${TB_CMD}
#    if [ $? -ne 0 ]; then
#	    exit ;
#    fi
#    touch ${TB_FINISHED}
#fi

touch ${ALL_FINISHED}
