#!/bin/bash

fPrintUsage()
{
	echo "USAGE:"
	echo "$0 config_file"
}

if [ $# -ne 1 ]; then
        fPrintUsage;
        exit 1;
fi

CONFIG_FILE=$1

### parse arguments
while read LINE ; do
    eval "${LINE}"
done < ${CONFIG_FILE}

echo ${PROJECT}
echo ${THREADS}
echo ${NUM_ITER}
echo ${GENOME_SIZE}

PROJECT_DIR="$(pwd)/${PROJECT}"

if [ ! -d ${PROJECT_DIR} ]; then
        mkdir -p ${PROJECT_DIR};
fi

### consensus

CNS_WRK_DIR="${PROJECT_DIR}/1-consensus"
if [ ! -d ${CNS_WRK_DIR} ]; then
        mkdir -p ${CNS_WRK_DIR};
fi
CNS_CMD="oc2cns.sh ${CNS_WRK_DIR} ${CONFIG_FILE}"
oc2cmd.sh ${CNS_CMD}
if [ $? -ne 0 ]; then
        exit 1;
fi

### extract cns reads
if [ "${CNS_WRK_DIR}/cns_iter${NUM_ITER}/cns.fasta" -nt "${CNS_WRK_DIR}/cns_final.fasta" ]; then
    echo $CNS_OUTPUT_COVERAGE
    echo $GENOME_SIZE
    if [ ! $GENOME_SIZE = "" ] && [ ! $CNS_OUTPUT_COVERAGE = "" ]; then
        CMD="oc2elr ${CNS_WRK_DIR}/cns_iter${NUM_ITER}/cns.fasta $GENOME_SIZE $CNS_OUTPUT_COVERAGE ${CNS_WRK_DIR}/cns_final.fasta"
    else
        CMD="mv ${CNS_WRK_DIR}/cns_iter${NUM_ITER}/cns.fasta ${CNS_WRK_DIR}/cns_final.fasta"
    fi
     
    oc2cmd.sh ${CMD}
    if [ $? -ne 0 ]; then
        exit 1;
    fi
fi
### trim bases

export OSA_TRIM_OVLP_OPTIONS=${TRIM_OVLP_OPTIONS}
TB_WRK_DIR="${PROJECT_DIR}/2-trim_bases"
TB_READS="${PROJECT_DIR}/trimReads.fasta"
#CNS_READS="${CNS_WRK_DIR}/cns_iter${NUM_ITER}/cns.fasta"
CNS_READS="${CNS_WRK_DIR}/cns_final.fasta"
TB_CMD="oc2trim_bases.sh ${TB_WRK_DIR} ${CNS_READS} ${THREADS} 0.09 1 1 500 ${TB_READS}"
oc2cmd.sh ${TB_CMD}
if [ $? -ne 0 ]; then
	exit 1;
fi

### assembly

export OSA_ASM_OVLP_OPTIONS=${ASM_OVLP_OPTIONS}
ASM_WRK_DIR="${PROJECT_DIR}/3-assembly"
ASM_READ_LIST="${ASM_WRK_DIR}/asm_read_list.txt"
ASM_READS="${ASM_WRK_DIR}/asm_reads.fasta"
ASM_PM_DIR="${ASM_WRK_DIR}/pm_dir"
if [ ! -d ${ASM_WRK_DIR} ]; then
	mkdir -p ${ASM_WRK_DIR}
fi
echo ${TB_READS} > ${ASM_READ_LIST}
ASM_PM="${ASM_WRK_DIR}/pm.m4"
ASM_PM_CMD="oc2asmpm.sh ${ASM_PM_DIR} ${ASM_READ_LIST} ${ASM_PM} -t ${THREADS}"
oc2cmd.sh ${ASM_PM_CMD}
if [ $? -ne 0 ]; then
	exit 1;
fi
