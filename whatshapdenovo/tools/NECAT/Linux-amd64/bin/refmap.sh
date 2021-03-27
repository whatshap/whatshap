
WRK_DIR=/share/home/chuanlex/chenying/test/arab2/ara/5-align_contigs/rm_dir
NTHREADS=20
NODE_ID=0
NUM_NODES=1
READ_LSIT=/share/home/chuanlex/chenying/test/arab2/ara/1-consensus/raw_reads/raw_read_list.txt
REFERENCE=/share/home/chuanlex/chenying/test/arab2/ara/4-fsa/contigs.fasta
OUTPUT=/share/home/chuanlex/chenying/test/arab2/ara/5-align_contigs/rawread2ctg.m4a

### STEP 1: pack reads, one onde
MKDB_CMD="oc2mkdb ${WRK_DIR} ${READ_LIST}"
${MKDB_CMD}

### STEP 2: referecence mapping, multiple nodes

RM_CMD="oc2rm_worker -t ${NTHREADS} ${WRK_DIR} ${REFERENCE} ${OUTPUT} -mn ${NODE_ID} ${NUM_NODES}"

### STEP 2: reference mapping, one node
RM_CMD="oc2rm_worker -t ${NTHREADS} ${WRK_DIR} ${REFERENCE} ${OUTPUT}"

