#!/bin/bash

PARSE_ARG_SUCCESS=1
WRK_DIR=""
CNS_OPTIONS="${ONTCNS_CNS_OPTIONS}"
CAN_OPTIONS="${ONTCNS_CAN_OPTIONS}"
ONT_READ_LIST=""
READ_LIST=""
CANDIDATES=""
CNS_RESULTS=""
UNCNS_RESULTS=""
PACK_WORKER="oc2mkdb"
CAN_WORKER="oc2pm"
PCAN_WORKER="oc2pcan"
CNS_WORKER="oc2cns"
PACKED_DATA_DIR=""
NTHREADS=""

fPrintUsage()
{
    echo "USAGE:"
    echo "$0 -w wrk_dir -r results -u uncorrected_results -n ont_read_list -t threads -c cleanup"
}

fValidateOptions()
{
    if [ ! -n "${WRK_DIR}" ]; then
        echo "Working Folder (-w) Is Not Set!"
        PARSE_ARG_SUCCESS=0;
    fi
    if [ ! -n "${READ_LIST}" ]; then
        echo "Read List (-n) Is Not Set!"
        PARSE_ARG_SUCCESS=0;
    fi
    if [ ! -n "${CNS_RESULTS}" ]; then
        echo "Corrected Results Is Not Set!"
        PARSE_ARG_SUCCESS=0;
    fi
    if [ ! -n "${UNCNS_RESULTS}" ]; then
        echo "Uncorrected Results Is Not Set!"
        PARSE_ARG_SUCCESS=0;
    fi
    if [ ! -n "${CLEANUP}" ]; then
	echo "Cleanup Is Not Set!"
	PARSE_ARG_SUCESS=0;
    fi
    if [ ! -n "${NTHREADS}" ]; then
	echo "Number of Threads Is Not Set!"
	PARSE_ARG_SUCCESS=0;
    fi
}

while getopts "w:n:r:u:c:t:" OPTION;
do
    case ${OPTION} in
        w)
            WRK_DIR=${OPTARG}
            ;;
        n)
            ONT_READ_LIST=${OPTARG}
            ;;
        r)
            CNS_RESULTS=${OPTARG}
            ;;
        u)
            UNCNS_RESULTS=${OPTARG}
            ;;
	c)
	    CLEANUP=${OPTARG}
	    ;;
	t)
	    NTHREADS=${OPTARG}
	    ;;
        ?)
            PARSE_ARG_SUCCESS=0
            break;
            ;;
    esac
done

READ_LIST=${ONT_READ_LIST}

fValidateOptions;

if [ ${PARSE_ARG_SUCCESS} -eq 0 ]; then
    echo "Parse Argument Failed!"
    fPrintUsage;
    exit 1;
fi

TMP_CNS_RESULTS="${WRK_DIR}/tmp_cns.fasta"
TMP_UNCNS_RESULTS="${WRK_DIR}/tmp_raw.fasta"

PACK_JOB_FINISHED="${WRK_DIR}/pac.finished"
CAN_JOB_FINISHED="${WRK_DIR}/can.finished"
PCAN_JOB_FINISHED="${WRK_DIR}/pcan.finished"
CNS_JOB_FINISHED="${WRK_DIR}/cns.finished"
if [ -f ${CNS_JOB_FINISHED} ]; then
    echo "Job (${WRK_DIR}) Has Been Finished, Skip It."
    exit 0;
fi

PACKED_DATA_DIR="${WRK_DIR}/PackedData"
echo "Packed Data Dir: ${PACKED_DATA_DIR}"
if [ ! -d ${PACKED_DATA_DIR} ]; then
    mkdir -p ${PACKED_DATA_DIR}
fi
CANDIDATES="${WRK_DIR}/cns_candidates.txt"

echo "Working Dir: ${WRK_DIR}"
echo "Candidate Options: ${CAN_OPTIONS}"
echo "Consensus Options: ${CNS_OPTIONS}"
echo "Read List: ${READ_LIST}"

### pack data
if [ -f ${PACK_JOB_FINISHED} ]; then
    echo "Job [PackData] Has Been Finished, Skip It."
else
    pwd
    PACK_CMD="${PACK_WORKER} ${PACKED_DATA_DIR} ${READ_LIST}"
    oc2cmd.sh ${PACK_CMD}
    if [ $? -ne 0 ]; then
        exit 1;
    fi
    touch ${PACK_JOB_FINISHED}
fi

### find candidates
if [ -f ${CAN_JOB_FINISHED} ]; then
    echo "Job [Find Candidates] Has Been Finished, Skip It."
else
    CAN_CMD="${CAN_WORKER} ${CAN_OPTIONS} ${PACKED_DATA_DIR} ${CANDIDATES}"
    oc2cmd.sh ${CAN_CMD}
    if [ $? -ne 0 ]; then
        exit 1;
    fi
    touch ${CAN_JOB_FINISHED}
fi

### partition candidates
if [ -f ${PCAN_JOB_FINISHED} ]; then
    echo "Job [Partition Candidates] Has Been Finished, Skip It."
else
    PCAN_CMD="${PCAN_WORKER} -t ${NTHREADS} ${PACKED_DATA_DIR} ${CANDIDATES}"
    oc2cmd.sh ${PCAN_CMD}
    if [ $? -ne 0 ]; then
        exit 1;
    fi
    touch ${PCAN_JOB_FINISHED}
fi

### delete candidates after partition
if [ ${CLEANUP} -eq 1 ]; then
    rm -f ${CANDIDATES}
fi

### consensus
CNS_CMD="${CNS_WORKER} ${CNS_OPTIONS} ${PACKED_DATA_DIR} ${CANDIDATES} ${TMP_CNS_RESULTS} ${TMP_UNCNS_RESULTS}"
oc2cmd.sh ${CNS_CMD}
if [ $? -ne 0 ]; then
    exit 1;
fi

### reorder corrected results
REORDER_CMD="oc2ReorderCnsReads ${TMP_CNS_RESULTS} ${TMP_UNCNS_RESULTS} ${CNS_RESULTS} ${UNCNS_RESULTS}"
oc2cmd.sh ${REORDER_CMD}
if [ $? -ne 0 ]; then
	exit 1;
fi
touch ${CNS_JOB_FINISHED}

### delete candidate partitions after consensus
if [ ${CLEANUP} -eq 1 ]; then
    rm -f ${CANDIDATES}.p*
fi

### delete packed data dir
if [ ${CLEANUP} -eq 1 ]; then
    rm -rf ${PACKED_DATA_DIR}
fi

### delete tmp corrected results
if [ ${CLEANUP} -eq 1 ]; then
    rm -f ${TMP_CNS_RESULTS}
    rm -f ${TMP_UNCNS_RESULTS}
fi
