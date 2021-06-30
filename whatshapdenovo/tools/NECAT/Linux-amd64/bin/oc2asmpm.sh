#!/bin/bash

fPrintUsage()
{
  echo "USAGE:"
  echo "$0 wrk_dir read_list output [map options]"
}


if [ $# -lt 3 ]; then
  fPrintUsage;
  exit 1;
fi

MKDB_WORKDER=oc2mkdb
PM_WORKER=oc2asmpm
WRK_DIR=$1
READ_LIST=$2
OUTPUT=$3
MAP_OPTIONS=""
ASM_OVLP_OPTIONS=${OSA_ASM_OVLP_OPTIONS}

ALL_FINISHED="${WRK_DIR}/all.finished"
if [ -f ${ALL_FINISHED} ]; then
    echo Job ${ALL_FINISHED} Has Been Finished. Exit Normally.
    exit 0;
fi

PAC_FINISHED="${WRK_DIR}/pac.finished"

i=1;
for opt in $*
do
  if [ $i -gt 3 ]; then
    MAP_OPTIONS="${MAP_OPTIONS} $opt"
  fi
  let i=i+1
done
MAP_OPTIONS="${MAP_OPTIONS} ${ASM_OVLP_OPTIONS}"

echo ${WRK_DIR}
echo ${READ_LIST}
echo ${OUTPUT}
echo ${MAP_OPTIONS}

if [ ! -f ${PAC_FINISHED} ]; then
	MKDB_CMD="${MKDB_WORKDER} ${WRK_DIR} ${READ_LIST}"
	oc2cmd.sh ${MKDB_CMD}
	if [ $? -ne 0 ]; then
  		exit 1;
	fi
	touch ${PAC_FINISHED}
fi

RINFO="${WRK_DIR}/reads_info.txt"
if [ ! -f ${RINFO} ]; then
  echo File ${RINFO} does not exist!
  exit 1;
fi

NVOLS=`cat ${RINFO} | awk '{print $1}' `
echo Number of volumes: $NVOLS

for((i=0;i<${NVOLS};i=i+1))
do
  RESULT="${WRK_DIR}/pm_result_$i"
  PM_FINISHED="${RESULT}.finished"
  if [ ! -f ${PM_FINISHED} ]; then
  	PM_CMD="${PM_WORKER} ${MAP_OPTIONS} ${WRK_DIR} $i ${RESULT}"
  	oc2cmd.sh ${PM_CMD}
  	if [ $? -ne 0 ]; then
    		exit 1;
  	fi
  	touch ${PM_FINISHED}
  fi
done

echo "mapping finish"

if [ -f ${OUTPUT} ]; then
	rm -f ${OUTPUT}
fi

for((i=0;i<${NVOLS};i=i+1))
do
  RESULT="${WRK_DIR}/pm_result_$i"
  cat ${RESULT} >> ${OUTPUT}
done

touch ${ALL_FINISHED}
