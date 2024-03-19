#!/bin/bash

fPrintUsage()
{
	echo "USAGE:"
	echo "$0 wrk_dir input cpu_threads output [options]"
}

if [ $# -lt 4 ]; then
	fPrintUsage;
	exit 1;
fi

MAP_OPTIONS=""

i=1
MAP_OPTIONS=""
for OPT in $*
do
	if [ $i -gt 4 ]; then 
		MAP_OPTIONS="${MAP_OPTIONS} $OPT"
	fi
let i=i+1
done

WRK_DIR=$1
INPUT=$2
NTHREADS=$3
OUTPUT=$4

MKVOL="v2mkvol"
MKVOL_JOB_FINISHED="${WRK_DIR}/mk_vol.finished"
PM="v2asmpm"

ALL_FINISHED="${WRK_DIR}/all.finished"

if [ ! -d ${WRK_DIR} ]; then
	mkdir -p ${WRK_DIR}
fi

if [ -f ${ALL_FINISHED} ]; then
	echo "All Jobs Have Been Finished. Exit Normally."
	exit 0;
fi

if [ -f ${MKVOL_JOB_FINISHED} ]; then
	echo "JOB mk_vol.finished has been done. Skip it."
else
	MK_VOL_CMD="${MKVOL} ${WRK_DIR} ${INPUT}"
	${MK_VOL_CMD}
	if [ $? -ne 0 ]; then
		echo "Fail at running (${MK_VOL_CMD})"
		exit 1;
	fi
	touch ${MKVOL_JOB_FINISHED}
fi

NVOL_FILE="${WRK_DIR}/num_volumes.txt"
if [ ! -f ${NVOL_FILE} ]; then
	echo "Fail to retrieve number of vlumes, file ${NVOL_FILE} not found!"
	exit 1;
fi
NVOL=`cat ${NVOL_FILE} | awk '{print $1}'`
echo "Number of volumes: ${NVOL}"

if [ -f ${OUTPUT} ]; then
	rm -f ${OUTPUT}
fi

for((i=1;i<=${NVOL};i=i+1))
do
	VOL_PM_FINISHED="${WRK_DIR}/pm_${i}.finished"
	echo ${VOL_PM_FINISHED}
	if [ -f ${VOL_PM_FINISHED} ]; then
		echo "Job ${VOL_PM_FINISHED} has been done. Skipt it."
	else
		VOL_PM_CMD="${PM} -P${WRK_DIR} -T${NTHREADS} -S${i} -E${NVOL} ${MAP_OPTIONS}"
		echo "Running ${VOL_PM_CMD}"
		${VOL_PM_CMD}
		if [ $? -ne 0 ]; then
			echo "Fail at running (${VOL_PM_CMD})"
			exit 1;
		fi
		cat ${WRK_DIR}/${i}_*.r >> ${OUTPUT}
		touch ${VOL_PM_FINISHED}
	fi
done

for((i=1;i<=${NVOL};i=i+1))
do
	rm -f ${WRK_DIR}/${i}_*.r
done

touch ${ALL_FINISHED}
