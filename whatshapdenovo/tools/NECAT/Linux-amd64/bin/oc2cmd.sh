#!/bin/bash

i=0;
for Q in $*;
do
    if [ $i -eq 0 ]; then
        CMD=$Q
    else
        CMD="${CMD} $Q"
    fi
    let i=i+1
done

if [ ! -n "${CMD}" ]; then
    echo Empty command!
    exit 1;
fi

echo "[$(date)] (${CMD}) STARTS"
${CMD}
if [ $? -ne 0 ]; then
    echo "Failed at Running (${CMD})"
    exit 1;
fi
echo "[$(date)] (${CMD}) FINISHES"
