getOsMachineType() {
    OSTYPE=`uname`
    MACHINETYPE=`uname -m`

    if [ ${MACHINETYPE} == "x86_64" ]; then
        MACHINETYPE="amd64"
    fi

    if [ ${MACHINETYPE} == "Power Macintosh" ]; then	
        MACHINETYPE="ppc"
    fi

    if [ ${OSTYPE} == "SunOS" ]; then
        MACHINETYPE=`uname -p`
        if [ ${MACHINETYPE} == "sparc" ]; then
            if [ `/usr/bin/isainfo -b` == "64" ]; then
                MACHINETYPE=sparc64
            else
                MACHINETYPE=sparc32
            fi
        fi
    fi
    echo ${OSTYPE}-${MACHINETYPE}
}

checkReturn() {
    if [ $? != 0 ]; then
        echo $1
        exit 1;
    fi
}

CONFIG_FILE=$1
if [ "$2" =  "" ]; then
  steps=("0" "1" "2" "3")
else
  steps=(${2//,/ })
fi

basepath=$(cd `dirname $0`; pwd)
echo $basepath
PATH=$basepath:$basepath:$PATH

### parse arguments
while read line ; do
    eval "${line}"
done < ${CONFIG_FILE}

for step in ${steps[@]}
do
    if [[ ${step} == "0" ]];then
        ontsa.sh ${CONFIG_FILE}
        checkReturn "Failed to run ontsa.sh" 
    fi

    dir_fsa=${PROJECT}/4-fsa
    mkdir ${dir_fsa} -p
    if [[ ${step} == "1" ]];then
        fsa_ol_filter  ${PROJECT}/3-assembly/pm.m4 ${dir_fsa}/filter.m4 --thread_size=${THREADS} --genome_size=${GENOME_SIZE} --output_directory=${dir_fsa} ${FSA_OL_FILTER_OPTIONS}
        checkReturn "Failed to run fsa_ol_filter"
        fsa_assemble ${dir_fsa}/filter.m4  --read_file=${PROJECT}/trimReads.fasta --output_directory=${dir_fsa}  --thread_size=${THREADS} ${FSA_ASSEMBLE_OPTIONS}
        checkReturn "Failed to run fsa_assemble"
    
        echo "The contig file: ${PROJECT}/4-fsa/contigs.fasta"
        fsa_rd_stat ${PROJECT}/4-fsa/contigs.fasta --thread_size=${THREADS}
        checkReturn "Failed to run fsa_rd_stat contifs"
    fi

    dir_align_contigs=${PROJECT}/5-align_contigs

    mkdir ${dir_align_contigs} -p
        if [[ ${step} == "2" ]];then
        lineno=0
        while read line; do
            rawread_file=${line}
            oc2rm -t ${THREADS} ${rawread_file} ${PROJECT}/4-fsa/contigs.fasta  ${dir_align_contigs}/rawread2ctg.m4a.$lineno
            checkReturn "Failed to run oc2rm rawreads contigs"
            lineno=$((lineno + 1))
        #done < ${ONT_READ_LIST}
        done < ${PROJECT}/1-consensus/raw_reads/raw_read_list.txt # The garbage reads have been removed
        cat  ${dir_align_contigs}/rawread2ctg.m4a.* >  ${dir_align_contigs}/rawread2ctg.m4a
        oc2ctgpm.sh ${dir_align_contigs}/temp  ${THREADS} ${dir_fsa}/contigs.fasta ${dir_align_contigs}/ctg2ctg.m4a
        checkReturn "Failed to run oc2ctgpm.sh contigs"
    fi


    dir_bridge_contigs=${PROJECT}/6-bridge_contigs
    mkdir ${dir_bridge_contigs} -p
    if [[  ${step} == "3" ]];then
        fsa_ctg_bridge  ${ONT_READ_LIST} ${dir_fsa}/contigs.fasta ${dir_align_contigs}/rawread2ctg.m4a ${dir_bridge_contigs}/bridged_contigs.fasta --ctg2ctg_file=${dir_align_contigs}/ctg2ctg.m4a --thread_size=${THREADS} --output_directory=${dir_bridge_contigs} ${FSA_CTG_BRIDGE_OPTIONS}        
        checkReturn "Failed to run fsa_ctg_bridge"

        echo "The contig file: ${dir_fsa}/contigs.fasta"
        fsa_rd_stat ${dir_fsa}/contigs.fasta --thread_size=${THREADS}
        checkReturn "Failed to run fsa_rd_stat contigs"

        echo "The briged contig file: ${dir_bridge_contigs}/bridged_contigs.fasta"
        fsa_rd_stat ${dir_bridge_contigs}/bridged_contigs.fasta --thread_size=${THREADS}
        checkReturn "Failed to run fsa_rd_stat bridged-contigs"
    fi

done


