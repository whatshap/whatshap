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

checkRetCode() {
    local retCode msg
    retCode=$1
    msg=$2

    if [ $retCode != 0 ]; then
        printf "retCode=%d, %s" $retCode, $msg
        exit 1;
    fi

}

delSavedRetCode() {
    local finished
    finished=$1
    rm $finished -f
}

getSavedRetCode() {
    local finished
    finished=$1
    if [ -e $finished ]; then
        cat $finished
    else
        echo 255
    fi
}

checkAndSaveRetCode() {
    local retCode finished msg
    retCode=$1
    finished=$2
    msg=$3

    echo $retCode > $finished
    if [ $retCode != 0 ]; then
        printf "retCode=%d, %s\n" $retCode "$msg"
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
        checkRetCode $?  "Failed to run ontsa.sh" 
    fi

    dir_fsa=${PROJECT}/4-fsa
    mkdir ${dir_fsa} -p
    if [[ ${step} == "1" && 
          ( ${PROJECT}/3-assembly/pm.m4 -nt  ${dir_fsa}/contigs.fasta  || 
            `getSavedRetCode "${dir_fsa}/finished"` != 0 ) ]];then
        
        delSavedRetCode "${dir_fsa}/finished"    

        fsa_ol_filter  ${PROJECT}/3-assembly/pm.m4 ${dir_fsa}/filter.m4 --thread_size=${THREADS} --genome_size=${GENOME_SIZE} --output_directory=${dir_fsa} ${FSA_OL_FILTER_OPTIONS}
        checkAndSaveRetCode $?  "${dir_fsa}/finished" "Failed to run fsa_ol_filter"
        fsa_assemble ${dir_fsa}/filter.m4  --read_file=${PROJECT}/trimReads.fasta --output_directory=${dir_fsa}  --thread_size=${THREADS} ${FSA_ASSEMBLE_OPTIONS}
        checkAndSaveRetCode $?  "${dir_fsa}/finished" "Failed to run fsa_assemble"

    else
        if [[ ${step} == "1" ]]; then
            echo "Outputs is newer, so skip assembling step"
        fi
    fi

    dir_align_contigs=${PROJECT}/5-align_contigs

    filtered_rawreads_list=${PROJECT}/1-consensus/raw_reads/raw_read_list.txt
    mkdir ${dir_align_contigs} -p
    if [[ ${step} == "2"  && 
          ( ${dir_fsa}/contigs.fasta -nt  ${dir_align_contigs}/rawread2ctg.m4a ||
            ${dir_fsa}/contigs.fasta -nt  ${dir_align_contigs}/ctg2ctg.m4a || 
            `getSavedRetCode "${dir_align_contigs}/finished"` != 0 ) ]]; then

        delSavedRetCode "${dir_align_contigs}/finished"    

        lineno=0
        while read line; do
            rawread_file=${line}
            if [[ ${dir_fsa}/contigs.fasta -nt  ${dir_align_contigs}/rawread2ctg.m4a.$lineno ||
                  `getSavedRetCode "${dir_align_contigs}/finished.r2c.$lineno"` != 0 ]]; then
                  
                delSavedRetCode "${dir_align_contigs}/finished.r2c.$lineno"    
                oc2rm -t ${THREADS} ${rawread_file} ${dir_fsa}/contigs.fasta  ${dir_align_contigs}/rawread2ctg.m4a.$lineno
                checkAndSaveRetCode $? "${dir_align_contigs}/finished.r2c.$lineno" "Failed to run oc2rm rawreads contigs"
            fi
            lineno=$((lineno + 1))

        done < ${filtered_rawreads_list} # The garbage reads have been removed

        if [[ ${dir_fsa}/contigs.fasta -nt  ${dir_align_contigs}/rawread2ctg.m4a ||
              `getSavedRetCode "${dir_align_contigs}/finished.r2c"` != 0 ]]; then
              
            delSavedRetCode "${dir_align_contigs}/finished.r2c"   
            cat  ${dir_align_contigs}/rawread2ctg.m4a.* >  ${dir_align_contigs}/rawread2ctg.m4a
            checkAndSaveRetCode $? "${dir_align_contigs}/finished.r2c" "Failed to run cat rawread2ctg.m4a.* > rawread2ctg.m4a"

            if [[ $CLEANUP == 1 ]]; then
                rm -f ${dir_align_contigs}/rawread2ctg.m4a.*
                rm -f ${dir_align_contigs}/finished.r2c.*
            fi
        fi
        
        if [[ ${dir_fsa}/contigs.fasta -nt  ${dir_align_contigs}/ctg2ctg.m4a ||
              `getSavedRetCode "${dir_align_contigs}/finished.c2c"` != 0 ]]; then
                
            delSavedRetCode "${dir_align_contigs}/finished.c2c"    
            oc2ctgpm.sh ${dir_align_contigs}/temp  ${THREADS} ${dir_fsa}/contigs.fasta ${dir_align_contigs}/ctg2ctg.m4a
            checkAndSaveRetCode $? "${dir_align_contigs}/finished.c2c" "Failed to run oc2ctgpm.sh contigs"
        fi
        
        checkAndSaveRetCode 0 "${dir_align_contigs}/finished" "Failed to run aligning contigs"

    else
        if [[ ${step} == "2" ]]; then
            echo "Outputs is newer, so skip contig-aligning step"
        fi
    fi
 

    dir_bridge_contigs=${PROJECT}/6-bridge_contigs
    mkdir ${dir_bridge_contigs} -p
    if [[  ${step} == "3" &&  
         ( ${dir_align_contigs}/rawread2ctg.m4a -nt ${dir_bridge_contigs}/bridged_contigs.fasta ||
           `getSavedRetCode "${dir_bridge_contigs}/finished"` != 0 ) ]];then

       delSavedRetCode "${dir_bridge_contigs}/finished"    

       fsa_ctg_bridge  ${filtered_rawreads_list} ${dir_fsa}/contigs.fasta ${dir_align_contigs}/rawread2ctg.m4a ${dir_bridge_contigs}/bridged_contigs.fasta --ctg2ctg_file=${dir_align_contigs}/ctg2ctg.m4a --thread_size=${THREADS} --output_directory=${dir_bridge_contigs} ${FSA_CTG_BRIDGE_OPTIONS}        
       checkAndSaveRetCode $? "${dir_bridge_contigs}/finished" "Failed to run fsa_ctg_bridge"

    else
       
        if [[ ${step} == "3" ]]; then
            echo "Outputs is newer, so skip contig-bridging step"
        fi
    fi

done


if [[ -e ${dir_fsa}/contigs.fasta ]];then
    echo "The contig file: ${dir_fsa}/contigs.fasta"
    fsa_rd_stat ${dir_fsa}/contigs.fasta --thread_size=${THREADS}
    checkRetCode $? "Failed to run fsa_rd_stat contigs"
fi


if [[ -e ${dir_bridge_contigs}/bridged_contigs.fasta ]];then
    echo "The briged contig file: ${dir_bridge_contigs}/bridged_contigs.fasta"
    fsa_rd_stat ${dir_bridge_contigs}/bridged_contigs.fasta --thread_size=${THREADS}
    checkRetCode $? "Failed to run fsa_rd_stat bridged-contigs"
fi

