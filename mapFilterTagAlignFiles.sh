#!/bin/bash

if [[ "$#" -lt 2 ]]
then
    echo 'batch filters tagAlign files using mappability tracks'  1>&2
    echo 'USAGE:'  1>&2
    echo "$(basename $0) [FromDir] [ToDir]"  1>&2
    echo "e.g. $(basename $0) /dir/from /dir/to" 1>&2
    exit 1
fi

IDIR=$(echo $1 | sed -r 's:/$::g')
ODIR=$(echo $2 | sed -r 's:/$::g')
JOBGROUPNAME='/compbio/short/filtTagAlign'

# Check them input directory exists

if [[ ! -d ${IDIR} ]]
then
    echo "From directory ${IDIR} does not exist" 1>&2
    exit 1
fi

# Make output directory if it doesnt exist
[[ ! -d ${ODIR} ]] && mkdir ${ODIR}

# Replicate input directory structure in output directory
for dname in $(find ${IDIR} -type d -name '*')
do
    newdname=$(echo ${dname} | sed -r 's:'"${IDIR}"':'"${ODIR}"':g')
    if [[ -z ${newdname} ]]
    then
	continue
    fi
    [[ ! -d ${newdname} ]] && mkdir ${newdname}
done

count=0
for fname in $(find "${IDIR}" -name '*.tagAlign.gz')
do
    SDIR="${SEQDIR}/encodeHg19Male"
    UDIR="${UMAPDIR}/encodeHg19Male/globalmap_k20tok54"
    
  # Set output and log file name
    OFNAME=$(echo ${fname} | sed -r -e 's:'"${IDIR}"':'"${ODIR}"':g' -e 's/\.tagAlign\.gz$/.filt.tagAlign/g')
    LOGFILE=$(echo ${fname} | sed -r -e 's:'"${IDIR}"':'"${ODIR}"':g' -e 's/\.tagAlign\.gz$/.filt.tagAlign.logfile/g')
    
    if [[ -e "${OFNAME}" || -e "${OFNAME}.gz" || -e "${LOGFILE}" ]]
    then
	echo "Skipping $(basename ${fname})" >> skipped.log
	continue
    fi
     
  # Check datasets in job queue
  #  while [[ $(bjobs -g ${JOBGROUPNAME} 2> /dev/null | wc -l) -gt 50 ]]
  #  do
  #	echo 'Waiting..'
  #	sleep 30s
  #  done
  sleep 3s
  
  # Create submit script
    submitScriptName="submit${RANDOM}.sh"
    echo '#!/bin/bash' > ${submitScriptName}
    echo "export TMP=${TMP}/tmp_${RANDOM}${RANDOM}" >> ${submitScriptName}
    echo 'mkdir ${TMP}' >> ${submitScriptName}
    echo 'export MCR_CACHE_ROOT=${TMP}' >> ${submitScriptName}
    echo "filterUniqueReads -s=${SDIR} -u=${UDIR} -v=${LOGFILE} ${fname} | grep -v 'Warning' | gzip -c > ${OFNAME}.gz" >> ${submitScriptName}
    echo 'rm -rf ${TMP}' >> ${submitScriptName}
    
    JOBNAME="filt$(basename ${fname})"
    #bsub -P compbiofolk -g ${JOBGROUPNAME} -J${JOBNAME} -o out.log -e out.log < ${submitScriptName}
    bsub -P compbiofolk -q compbio-week -J${JOBNAME} -M 16 -R "rusage[mem=16]" -o ${LOGFILE} -e ${LOGFILE} < ${submitScriptName}
    rm ${submitScriptName}
done
