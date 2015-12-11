
########################################################################################################
### ChromHMM chromatin state calling 
########################################################################################################

cell_type=$1;

MEM=8;

#MODEL="/broad/compbio/anshul/projects/roadmap/segmentations/models/core_K27ac/parallel/set1/final/model_18_core_K27ac.txt";
#MODELNAME="observed_aux";
#NUMSTATES=18;
#MODEL="/broad/compbio/anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/final/model_15_coreMarks.txt";
MODEL="/broad/compbio/anshul/projects/roadmap/chromhmmSegmentations/ChmmModels/coreMarks/parallel/set2/final/model_15_coreMarks.txt";
MODELNAME="observed";
NUMSTATES=15;

# In/Output dir
FILEDIR="ENCODE3_state_calls/${cell_type}";
DATDIR="ENCODE3_data/${cell_type}";
mkdir -p ${FILEDIR}

#==================================
# Create a file with pointers to tagAlign files to use with ChromHMM
#==================================
#epitopes="H3K4me1 H3K4me3 H3K27me3 H3K36me3 H3K9me3 H3K27ac";
#epitopes="H3K4me1 H3K4me3 H3K27me3 H3K36me3 H3K9me3";
epitopes="H3K4me3 H3K4me1 H3K27me3 H3K9me3 H3K36me3";
rm -f ${FILEDIR}/metaDataFileList.txt
for epitope in ${epitopes}; do
  echo -e "${cell_type}\t${epitope}\tFINAL_${cell_type}_${epitope}.tagAlign.gz\tFINAL_${cell_type}_WCE.tagAlign.gz" >> \
    ${FILEDIR}/metaDataFileList.txt
done

#==================================
# Binarize tagAlign files
#==================================
#FRAGLEN_CONS=`cat ${DATDIR}/FINAL_*.qc | awk '{print $3}' | sed -r 's/,[^\t]+//g' | sort -n | uniq -c | sort -n | tail -n 1 | awk '{print $2}'`;
FRAGLEN_CONS=200;
#ChromHMM_scripts/chmm_binarize.sh ${FILEDIR}/metaDataFileList.txt $((${FRAGLEN_CONS}/2)) ${FILEDIR} ~/hg19.genome ${MEM} ${DATDIR} ${DATDIR}
java -mx$((MEM * 1024))M -jar $CHMM BinarizeBed -c ${DATDIR} -n $((${FRAGLEN_CONS}/2)) ~/hg19.genome ${DATDIR} ${FILEDIR}/metaDataFileList.txt ${FILEDIR}

#==================================
# Call chromatin states
#==================================
FILELIST="${FILEDIR}/filelist.txt";
find ${FILEDIR} -name "*_binary.txt" -exec basename {} \; | sort > ${FILELIST}
#ChromHMM_scripts/chmm_predict.sh ${FILEDIR} ${FILELIST} ${FILEDIR} CALLS ${MODEL} ${MEM}
java -mx$((MEM * 1024))M -jar $CHMM MakeSegmentation -b 200 -f ${FILELIST} -i CALLS -l ~/hg19.genome ${MODEL} ${FILEDIR} ${FILEDIR}

#ChromHMM_scripts/chmm_relabel.sh ${MODEL} ${FILEDIR} CALLS_REORDERED \
#  /broad/compbio/anshul/projects/roadmap/segmentations/models/core_K27ac/parallel/set1/final/colormap_18_core_K27ac.tab

mkdir -p ${FILEDIR}/STATEBYLINE
#ChromHMM_scripts/chmm_statesbyline.sh ${FILEDIR} ${FILELIST} ${FILEDIR} CALLS_PER_LINE ${MODEL} ${MEM}
java -mx$((MEM * 1024))M -jar $CHMM MakeSegmentation -b 200 -nobed -printstatesbyline -f ${FILELIST} \
  -i CALLS_PER_LINE -l ~/hg19.genome ${MODEL} ${FILEDIR} ${FILEDIR}
gzip -f ${FILEDIR}/STATEBYLINE/*_statebyline.txt

mkdir -p freqs/${MODELNAME};
zcat ${FILEDIR}/STATEBYLINE/${cell_type}_${NUMSTATES}_*_statebyline.txt.gz | grep "^[0-9]" | \
  sort -n | uniq -c | awk '{print $2, $1}' > freqs/${MODELNAME}/table_${cell_type}.txt

echo "Chromatin state frequencies:";
cat freqs/${MODELNAME}/table_${cell_type}.txt



