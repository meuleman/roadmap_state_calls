#!/bin/bash -l

###################################################
### ChromHMM chromatin state calling 
###################################################

cell_type=$1;
CELL_DIR=$2
CC_DIR=$3

MEM=8;
NUMSTATES=18;

if [[ "$NUMSTATES" == "18" ]] 
then 
    MODEL="/broad/compbio/anshul/projects/roadmap/segmentations/models/core_K27ac/parallel/set1/final/model_18_core_K27ac.txt";
    MODELNAME="observed_aux";
    epitopes="H3K4me1 H3K4me3 H3K27me3 H3K36me3 H3K9me3 H3K27ac";
else 
    # MODEL="/broad/compbio/anshul/projects/roadmap/segmentations/models/coreMarks/parallel/set2/final/model_15_coreMarks.txt";
    MODEL="/broad/compbio/anshul/projects/roadmap/chromhmmSegmentations/ChmmModels/coreMarks/parallel/set2/final/model_15_coreMarks.txt";
    MODELNAME="observed";
    NUMSTATES=15;
    epitopes="H3K4me3 H3K4me1 H3K27me3 H3K9me3 H3K36me3";
fi

#==================================
# Create a file with pointers to tagAlign files to use with ChromHMM
#==================================
#epitopes="H3K4me1 H3K4me3 H3K27me3 H3K36me3 H3K9me3 H3K27ac";
#epitopes="H3K4me1 H3K4me3 H3K27me3 H3K36me3 H3K9me3";
#epitopes="H3K4me3 H3K4me1 H3K27me3 H3K9me3 H3K36me3";

rm -f ${CC_DIR}/metaDataFileList.txt
for epitope in ${epitopes} 
do
    echo -e "${cell_type}\t${epitope}\tFINAL_${cell_type}_${epitope}.tagAlign.gz\tFINAL_${cell_type}_WCE.tagAlign.gz" >> ${CC_DIR}/metaDataFileList.txt
done

#==================================
# Binarize tagAlign files
#==================================
#FRAGLEN_CONS=`cat ${CELL_DIR}/FINAL_*.qc | awk '{print $3}' | sed -r 's/,[^\t]+//g' | sort -n | uniq -c | sort -n | tail -n 1 | awk '{print $2}'`;
FRAGLEN_CONS=200;
#ChromHMM_scripts/chmm_binarize.sh ${CC_DIR}/metaDataFileList.txt $((${FRAGLEN_CONS}/2)) ${CC_DIR} ~/hg19.genome ${MEM} ${CELL_DIR} ${CELL_DIR}
java -mx$((MEM * 1024))M -jar $CHMM BinarizeBed -c ${CELL_DIR} -n $((${FRAGLEN_CONS}/2)) ~/hg19.genome ${CELL_DIR} ${CC_DIR}/metaDataFileList.txt ${CC_DIR}

#==================================
# Call chromatin states
#==================================
FILELIST="${CC_DIR}/filelist.txt";
find ${CC_DIR} -name "*_binary.txt" -exec basename {} \; | sort > ${FILELIST}
#ChromHMM_scripts/chmm_predict.sh ${CC_DIR} ${FILELIST} ${CC_DIR} CALLS ${MODEL} ${MEM}
java -mx$((MEM * 1024))M -jar $CHMM MakeSegmentation -b 200 -f ${FILELIST} -i CALLS -l ~/hg19.genome ${MODEL} ${CC_DIR} ${CC_DIR}

#ChromHMM_scripts/chmm_relabel.sh ${MODEL} ${CC_DIR} CALLS_REORDERED \
#  /broad/compbio/anshul/projects/roadmap/segmentations/models/core_K27ac/parallel/set1/final/colormap_18_core_K27ac.tab

mkdir -p ${CC_DIR}/STATEBYLINE
#ChromHMM_scripts/chmm_statesbyline.sh ${CC_DIR} ${FILELIST} ${CC_DIR} CALLS_PER_LINE ${MODEL} ${MEM}
java -mx$((MEM * 1024))M -jar $CHMM MakeSegmentation -b 200 -nobed -printstatesbyline -f ${FILELIST} \
  -i CALLS_PER_LINE -l ~/hg19.genome ${MODEL} ${CC_DIR} ${CC_DIR}
gzip -f ${CC_DIR}/STATEBYLINE/*_statebyline.txt

mkdir -p freqs/${MODELNAME};
zcat ${CC_DIR}/STATEBYLINE/${cell_type}_${NUMSTATES}_*_statebyline.txt.gz | grep "^[0-9]" | \
  sort -n | uniq -c | awk '{print $2, $1}' > freqs/${MODELNAME}/table_${cell_type}.txt

echo "Chromatin state frequencies:";
cat freqs/${MODELNAME}/table_${cell_type}.txt



