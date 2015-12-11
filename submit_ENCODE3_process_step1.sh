#!/bin/bash

mkdir -p out;

IFS=$'\t'
while read ID cell_type epitope BAMfile reads_aligned reads_duplicate; do
  cell_type=${cell_type//\ /-};
  echo "$cell_type / $epitope";
  if [ "${cell_type}" != "cell_type" ]; then 
    bsub -q priority -P compbiofolk -J ENCODE3_process_step1_${ID}_${cell_type}_${epitope} \
         -oo out/output_ENCODE3_process_step1_${ID}_${cell_type}_${epitope}.out \
         -R "rusage[mem=8]" ./code_ENCODE3_process_step1.sh ${ID} ${cell_type} ${epitope} ${BAMfile};
  fi;
done < track_info.txt


