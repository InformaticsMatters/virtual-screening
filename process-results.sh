#!/bin/bash

REPORT_FIELDS='SCORE,SCORE.INTRA,SCORE.INTER,SCORE.norm,SCORE.INTRA.norm,SCORE.INTER.norm,enum_code,enum_smi'
: "${SORT_FIELD:=SCORE.INTER}"
: "${BASE_NAME:=results_rdock}"

#echo "Generating full report"
#sdreport -c$REPORT_FIELDS ${BASE_NAME}.sdf > ${BASE_NAME}.csv
echo "Reduce to 1 pose per ligand"
sdsort -n -s -f$SORT_FIELD -idstd_smi ${BASE_NAME}.sdf | sdfilter -f'$_COUNT == 1' -sstd_smi > ${BASE_NAME}_1poseperlig_${SORT_FIELD}.sdf
echo "Generating one pose per ligand report"
sdreport -c$REPORT_FIELDS ${BASE_NAME}_1poseperlig_${SORT_FIELD}.sdf > ${BASE_NAME}_1poseperlig_${SORT_FIELD}.csv
echo "Sorting by $SCORE_FIELD"
sdsort -n -f$SORT_FIELD ${BASE_NAME}_1poseperlig_${SORT_FIELD}.sdf > ${BASE_NAME}_1poseperlig_${SORT_FIELD}_sorted.sdf
echo "Generating top 50 SDF"
sdfilter -f'$_REC <= 50' ${BASE_NAME}_1poseperlig_${SORT_FIELD}_sorted.sdf > ${BASE_NAME}_1poseperlig_${SORT_FIELD}_sorted_top50.sdf
echo "Generating top 50 report"
sdreport -c$REPORT_FIELDS ${BASE_NAME}_1poseperlig_${SORT_FIELD}_sorted_top50.sdf > ${BASE_NAME}_1poseperlig_${SORT_FIELD}_sorted_top50.csv
