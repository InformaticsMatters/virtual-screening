#!/bin/bash
# arg1: file to split
# arg2: chunk size
# arg3: base filename for outputs

if [ ${1##*.} == 'gz' ]; then
  zcat $1 | sdsplit -$2 -omols_part_
else
  sdsplit -$2 -omols_part_ $1
fi

for f in mols_part_*.sd; do
  n=${f:10:-3}
  if [ ${#n} == 1 ]; then
    mv $f $3_0000${n}.sdf
  elif [ ${#n} == 2 ]; then
    mv $f $3_000${n}.sdf
  elif [ ${#n} == 3 ]; then
    mv $f $3_00${n}.sdf
  elif [ ${#n} == 4 ]; then
    mv $f $3_0${n}.sdf
  else
    mv $f $3_${n}.sdf
  fi
done