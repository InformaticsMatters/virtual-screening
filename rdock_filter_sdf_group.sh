#!/bin/bash

set -e

if [[ $# -lt 4 || $# -gt 6 ]]; then
  echo "Usage: rdock_filter_sdf_group.sh infile.sdf outfile.sdf sort_field sort_descending [top] [group_by_field]"
  exit 1
fi

# args:
INFILE=$1
OUTFILE=$2
SORT_FIELD=$3      # e.g. SCORE
SORT_DESCENDING=$4 # true or false
if [ $# -gt 4 ]; then
  TOP=$5           # e.g 100
else
  TOP="0"
fi
if [ $# -gt 5 ]; then
  GROUP_BY_FIELD=$6  # e.g. std_smi
else
  GROUP_BY_FIELD=_TITLE1
fi

if "$SORT_DESCENDING"; then SORT_FLAG="-r"; else SORT_FLAG=""; fi

DIR=$(dirname $OUTFILE)
mkdir -p $DIR

if [ $TOP -gt 0 ]; then
  echo "Starting grouped sort with top $TOP ..."
  filt="\$_REC <= $TOP"
  sdsort -n -s -f"$SORT_FIELD" -id"$GROUP_BY_FIELD" $SORT_FLAG "$INFILE" |
    sdfilter -f'$_COUNT == 1' -s"$GROUP_BY_FIELD" |
    sdsort -n -f"$SORT_FIELD" $SORT_FLAG |
    sdfilter -f"$filt"> "$OUTFILE"
else
  echo "Starting grouped sort ..."
  sdsort -n -s -f"$SORT_FIELD" -id"$GROUP_BY_FIELD" $SORT_FLAG "$INFILE" |
    sdfilter -f'$_COUNT == 1' -s"$GROUP_BY_FIELD" |
    sdsort -n -f"$SORT_FIELD" $SORT_FLAG > "$OUTFILE"
fi

IN_COUNT=$(fgrep -c '$$$$' $INFILE)
OUT_COUNT=$(fgrep -c '$$$$' $OUTFILE)
echo "Filtered $IN_COUNT records down to $OUT_COUNT"
