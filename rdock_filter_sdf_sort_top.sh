#!/bin/bash

set -e

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: rdock_filter_sdf_sort_top.sh infile.sdf outfile.sdf sort_field sort_descending [top]"
  exit 1
fi

# args:
INFILE=$1
OUTFILE=$2
SORT_FIELD=$3      # e.g. SCORE
SORT_DESCENDING=$4 # true or false
TOP=$5             # e.g. 100

if "$SORT_DESCENDING"; then SORT_FLAG="-r"; else SORT_FLAG=""; fi

DIR=$(dirname $OUTFILE)
mkdir -p $DIR

echo "Starting sort ..."

if [ $# -eq 5 ] && [ $TOP -gt "0" ]; then
  filt="\$_REC <= $TOP"
  echo "Keeping top $TOP"
  sdsort -n -f"$SORT_FIELD" $SORT_FLAG "$INFILE" | sdfilter -f"$filt" > "$OUTFILE"
else
  echo "No top - keeping all"
  sdsort -n -f"$SORT_FIELD" $SORT_FLAG "$INFILE" > "$OUTFILE"
fi

IN_COUNT=$(fgrep -c '$$$$' $INFILE)
OUT_COUNT=$(fgrep -c '$$$$' $OUTFILE)
echo "Filtered $IN_COUNT records down to $OUT_COUNT"
