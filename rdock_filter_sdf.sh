#!/bin/bash

set -e

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: filter_sdf.sh infile.sdf outfile.sdf sort_field sort_descending [group_by_field]"
  exit 1
fi

# args:
INFILE=$1
OUTFILE=$2
SORT_FIELD=$3      # e.g. SCORE
SORT_DESCENDING=$4 # true or false
if [ $# -eq 5 ]; then
  GROUP_BY_FIELD=$5  # e.g. std_smi
else
  GROUP_BY_FIELD=_TITLE1
fi

if "$SORT_DESCENDING"; then SORT_FLAG="-r"; else SORT_FLAG=""; fi

DIR=$(dirname $OUTFILE)
mkdir -p $DIR

sdsort -n -s -f"$SORT_FIELD" -id"$GROUP_BY_FIELD" $SORT_FLAG "$INFILE" |
  sdfilter -f'$_COUNT == 1' -s"$GROUP_BY_FIELD" |
  sdsort -n -f"$SORT_FIELD" $SORT_FLAG > "$OUTFILE"
