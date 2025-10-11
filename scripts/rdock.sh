#!/bin/bash

while getopts ":r:p:a:n:c:m:o:" opt; do
  case $opt in
    r)
      RECEPTOR=$OPTARG
      echo "setting RECEPTOR to $RECEPTOR"
      ;;
    p)
      PRM=$OPTARG
      echo "setting PRM to $PRM"
    ;;
    n)
      COUNT=$OPTARG
      echo "setting COUNT to $COUNT"
    ;;
    a)
      AS=$OPTARG
      echo "setting AS to $AS"
      ;;
    c)
      MOLS=$OPTARG
      echo "setting MOLS to $MOLS"
      ;;
    m)
      MODE=$OPTARG
      echo "setting MODE to $MODE"
      ;;
    o)
      OUTFILE=$OPTARG
      echo "setting OUTFILE to $OUTFILE"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [ -z ${AS+x} ]; then
  echo "active site file (.as) not specified"
  exit 1
elif
  [[ "$AS" == 'docking.as' ]]; then
    echo "Using supplier docking.as file";
  else
    echo "Using supplied .as file $AS";
    cp $AS docking.as
  fi

if [ -z ${PRM+x} ]; then
  echo "generating .prm file"
  echo -e "RBT_PARAMETER_FILE_V1.00\\nTITLE rDockdocking\\nRECEPTOR_FILE @@RECEPTOR_FILE@@\\nRECEPTOR_FLEX 3.0\\n" > docking.prm
elif [[ "$PRM" == 'docking.prm' ]]; then
  echo "Using supplier docking.prm file"
else
  echo "Using supplied .prm file $PRM"
  cp $PRM docking.prm
fi

sed -i "s/@@RECEPTOR_FILE@@/$RECEPTOR/" docking.prm


echo "running rdock"
rbdock -r 'docking.prm' -p $MODE.prm -n 1 -i $MOLS -o $OUTFILE > $OUTFILE.log
mv $OUTFILE.sd $OUTFILE.sdf