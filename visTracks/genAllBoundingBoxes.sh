#!/bin/bash
BASEDIR=$1
if [ "$2" == "--force" ]; then
  FORCE=1
fi
if [[ -z "$BASEDIR" || ! -d $BASEDIR ]]; then
  echo "usage: $0 <basedir> [--force]"
  exit 1
fi

TRACESFILE="tracktraces.dat"
GENSCRIPT="./genBoundingBoxes.sh"
for dir in $BASEDIR/*; do
  if [ ! -d $dir ]; then continue; fi
  NBT=`ls $dir/*.png.txt | wc -l` # number of track files
  NBB=`ls $dir/*.png.txt.bb | wc -l` # number of bounding-box-files
  if [[ -f $dir/$TRACESFILE && "$FORCE" != "1" && "$NBT" == "$NBB" ]]; then continue; fi
  #echo "processing $dir"
  $GENSCRIPT $dir
done
