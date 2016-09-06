#!/bin/bash
BASEDIRDEFAULT="../KogmoImageWiseProcessing/tracks"
BASEDIR=$BASEDIRDEFAULT
if [ -n "$1" ]; then
  BASEDIR="$1"
fi

if [ ! -d "$BASEDIR" ]; then
  echo "usage: $0 <BASEDIR=$BASEDIRDEFAULT>"
  exit 1
fi


#TRACESFILE=`mktemp`
TRACESFILE="$BASEDIR/tracktraces.dat"
echo "processing $BASEDIR"
sleep 1
echo "##################  extracting track-traces  ##################"
./xtractTrackTraces.sh $BASEDIR > $TRACESFILE
echo "################## extracting bounding boxes ##################"
./xtractBoundingBoxes.py $BASEDIR $TRACESFILE
