#!/bin/bash
BASEINDIR="trackingResults"
BASEOUTDIR="trackingScreenshots"
if [ "$1" == "--force" ]; then
  FORCE=1
fi
if [[ "$1" == "--help" ]]; then
  echo "usage: $0 [--force]"
  exit 0
fi
echo "FORCE=$FORCE"

for indir in $BASEINDIR/2012*; do
  if [ ! -d $indir ]; then continue; fi
  RCMD="s/${BASEINDIR}/${BASEOUTDIR}/"
  outdir=`echo $indir|sed -e $RCMD`
  if [[ -d $outdir && "$FORCE" != "1" ]]; then echo "skipping $indir"; continue; fi
  echo "processing $indir -> $outdir"
  rm -rf $outdir
  mkdir $outdir
  ./trackVisualizer $indir --pilotBuffer ~/q7_seq.cam --autogenScreenshots $outdir > /dev/null
done
