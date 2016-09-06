#!/bin/bash
if [ -n "$1" ]; then
  BASEDIR="$1"
fi

if [ ! -d "$BASEDIR" ]; then
  echo "usage: $0 <BASEDIR>"
  exit 1
fi

# files hold in each line: [TrackID] Age LastUp R G B state12D covar12x12 nbVIPPts nbAllPts Nx(surfaceinfo)
# index                    1         2   3      4 5 6 7        19
{
for f in ${BASEDIR}/*.png.txt; do
  # assume world track [1] is always in first line --> store pose
  worldpose=`head -n 1 $f | awk '{printf "%s %s %s %s %s %s",$7,$8,$9,$10,$11,$12}'`
  # for each line in file extract ID Age pose and stored worldpose
  awkcmd="{printf \"%s %s %s %s %s %s %s %s ${worldpose}\n\",\$1,\$2,\$7,\$8,\$9,\$10,\$11,\$12}"
  cat $f | awk "$awkcmd" | sed -e 's/\[//g' -e 's/\]//g'  # remove [] around ID
done
} | sort -k 1,1 -k 2,2 -n #and sort by ID then by age
