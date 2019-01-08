#!/bin/bash
datafolder="$1" #"$inputFolder/data"
startIdx="$2"
endIdx="$3"
timestampsfile="$4" #"$inputFolder/timestamps.txt"
outputFolder="$5"
imgConfig="$6"
convertExe="./kitti_convert_velodyne"

if [ ! -f "$convertExe" ]; then
  echo "$convertExe does not exist"
  exit 1
fi

if [[ -z "$datafolder" || -z "$startIdx" || -z "$endIdx" || -z "$timestampsfile" || -z "$outputFolder"  || -z "$imgConfig" ]]; then
  echo "parameters not specified. usage:"
  echo "$0 [velodyne_points_folder] [start_index] [end_index] [timestamps_file] [output_folder] [image_config_file]"
  exit 2
fi

if [[ ! -d "$outputFolder" ]]; then
  echo "output folder does not exist"
  exit 3
fi

if [[ ! -f "$imgConfig" ]]; then
  echo "image_config_file does not exist"
  exit 4
fi

if [ "${endIdx}" -lt "${startIdx}" ]; then
  echo "start_index must be smaller than end_index"
  exit 5
fi

if [[ ! -d "$datafolder" || ! -f "$timestampsfile"  || ! -f "$datafolder/$startIdx.bin" ]]; then
  echo "input specification is invalid. $timestampsfile and $datafolder/$startIdx.bin must exist"
  exit 6
fi

formatlength=${#startIdx}
for index in $(seq -f "%0${formatlength}g" $startIdx $endIdx)
do
  sourcefile="$datafolder/$index.bin"
  if [ ! -f "$sourcefile" ]; then
    echo "WARNING: file $sourcefile does not exist. skipping."
    continue
  fi
  targetfile="$outputFolder/scan$index.png"
  "${convertExe}" "${sourcefile}" "${targetfile}" "${imgConfig}"
done

cp "${imgConfig}" "${outputFolder}/img.cfg"

targetTimestampsfile="${outputFolder}/imu.cfg"
cat "imuheader4mapping.txt" > "${targetTimestampsfile}"
#cat "${timestampsfile}" >> "${targetTimestampsfile}"
sed -n $((startIdx+1)),$((endIdx+1))p "${timestampsfile}" >> "${targetTimestampsfile}"
