kittiTrackingBase="$1"
outputBase="$2"
timestampsfile="$3"
imgConfig="$4"
convertRawScript="./convertKittiRawToMapping.sh"

if [ ! -f "$convertRawScript" ]; then
  echo "$convertRawScript does not exist"
  exit 1
fi

if [[ -z "$kittiTrackingBase" || -z "$outputBase"  || -z "$timestampsfile"  || -z "$imgConfig" ]]; then
  echo "parameters not specified. usage:"
  echo "$0 [kitti_tracking_base_folder] [output_base_folder] [timestamps_file] [image_config_file]"
  exit 2
fi

if [[ ! -d "$outputBase" ]]; then
  echo "output folder does not exist"
  exit 3
fi

if [[ ! -f "$imgConfig" ]]; then
  echo "image_config_file does not exist"
  exit 4
fi

if [[ ! -f "$timestampsfile" ]]; then
  echo "timestamps_file does not exist"
  exit 5
fi

if [ ! -d "$kittiTrackingBase/training/velodyne/0000" ]; then
  echo "kitti_tracking_base must contain a /training/velodyne/0000 subfolder"
  exit 6
fi

for mode in training testing; do
  for seqPath in $kittiTrackingBase/$mode/velodyne/*; do
    seqNb=$(basename $seqPath)
    outdir="$outputBase/png_mapping_tracking_${mode}_${seqNb}"
    startIndex=`ls $seqPath | head -n 1 | sed -e 's/\.bin//'`
    endIndex=`ls $seqPath | tail -n 1 | sed -e 's/\.bin//'`
    echo "processing $seqPath -> $outdir"
    mkdir "$outdir"
    "$convertRawScript" "$seqPath" "$startIndex" "$endIndex" "$timestampsfile" "$outdir" "$imgConfig" > "png_mapping_tracking_${mode}_${seqNb}.log" 
  done
done
exit 0


