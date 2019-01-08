echo this file is outdated
exit 1
d="$1" #2011_10_03/2011_10_03_drive_0042_sync
idxS="$2"
idxE="$3" #0000001100
onb="$4" #02
./convertKittiRawToMapping.sh "/lhome/fmoosma/sequences/kitti/kitti_raw/${d}/velodyne_points" "$idxS" "$idxE" "/lhome/fmoosma/sequences/kitti/png_mapping_odometry_$onb" "/lhome/fmoosma/code/FranksVelodyneAlgos/exampleData/DurlacherTor/img.cfg"
