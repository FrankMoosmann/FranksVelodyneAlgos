#!/usr/bin/python
# encoding: utf-8

from __future__ import division, print_function

import sys
import csv
import pprint
import glob
import math

sys.path.append('/home/moosmann/code/python/libTransformation/')
from transformations import *
from myTransformations import *

# function checks if move from center-position "c0_W" to "ci_W" (in world-coordinates) exceeds "minMove"
# if yes, the normalized move vector is calculated, rotated into track-CS using "H0_OW"
# and stored in global variable "moves"
def trydetect(currid, ci_W, c0_W, H0_OW, minMove):
  global moves  # dictionary from trackID to move-vector
  move_W_3D = (ci_W - c0_W)[:3]
  dist = numpy.linalg.norm(move_W_3D)
  if (dist > minMove):
    sys.stdout.write('%s(%f)...' % (currid, dist))
    #print(rowIdx, currid, srow[1], move_W_3D, "-->", dist, " | ", movenorm_O_3D )
    #print("     DETECT  ", movenorm_O_3D, "-->", dist)
    move_O_3D = numpy.dot(H0_OW[0:3,0:3],move_W_3D) # rotate only
    movenorm_O_3D = move_O_3D/numpy.linalg.norm(move_O_3D)
    moves[int(currid)] = movenorm_O_3D # store move
    return 1
  return 0

# detects outliers in a sequence, where maxOutlierRatio is a maximum
# portion of outliers on each side in the interval (0,0.5) 
# returns (lowVal, highVal), the entries that mark the valid data-range
def detectOutlier(seq, maxOutlierRatio):
  # return seq[0], seq[len(seq)-1] # no outliers
  seq.sort()
  nbPts=len(seq) # should be equivalent to row[163], check??
  QLOW=int(nbPts*maxOutlierRatio)
  QHIG=int(nbPts*(1-maxOutlierRatio))
  #DRANGE=min(1,int(QLOW/2))
  DRANGE=QLOW
  DETECT=0.3
  vlow=seq[0]
  for l in range(1,QLOW):
    okdiff= (seq[l+DRANGE]-seq[l])
    totdiff=(seq[l+DRANGE]-seq[l-1])
    if (totdiff > 0.0) and (okdiff/totdiff < DETECT):
      #print("detect at",l,"of",QLOW,"of",nbPts,"with diffs",okdiff,totdiff)
      vlow=seq[l]
      break
  vhig=seq[nbPts-1]
  for l in range(nbPts-2,QHIG,-1):
    okdiff= (seq[l]  -seq[l-DRANGE])
    totdiff=(seq[l+1]-seq[l-DRANGE])
    if (totdiff > 0.0) and (okdiff/totdiff < DETECT):
      #print("detect at",l,"of",QLOW,"of",nbPts,"with diffs",okdiff,totdiff)
      vhig=seq[l]
      break
  return vlow,vhig


# parse command line
if len(sys.argv) < 2:
  sys.exit('Usage: %s track-directory' % sys.argv[0])

basedir=sys.argv[1]
if len(sys.argv) > 2:
  tracesfile=sys.argv[2]
else:
  tracesfile=basedir+'/tracktraces.dat'
bbext=".bb"

if not os.path.isdir(basedir):
  sys.exit('Error: directory %s not found' % basedir)
  
if not os.path.isfile(tracesfile):
  sys.exit('Error: file %s not found' % tracesfile)


# global variables
origin, xaxis, yaxis, zaxis = [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
OPT_MOVE=4.0 # optimally calculate move direction after this number of travel distance
MIN_MOVE=0.5 # if moved less, calculate anyway but skip calulation if below this number
moves = {} # {1: [1, 0, 0, 1]} # dictionary from trackID to move-vector

# determine main direction of travel for each track-ID relative to the object coordinate system
# parse tracktraces which holds a list of track details with following information
# ID age rx ry rz tx ty tz wrx wry wrz wtx wty wtz (w=worldtrack) --> transformation by t->rz->rz->rx 
# the list is sorted by ID and age, hence the details for a specific track can be parsed in sequence
sys.stdout.write('determining main direction of travel...')
tracks = csv.reader(open(tracesfile, 'r'), delimiter=' ')
oldid=0
rowIdx=0
detect=0
nbData=0
for srow in tracks:
  #print(srow)
  rowIdx += 1
  currid=srow[0] # string
  if (currid == "1"): # skip world track
    detect=1
    continue
  if (currid != oldid) and (detect == 0):
    # if detection for old ID did not succeed with optimal move-distance (see below),
    # try to detect with minimum move-distance using variable stated from last loop
    detect = trydetect(oldid, ci_W, c0_W, H0_OW, MIN_MOVE)
  if (currid != oldid) and (oldid != 0):
    # print move in WC for statistical reasons
    totDist = numpy.linalg.norm((ci_W - c0_W)[:3])
    sys.stdout.write('%s[%f|%d]...' % (oldid, totDist, nbData))
    # reset data count
    nbData=0
  row = [float(x) for x in srow] # convert to float
  # prepare transformations
  Hi_SO = YawPitchRollXYZ_2_HTM    (row[4], row[3], row[2], row[5], row[6], row[7])
  Hi_WS = YawPitchRollXYZ_2_HTM_inv(row[10], row[9], row[8], row[11], row[12], row[13])
  Hi_WO = numpy.dot(Hi_WS,Hi_SO) # object CS to world CS
  ci_W = numpy.dot(Hi_WO,[0,0,0,1]) # object center in world CS
  nbData = nbData + 1
  if (currid != oldid):
    #print(rowIdx, currid, srow[1])
    # TODO: not only keep first center, but whole sequence and calculate direction of travel for each frame independently using past sequence
    c0_W = ci_W
    H0_OW = numpy.linalg.inv(Hi_WO)
    detect = 0
  else:
    if (detect == 0):
      # try to detect with optimal move-distance
      detect = trydetect(currid, ci_W, c0_W, H0_OW, OPT_MOVE)
  oldid=currid
totDist = numpy.linalg.norm((ci_W - c0_W)[:3])
sys.stdout.write('%s[%f|%d]...' % (oldid, totDist, nbData))
print('done')
#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(moves)
#print(moves.keys())

# determine bounding boxes in each frame for each track wrt. the sensor coordinate system
# hence, parse all frame-details-files in directory and generate bounding-box-files
# files hold in each line: [TrackID] Age LastUp R G B state12D covar12x12 nbVIPPts nbAllPts Nx(surfaceinfo=points+normal+NC in sensor-CS)
# index                    0         1   2      3 4 5 6        18         162      163      164
fileIdx=0
print("")
fstr='\rcalculating bounding box for each frame...%d%%'
trackfiles=sorted(glob.glob( os.path.join(basedir, '*.txt')))
# loop over trackfiles
for tracksfile in trackfiles:
  sys.stdout.write(fstr % int(100*fileIdx/len(trackfiles)))
  sys.stdout.flush()
  fileIdx += 1
  #if (fileIdx > 5):
  #  break
  bbfile=tracksfile+'.bb'
  #print("processing ",tracksfile," -> ",bbfile)
  trackdetails = csv.reader(open(tracksfile, 'r'), delimiter=' ')
  bbout = open(bbfile, 'w')
  # loop over tracks in trackfile
  for srow in trackdetails:
    srow[0] = srow[0].replace('[','').replace(']','') # remove brackets around ID
    currid=srow[0] # string
    #print(currid)
    row = [float(x) for x in srow] # convert to float
    # determine direction of motion wrt sensor CS
    if (moves.has_key(int(currid)) and (len(row) > 164)):
      Hi_SO = YawPitchRollXYZ_2_HTM(row[8], row[7], row[6], row[9], row[10], row[11])
      movenorm_O = moves[int(currid)]
      movenorm_S = numpy.dot(Hi_SO[0:3,0:3],movenorm_O)
      alpha = math.atan2(movenorm_S[1],movenorm_S[0])
      Ralpha = rotation_matrix(alpha, zaxis)
      RalphaI = numpy.linalg.inv(Ralpha)
      xpos = []
      ypos = []
      zpos = []
      if (len(srow) != 164+(row[163]*7)):
        print("ERROR: len=",len(srow)," but [163]=",row[163]," --> last=",164+(row[163]*7))
      for pti in range(0,int(row[163])):
        idx = 164+(pti*7)
        pos = [row[idx], row[idx+1], row[idx+2], 1] #homogeneous point in sensorCS
        pos = numpy.dot(RalphaI,pos)
        xpos.append(pos[0])
        ypos.append(pos[1])
        zpos.append(pos[2])
      xlow,xhigh = detectOutlier(xpos, 0.35)
      ylow,yhigh = detectOutlier(ypos, 0.35)
      zlow,zhigh = detectOutlier(zpos, 0.35)
      # TODO: detect outlier not independently for each dimension but together
      minpos=numpy.array([xlow, ylow, zlow])
      maxpos=numpy.array([xhigh, yhigh, zhigh])
      extents=(maxpos-minpos)/2.0
      center=(maxpos+minpos)/2.0
      center=numpy.dot(Ralpha[0:3,0:3],center)#rotate back center
      #bbout.write('%s %f \n' % (currid % minpos[0]) )
      bbout.write('%f %f %f %f %f %f %f \n' % (center[0], center[1], center[2], alpha, extents[0], extents[1], extents[2]) )
    else:
      bbout.write('\n') # empty line
  bbout.close()
print("")
