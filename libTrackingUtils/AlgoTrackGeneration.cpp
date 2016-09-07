#include "AlgoTrackGeneration.hpp"

#include <cfloat>
#include <climits>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <map>
#include <boost/foreach.hpp>
#include <boost/array.hpp>
#include <boost/filesystem.hpp>
#include <boost/typeof/std/utility.hpp>

#include "AlgoFrameFeatures.hpp"
#include "AlgoSegmentation.hpp"
#include "ParameterHeap.hpp"
#include "CsvReader.hpp"
#include "Classifier.hpp"
#include "SvmClassifier.hpp"
#include "FeatureVector.hpp"
#include "FeatureMaxRangeNormalizer.hpp"

using namespace std;
using namespace matrixTools;
using namespace boost::filesystem;

// copies tracks from last frame to current frame and calls kalmanfilter.predict()
void predictTracks(FrameSPtr lf, FrameSPtr cf)
{
  ParameterHeap* params = ParameterHeap::get();
  const bool deepCopyTracks = params->deepCopyTracks;
  //const double trkDelVarThresh = params->trkDelVarThresh;
  const double trkDelUpCntThresh = params->trkDelUpCntThresh;
  const double trkDelMaxDist = params->trkDelMaxDist;
  const unsigned int minPtCnt = params->segMinSizePt;

  ///////////////////////////////////////////////////////////////
  // make copy of the short-term-tracks of lastFrame and predict
  ///////////////////////////////////////////////////////////////
  assert(cf->shortTermTracks.size() == lf->shortTermTracks.size() && "predictTracks: shortTermTracks have different buffer sizes!");
  const unsigned int nbVStep = cf->sttValidationSteps;
  vector<LidarFrame::ExtTrackList> &cfStt = cf->shortTermTracks;
  vector<LidarFrame::ExtTrackList> &lfStt = lf->shortTermTracks;
  unsigned int copiedSTTrackCount = 0;
  // don't copy buffer [>=nbVStep] since these are contained in long-term-tracks --> handled below
  for (unsigned int bi=0; bi<nbVStep; ++bi) {
    cfStt[bi].linksT1 = lfStt[bi].linksT1; // graph is on indices which didn't change, so just re-use it
    cfStt[bi].linksT2 = lfStt[bi].linksT2; // graph is on indices which didn't change, so just re-use it
    cfStt[bi].lttLinks = lfStt[bi].lttLinks; // graph is on indices which didn't change, so just re-use it
    cfStt[bi].tracks.clear();
    BOOST_FOREACH(PointCloudTrack::SPtr tsp, lfStt[bi].tracks) {
      if (deepCopyTracks)
        tsp = PointCloudTrack::SPtr(new PointCloudTrack(*tsp.get()));
      tsp->predict();
      cfStt[bi].tracks.push_back(tsp);
      ++copiedSTTrackCount;
    } // for each track in current buffer
    assert((cfStt[bi].tracks.size() == lfStt[bi].tracks.size()) && "predictTracks: corresponding shortTermTracks buffers have different track counts!");
  } // for each buffer

  /////////////////////////////////////////////////////////////////////////////////
  // make copy of the long-term-tracks of lastFrame and predict + decide on delete
  /////////////////////////////////////////////////////////////////////////////////
  unsigned int copiedTrackCount = 0;
  unsigned int totalTrackCount = 0;
  cf->longTermTracks.clear();
  cfStt[nbVStep].tracks.clear();
  cfStt[nbVStep].tracks.resize(lfStt[nbVStep].tracks.size()); // must be same size - initializes with NULL pointers
  cfStt[nbVStep+1].tracks.clear();
  cfStt[nbVStep+1].tracks.resize(lfStt[nbVStep+1].tracks.size()); // must be same size - initializes with NULL pointers
  BOOST_FOREACH(PointCloudTrack::SPtr tsp, lf->longTermTracks) {
    PointCloudTrack::SPtr tspNew = tsp;
    if (deepCopyTracks)
      tspNew = PointCloudTrack::SPtr(new PointCloudTrack(*tsp.get()));
    ++totalTrackCount;
    tspNew->removeFarPoints(trkDelMaxDist);
    DMatrix tCovar; tspNew->getCovar(tCovar);
    double mPVar = max(max(tCovar(3,3),tCovar(4,4)),tCovar(5,5)); // use maximum of position variances
    bool keep = true;
    if (tsp != lf->getLttWorldTrack()) { // always keep world-track
      //keep = keep && (mPVar < trkDelVarThresh);                       // 1) delete if high covariance
      keep = keep && (tspNew->getLastUpdateCounter() < trkDelUpCntThresh); // 2) delete if high noUpdateCounter
      keep = keep && (tspNew->getPointCount() > minPtCnt);                 // 3) delete if track too far away (far points are deleted above)
    }
    if (!keep) {
      cout << endl << "deleting track[" << tspNew->getUID() << "] with covar " << mPVar << ", last update at -" << tspNew->getLastUpdateCounter() << ", point count " << tspNew->getPointCount() << " and age " << tspNew->getAge() << flush;
      cf->removedTracks.push_back(tspNew);
    } else {
      tspNew->predict();
      cf->longTermTracks.push_back(tspNew);
      ++copiedTrackCount;
      // adjust stt[nbVStep] which is a multi-subset of long-term-tracks
      for (unsigned int i=0; i<lfStt[nbVStep].tracks.size(); ++i) {
        if (lfStt[nbVStep].tracks[i] == tsp)
          cfStt[nbVStep].tracks[i] = tspNew;
      }
      // adjust stt[nbVStep+1] which is a multi-subset of long-term-tracks
      for (unsigned int i=0; i<lfStt[nbVStep+1].tracks.size(); ++i) {
        if (lfStt[nbVStep+1].tracks[i] == tsp)
          cfStt[nbVStep+1].tracks[i] = tspNew;
      }
    }
  } // for each track in current buffer
  cout << endl << "copied " << copiedSTTrackCount << " short-term-tracks and " << copiedTrackCount << " out of " << totalTrackCount << " long-term-tracks from last frame into new frame" << flush;
}

DVector createInitState(DVector point3d, double egoVelocity) {
  DVector initState(12);
  initState(0) = 0.0; // yaw(RAD)
  initState(1) = 0.0; // pitch(RAD)
  initState(2) = 0.0; // roll(RAD)
  // Position: take coordinates of 1 point as Coordinate-System(CS) base, set orientation parallel to ego-CS
  DVectorRange(initState,ublas::range(3,6)) = point3d; // x,y,z
  // Speed: Assume zero ground-level-speed. As Orientation of CS is parallel to car and tracking is performed relative to ego-car, set only x-component
  initState(6) = 0.0; // roll-rate
  initState(7) = 0.0; // pitch-rate
  initState(8) = 0.0; // yaw-rate
  initState(9) = -egoVelocity; // x-velocity = -ego-velocity
  initState(10) = 0.0; // y-velocity
  initState(11) = 0.0; // z-velocity
  return initState;
}

double mlog(double x) {
  return log(max(1.0,1.0+x));
}

FeatureMaxRangeNormalizer* normalizerFromFile(std::string filename)
{
  FeatureMaxRangeNormalizer *normalizer = new FeatureMaxRangeNormalizer();
  try {
    ifstream infile(filename.c_str());
    string line;
    getline(infile,line);  // comment line
    getline(infile,line);
    istringstream iss(line);
    iss >> (*normalizer);
    cout << " loaded normalizer" << flush;
  } catch (std::exception &e) {
    delete normalizer;
    cout << " loading normalizer failed: " << e.what() << flush;
    normalizer = new FeatureMaxRangeNormalizer();
  }
  return normalizer;
}

void decideMerge(FrameSPtr frame, TrackRegistration::SPtr reg, bool coutOnly)
{
  ParameterHeap* params = ParameterHeap::get();
  const unsigned int nbVStep = frame->sttValidationSteps;
  vector<LidarFrame::ExtTrackList> &stt = frame->shortTermTracks;
  LidarFrame::ExtTrackList &buffer = stt[nbVStep-1]; // buffer to merge
  unsigned int tCount = buffer.tracks.size();
  const double fMaxDist = 2.0;
  const unsigned int fMaxAge = 20;

  if ((frame->longTermTracks.size() == 0) && (tCount == 1)) { // only one track -> world track -> always keep
    buffer.tracks[0]->mergeDecision = PointCloudTrack::MDKeep;
    return;
  }

  // classification / scoring parameters
  FeatureNormalizer *stage1normalizer = normalizerFromFile(params->stage1Normalizer);
  FeatureNormalizer *stage2normalizer = normalizerFromFile(params->stage2Hyperplane);
  Classifier        *stage1classifyer = NULL;
  double             stage2b;
  vector<double>     stage2w;
  if (!coutOnly) { // load parameters
    try {
      //stage1classifyer = new ExtraTreeClassifier("stage1_classifier");
      stage1classifyer = new SvmClassifier(params->stage1Classifier);
      cout << " loaded 1st-stage-classifier" << flush;
    } catch (std::exception &e) {
      cout << " loading 1st-stage-classifier failed: " << e.what() << flush;
      stage1classifyer = NULL;
    }
    // load 2nd-stage hyperplane parameters
    try {
      CsvReader csv(params->stage2Hyperplane,"\t ,;","#");
      stage2b = csv.get<double>(1, 0);
      for (unsigned int i=0; i<csv.fieldCount(2); ++i)
        stage2w.push_back(csv.get<double>(2, i));
      cout << " loaded 2nd-stage-hyperplane" << flush;
    } catch (const exception &e) {
      cout << " loading 2nd-stage-hyperplane failed: " << e.what() << flush;
      stage2b = 1.7;
      boost::array<double,1*37> user_w = { { 34.1296, -12.2061, -0.797605, -1.02433, -28.114, 2.17327, -0.616866, -0.0868315, -0.368506, 0, -0.341044, -5.23218, 0.517116, 0.0, -0.0766394, -0.152705, 0.797345, 1.48914, -1.06371, -0.181474, -0.0757976, -0.603907, -0.0195113, -8.34594, 0, -0.43991, 0.21694, 2.55696, 0, -8.89441, 0.502779, -0.352475, 0, -7.21229, 1.66212, 1.61145, 0.0 } };

      BOOST_FOREACH(double v, user_w) stage2w.push_back(v);
    }
  }
  const unsigned int NBF = stage2w.size();

  if (coutOnly)
    cout << endl << "informing about merge...";
  else
    cout << endl << "deciding about merge...";

  // prepare loop
  ofstream tdataout1;
  ofstream tdataout2;
  if (coutOnly) {
    tdataout1.open("tdata1.txt", ios::out | ios::app);
    tdataout1 << std::scientific;
    tdataout1 << std::setprecision(8);
    tdataout2.open("tdata2.txt", ios::out | ios::app);
    tdataout2 << std::scientific;
    tdataout2 << std::setprecision(8);
  }

  // MAIN LOOP: loop over all tracks to be decided on
  reg->usePixMask(false);
  for (unsigned int idx = 0; idx < tCount; ++idx) {
    cout << "." << flush;
    PointCloudTrack::SPtr track = buffer.tracks[idx];
    assert((track->getAge() >= 2) && "decideMerge(): track-age too small");
    assert((track->getAge()+1 == track->histAllStates.size()) && "decideMerge(): track-history-size doesn't match age");
    assert((fabs(ublas::norm_2(track->getState() - track->histAllStates.back().x)) <= 0.00001) && "decideMerge(): states don't match");
    DVector origstate6D = DVectorRange(track->histAllStates.front().x, ublas::range(0,6));
    DVector currentstate6D = DVectorRange(track->histAllStates.back().x, ublas::range(0,6));
    bool currentIcpSuccess = track->histAllStates.back().wasMeasured;
    BOOST_AUTO(hi0, track->histAllStates.rbegin());
    BOOST_AUTO(hi1, hi0++);
    BOOST_AUTO(hi2, hi1++);
    int nbIcpFails = 3 - (int)hi0->wasMeasured - (int)hi1->wasMeasured - (int)hi2->wasMeasured;
    DVector stateDiffDiff = hi0->x - 2*hi1->x + hi2->x; // (hi0 - hi1) - (hi1 - hi2)
    pair<double,double> errOrig = reg->calculateError(track, track->getState()); // Pt-Pl and Pt-Pt
    double errOrigPtPl = min(fMaxDist, errOrig.first);
    double errOrigPtPt = min(fMaxDist, errOrig.second);
    double moveInMostNormalDir;
    BOOST_AUTO(moveHist1, track->getMotion2NormalHistogram<4>(origstate6D, currentstate6D, &moveInMostNormalDir)); // ICP-measured move

    // pre-calc features for all links
    LidarFrame::TrackLinks links = frame->getMergedSTTLinks(nbVStep-1, idx, true); // include world track
    if (links.size() == 0)
      throw runtime_error("no links but should at least be linked to world");
    unsigned int totalLinkStrength = 0;
    unsigned int minLinkStrength = UINT_MAX;
    unsigned int maxLinkStrength = 0;
    bool oneIsMoving = false;
    vector<DVector> linkStateComp6D;
    vector<bool>   linkSuccfulICP; // whether the linked track was successfully registrated at the time this track was created (and linked)
    vector<double> linkErrPtPl;
    vector<double> linkErrPtPt;
    double minErrPtPl = DBL_MAX;
    double minErrPtPt = DBL_MAX;
    double maxErrPtPl = 0;
    double maxErrPtPt = 0;
    size_t minErrPtPlIdx = 0;
    size_t minErrPtPtIdx = 0;
    size_t maxLinkStrengthIdx = 0;
    for (int link = 0; link < (int)links.size(); ++link) {
      PointCloudTrack::SPtr  linkedTrack = links[link].track;
      unsigned int           linkStrength = links[link].linkStrength;
      assert((linkedTrack->getAge()+1 == linkedTrack->histAllStates.size()) && "render3D(): tc-history-size doesn't match age");
      assert((linkedTrack->getAge() > track->getAge()) && "render3D(): tc-history-size < track-history-size");
      //DVector ltOrigstate6D = DVectorRange(linkedTrack->histAllStates.front().x, ublas::range(0,6));
      //DVector ltCurrentstate6D = DVectorRange(linkedTrack->histAllStates.back().x, ublas::range(0,6));
      DVector tLinkedNow, t2linked;
      DMatrix RLinkedNow, R2linked;
      tLinkedNow = linkedTrack->histAllStates.back().t; // current state as R/t to transform to sensor-CS
      RLinkedNow = linkedTrack->histAllStates.back().R;
      unsigned int tAge = track->getAge();
      bool icpSucc = track->getRt2OtherTrkCS(R2linked, t2linked, *linkedTrack, tAge, tAge);
      DVector tcomp; DMatrix Rcomp;
      Rcomp = ublas::prod(RLinkedNow,R2linked);
      tcomp = ublas::prod(RLinkedNow,t2linked) + tLinkedNow;
      DVector xcomp6D; PointCloudTrack::Rt2state_T2E(Rcomp, tcomp, xcomp6D);
      // TODO: check: xcomp6D should be equal to tstate[-age]+(ltstate[curr]-ltstate[-age])
      // somehow does not seem to be fully equal, check!
//      DVector xcomp6Dalt =  track->histAllStates.front().x // creation state
//                           +linkedTrack->histAllStates[linkedTrack->histAllStates.size()-1].x // current state
//                           -linkedTrack->histAllStates[linkedTrack->histAllStates.size()-1-tAge].x; // state at "track"-creation
//      xcomp6Dalt.resize(6, true);
//      double xcompDiff = ublas::norm_2(xcomp6D-xcomp6Dalt);
//      if (xcompDiff > 1e-8) {
//        cerr << endl << "WARNING in decideMerge(): xcomp states for track [" << track->getUID() << "] are different:";
//        cerr << endl << "xcomp6D:    " << xcomp6D;
//        cerr << endl << "xcomp6Dalt: " << xcomp6Dalt;
//        cerr << flush;
//      }

      totalLinkStrength += linkStrength;
      minLinkStrength = min(minLinkStrength, linkStrength);
      maxLinkStrength = max(maxLinkStrength, linkStrength);
      if (linkStrength == maxLinkStrength) maxLinkStrengthIdx = link;
      if (linkedTrack->isMovingObject()) oneIsMoving = true;
      linkStateComp6D.push_back(xcomp6D);
      linkSuccfulICP.push_back(icpSucc);
      pair<double,double> errLink = reg->calculateError(track, xcomp6D); // Pt-Pl and Pt-Pt
      errLink.first  = min(fMaxDist, errLink.first);
      errLink.second = min(fMaxDist, errLink.second);
      linkErrPtPl.push_back(errLink.first);
      linkErrPtPt.push_back(errLink.second);
      minErrPtPl = min(minErrPtPl, errLink.first);
      maxErrPtPl = max(maxErrPtPl, errLink.first);
      minErrPtPt = min(minErrPtPt, errLink.second);
      maxErrPtPt = max(maxErrPtPt, errLink.second);
      if (minErrPtPl == errLink.first) minErrPtPlIdx = link;
      if (minErrPtPt == errLink.second) minErrPtPtIdx = link;
    }

    // gather features for first-stage classification (ignore|keen|merge)
    //vector<double> fVec1;
    FeatureVector fVec1;
    FeatureVector fVec1New;
    // [0]
    fVec1.push_back(currentIcpSuccess);
    fVec1.push_back(nbIcpFails);
    fVec1.push_back(mlog(norm_2(stateDiffDiff))); // characterizes if the two last relative movements were approximately the same
    fVec1.push_back(errOrigPtPl);
    fVec1.push_back(errOrigPtPt);
    fVec1.push_back(track->getLastUpdateCounter());
    fVec1.push_back(mlog(track->getPointCount()));
    fVec1.push_back(moveInMostNormalDir);
    fVec1.push_back(mlog(moveHist1[0])); // number of surfaces parallel to move ZONK!!! should have selected [3] = perpendicular to move
    fVec1.push_back(mlog(moveHist1[0]+moveHist1[1]));
    // [10]
    fVec1.push_back((moveHist1[0])/(double)(track->getPointCount()));
    fVec1.push_back((moveHist1[0]+moveHist1[1])/(double)(track->getPointCount()));
    fVec1.push_back(links.size()); // number of links
    fVec1.push_back(mlog(totalLinkStrength)); // summed linkStrength

    // [0]
    fVec1New.push_back(currentIcpSuccess);
    fVec1New.push_back(nbIcpFails);
    fVec1New.push_back(mlog(norm_2(stateDiffDiff))); // characterizes if the two last relative movements were approximately the same
    fVec1New.push_back(errOrigPtPl);
    fVec1New.push_back(errOrigPtPt);
    fVec1New.push_back(track->getLastUpdateCounter());
    fVec1New.push_back(mlog(track->getPointCount()));
    fVec1New.push_back(moveInMostNormalDir);
    fVec1New.push_back(mlog(moveHist1[3]));
    fVec1New.push_back(mlog(moveHist1[2]));
    // [10]
    fVec1New.push_back(mlog(moveHist1[1]));
    fVec1New.push_back(mlog(moveHist1[3]+moveHist1[2]));
    fVec1New.push_back((moveHist1[3])/(double)(track->getPointCount()));
    fVec1New.push_back((moveHist1[3]+moveHist1[2])/(double)(track->getPointCount()));
    fVec1New.push_back(links.size()); // number of links
    fVec1New.push_back(mlog(totalLinkStrength)); // summed linkStrength

    if (links.size() == 0) {
      for (size_t i=0; i<3*12; ++i) {
        fVec1.push_back(0);
        fVec1New.push_back(0);
      }
    } else {
      boost::array<size_t, 3> indices = {{ minErrPtPlIdx, minErrPtPtIdx, maxLinkStrengthIdx }};
      BOOST_FOREACH(size_t idx, indices) {
        // [14]/[26]/[38]
        fVec1.push_back(linkSuccfulICP[idx]);
        fVec1.push_back(mlog(links[idx].linkStrength));
        fVec1.push_back(links[idx].linkStrength/max(1.,(double)totalLinkStrength)); // avoid div by 0
        fVec1.push_back(links[idx].linkStrength/max(1.,(double)maxLinkStrength)); // avoid div by 0
        fVec1.push_back(links[idx].track->getLastUpdateCounter());
        fVec1.push_back(min(fMaxAge,links[idx].track->getAge()));
        fVec1.push_back(mlog(links[idx].track->getPointCount()));
        fVec1.push_back(mlog(links[idx].track->getPointCount()/(double)track->getPointCount())); // point count is never 0
        fVec1.push_back(linkErrPtPl[idx]);
        fVec1.push_back(linkErrPtPt[idx]);
        fVec1.push_back(mlog(linkErrPtPl[idx] / max(1e-9,errOrigPtPl)));
        fVec1.push_back(mlog(linkErrPtPt[idx] / max(1e-9,errOrigPtPt)));
        // ^^^ [25]/[37]/[49]

        // [16]/[28]/[40]
        fVec1New.push_back(linkSuccfulICP[idx]);
        fVec1New.push_back(mlog(links[idx].linkStrength));
        fVec1New.push_back(links[idx].linkStrength/max(1.,(double)totalLinkStrength)); // avoid div by 0
        fVec1New.push_back(links[idx].linkStrength/max(1.,(double)maxLinkStrength)); // avoid div by 0
        fVec1New.push_back(links[idx].track->getLastUpdateCounter());
        fVec1New.push_back(min(fMaxAge,links[idx].track->getAge()));
        fVec1New.push_back(mlog(links[idx].track->getPointCount()));
        fVec1New.push_back(mlog(links[idx].track->getPointCount()/(double)track->getPointCount())); // point count is never 0
        fVec1New.push_back(linkErrPtPl[idx]);
        fVec1New.push_back(linkErrPtPt[idx]);
        fVec1New.push_back(mlog(linkErrPtPl[idx] / max(1e-9,errOrigPtPl)));
        fVec1New.push_back(mlog(linkErrPtPt[idx] / max(1e-9,errOrigPtPt)));
        // ^^^ [27]/[39]/[51]
      }
    }

    // carry out classification or write log file
    if (coutOnly) { // output result
      tdataout1 << endl;
      for (size_t i=0; i<fVec1New.size(); ++i)
        tdataout1 << fVec1New[i] << "\t";
      switch (track->mergeDecision) {
        case PointCloudTrack::MDIgnore: tdataout1 << "-->\t" << "1"; break;
        case PointCloudTrack::MDKeep:   tdataout1 << "-->\t" << "2"; break;
        case PointCloudTrack::MDMerge:  tdataout1 << "-->\t" << "3"; break;
      }
    } else {
      if (params->trkFeatureVersion == 2)
        fVec1 = fVec1New;
      stage1normalizer->normalize(fVec1);
      if (stage1classifyer) {
        switch (stage1classifyer->classify(fVec1)) {
          case 1: track->mergeDecision = PointCloudTrack::MDIgnore; break;
          case 2: track->mergeDecision = PointCloudTrack::MDKeep; break;
          case 3: track->mergeDecision = PointCloudTrack::MDMerge; break;
        }
      } else {
        track->mergeDecision = PointCloudTrack::MDIgnore;
      }
      if (ParameterHeap::get()->disableTracking) { // only mapping => always merge
        track->mergeDecision = PointCloudTrack::MDMerge;
      }
    }
    if (track->mergeDecision != PointCloudTrack::MDKeep) {
      // gather features for all links and decide
      double maxScore = -DBL_MAX;
      PointCloudTrack::SPtr maxLink;
      vector<double> scoreBuff(links.size(),-DBL_MAX); // TODO(9): can be removed when everything is working
      for (int link = 0; link < (int)links.size(); ++link) {
        PointCloudTrack::SPtr  linkedTrack = links[link].track;
        unsigned int           linkStrength = links[link].linkStrength;

        // calculate move of "tc" within the same period of time, thus -p0+pcurr
        DVector xcomp6D = linkStateComp6D[link];
        DVector stateDiff = currentstate6D-xcomp6D;
        DSMatrix tCovar = DSMatrixRangeConst(track->getCovar(), ublas::range(0,6), ublas::range(0,6));
        DCMatrix tCovarI = invSym(tCovar);
        DSMatrix ltCovar = DSMatrixRangeConst(linkedTrack->getCovar(), ublas::range(0,6), ublas::range(0,6));
        DCMatrix ltCovarI = invSym(ltCovar);
        double stateDiffL2Norm = ublas::norm_2(stateDiff);
        double stateDiffMahaTNorm = ublas::inner_prod(stateDiff,ublas::prod(tCovarI,stateDiff));
        double stateDiffMahaLTNorm = ublas::inner_prod(stateDiff,ublas::prod(ltCovarI,stateDiff));
        // histograms[0] -> move perpendicular to normal
  //      BOOST_AUTO(moveHist2, track->getMotion2NormalHistogram<4>(origstate6D, xcomp6D)); // applied move of linked track
  //      BOOST_AUTO(moveHist3, track->getMotion2NormalHistogram<4>(currentstate6D, xcomp6D)); // difference

        FeatureVector fVec2;
        FeatureVector fVec2New;
        // [0]
        fVec2.push_back(links.size()); // number of links
        fVec2.push_back(mlog(linkStrength)); // 0 if world track
        fVec2.push_back(mlog(totalLinkStrength)); // summed linkStrength
        fVec2.push_back(linkStrength/max(1.,(double)totalLinkStrength));
        fVec2.push_back(linkStrength/max(1.,(double)maxLinkStrength));

        fVec2New.push_back(links.size()); // number of links
        fVec2New.push_back(mlog(linkStrength)); // 0 if world track
        fVec2New.push_back(mlog(totalLinkStrength)); // summed linkStrength
        fVec2New.push_back(linkStrength/max(1.,(double)totalLinkStrength));
        fVec2New.push_back(linkStrength/max(1.,(double)maxLinkStrength));

        // [5]
        fVec2.push_back(errOrigPtPl);
        fVec2.push_back(linkErrPtPl[link]);
        fVec2.push_back(mlog(linkErrPtPl[link] / max(1e-9,errOrigPtPl)));
        fVec2.push_back(mlog(linkErrPtPl[link] / max(1e-9,minErrPtPl)));
        fVec2.push_back(errOrigPtPt);
        fVec2.push_back(linkErrPtPt[link]);
        fVec2.push_back(mlog(linkErrPtPt[link] / max(1e-9,errOrigPtPt)));
        fVec2.push_back(mlog(linkErrPtPt[link] / max(1e-9,minErrPtPt)));

        fVec2New.push_back(errOrigPtPl);
        fVec2New.push_back(linkErrPtPl[link]);
        fVec2New.push_back(mlog(linkErrPtPl[link] / max(1e-9,errOrigPtPl)));
        fVec2New.push_back(mlog(linkErrPtPl[link] / max(1e-9,minErrPtPl)));
        fVec2New.push_back(errOrigPtPt);
        fVec2New.push_back(linkErrPtPt[link]);
        fVec2New.push_back(mlog(linkErrPtPt[link] / max(1e-9,errOrigPtPt)));
        fVec2New.push_back(mlog(linkErrPtPt[link] / max(1e-9,minErrPtPt)));

        // [13]
        fVec2.push_back(track->getLastUpdateCounter());
        fVec2.push_back(linkedTrack->getLastUpdateCounter());
        fVec2.push_back(mlog(track->getPointCount()));
        fVec2.push_back(mlog(linkedTrack->getPointCount()));
        fVec2.push_back(mlog(track->getPointCount() / max(1.,(double)linkedTrack->getPointCount())));

        fVec2New.push_back(track->getLastUpdateCounter());
        fVec2New.push_back(linkedTrack->getLastUpdateCounter());
        fVec2New.push_back(mlog(track->getPointCount()));
        fVec2New.push_back(mlog(linkedTrack->getPointCount()));
        fVec2New.push_back(mlog(track->getPointCount() / max(1.,(double)linkedTrack->getPointCount())));

        // [18]
        fVec2.push_back(linkedTrack->isMovingObject());
        fVec2.push_back(oneIsMoving); // 1 if one of the links is a moving object
        fVec2.push_back(mlog(stateDiffL2Norm));
        fVec2.push_back(mlog(stateDiffMahaTNorm));
        fVec2.push_back(mlog(stateDiffMahaLTNorm));

        fVec2New.push_back(linkedTrack->isMovingObject());
        fVec2New.push_back(oneIsMoving); // 1 if one of the links is a moving object
        fVec2New.push_back(mlog(stateDiffL2Norm));
        fVec2New.push_back(mlog(stateDiffMahaTNorm));
        fVec2New.push_back(mlog(stateDiffMahaLTNorm));

        // [23]
        double hist1val1 = (moveHist1[3]+moveHist1[2])==0 ? 0.0 : ((double)(moveHist1[1]+moveHist1[0]))/((double)(moveHist1[3]+moveHist1[2])); // relative number of surfaces perpendicular to move
        double hist1val2 = (moveHist1[3]==0) ? 0.0 : ((double)(moveHist1[2]+moveHist1[1]+moveHist1[0]))/((double)(moveHist1[3])); // relative number of surfaces perpendicular to move
        fVec2.push_back(mlog(moveHist1[0])); // number of surfaces perpendicular to move
        fVec2.push_back(mlog(moveHist1[0]+moveHist1[1]));
        fVec2.push_back(hist1val1);
        fVec2.push_back(hist1val2);
        fVec2.push_back((double)(moveHist1[0])/(double)(max(1,(int)track->getPointCount())));
        fVec2.push_back((double)(moveHist1[0]+moveHist1[2])/(double)(max(1,(int)track->getPointCount())));
        fVec2.push_back(moveInMostNormalDir);
        // ^^^ [29]

        // [23]
        double hist1val1New = (moveHist1[0]+moveHist1[1])==0 ? 0.0 : ((double)(moveHist1[2]+moveHist1[3]))/((double)(moveHist1[0]+moveHist1[1])); // relative number of surfaces perpendicular to move
        double hist1val2New = (moveHist1[0]==0) ? 0.0 : ((double)(moveHist1[1]+moveHist1[2]+moveHist1[3]))/((double)(moveHist1[0])); // relative number of surfaces perpendicular to move
        fVec2New.push_back(mlog(moveHist1[3])); // number of surfaces perpendicular to move
        fVec2New.push_back(mlog(moveHist1[2]));
        fVec2New.push_back(mlog(moveHist1[1]));
        fVec2New.push_back(mlog(moveHist1[3]+moveHist1[2]));
        fVec2New.push_back(hist1val1New);
        fVec2New.push_back(hist1val2New);
        fVec2New.push_back((double)(moveHist1[3])/(double)(max(1,(int)track->getPointCount())));
        fVec2New.push_back((double)(moveHist1[3]+moveHist1[2])/(double)(max(1,(int)track->getPointCount())));
        fVec2New.push_back(moveInMostNormalDir);
        // ^^^ [31]

        // decide upon feature vector or write log file
        if (coutOnly) {
          if (links.size() > 1) {
            tdataout2 << endl;
            for (size_t i=0; i<fVec2New.size(); ++i) {
              tdataout2 << fVec2New[i] << "\t";
            }
            int decision = (track->trackMaxMergeScore == linkedTrack) ? 1 : -1;
            tdataout2 << "-->\t" << decision;
          }
        } else {
          if (params->trkFeatureVersion == 2)
            fVec2 = fVec2New;
          assert(fVec2.size() == NBF && "SVM-dimensionality != feature dimensionality!!!");
          stage2normalizer->normalize(fVec2);
          double score = stage2b; // positive --> merge, negative --> don't merge
          for (unsigned int dim=0; dim<NBF; ++dim) {
            score += stage2w[dim] * fVec2[dim];
          }
          if ((links.size() == 1) && (track->getLastUpdateCounter() > 0))
            score = 1.0; // always link(=merge if no 1st stage), if only 1 link and last ICP step did not succeed
          scoreBuff[link] = score;
          if (score <= -DBL_MAX)
            throw runtime_error("score is -DBL_MAX");
          if ((score > maxScore) && ((score > 0.0) || (stage1classifyer))) { // merge if positive score. the link with "maxScore" wins
            maxScore = score;
            maxLink = linkedTrack;
          } else {
            cout << " ";
          }
        } // end: no output but decide
      } // end: loop over links
      if (!coutOnly) {
        track->trackMaxMergeScore = maxLink; // might be NULL (e.g. if decision==keep)
        if (!stage1classifyer) {
          // try to decide upon scores, if 1st-stage classifier not loaded
          if (maxLink.get() == NULL) {
            track->mergeDecision = PointCloudTrack::MDKeep;
          } else {
            if (track->getMergeTrackIcpSuccess())
              track->mergeDecision = PointCloudTrack::MDMerge;
            else
              track->mergeDecision = PointCloudTrack::MDIgnore;
          }
        } // end: first stage classifier not valid
      } // end: !coutOnly --> decide
      if ((track->mergeDecision != PointCloudTrack::MDKeep) && (track->trackMaxMergeScore.get() == NULL))
        throw runtime_error("link must be valid if track is not kept");
    } // end: second stage decision about linked track
    if ((track->mergeDecision != PointCloudTrack::MDKeep) && (track->trackMaxMergeScore.get() == NULL))
      throw runtime_error("link must be valid if track is not kept");
  } // end: loop over buffer.tracks
  reg->usePixMask(true);

  delete stage1classifyer;
  delete stage2normalizer;
  delete stage1normalizer;
  cout << "done" << flush;
}

bool trackAngleSort(const PointCloudTrack::SPtr ti, const PointCloudTrack::SPtr tj) {
  return atan2(ti->getState()[4], ti->getState()[3]) < atan2(tj->getState()[4], tj->getState()[3]);
}

PointCloudTrack::SPtr getTrack(const LidarFrame::TrackSPtrList &tracklist, const size_t idx) {
  size_t ci = 0;
  for (BOOST_AUTO(cti, tracklist.begin()); cti != tracklist.end(); ++cti, ++ci) {
    if (ci == idx)
      return *cti;
  }
  return PointCloudTrack::SPtr();
}

void mergeSplitCreateTracks(FrameSPtr frame, TrackRegistration::SPtr reg)
{
  // TODO (8): early-merge track as soon as ICP fails (speedup). may break projection-association??
  ParameterHeap* params = ParameterHeap::get();
  DSMatrix &initVariance = params->trkInitVariance;
  DSMatrix &predictionVariance = params->trkPredVariance;
  const bool useDynamicMapping = params->useDynamicMapping; // whether to add appearance (prerequisite for next flag)
  const bool useSimpleMapping = params->useSimpleMapping; // add appearance without checking neighbors
  const bool useOverwriteMapping = params->useOverwriteMapping; // replace appearance
  //const double projDistTolerance = params->projectionDistTolerance;
  const unsigned int nbVStep = frame->sttValidationSteps;
  vector<LidarFrame::ExtTrackList> &stt = frame->shortTermTracks;
  LidarFrame::TrackSPtrList &ltt = frame->longTermTracks;
  int hsize, vsize;
  frame->segments.getSize(hsize,vsize);

  if (stt[0].tracks.size() > 0) {
    if (stt[0].tracks[0]->getAge() == 0) {
      cout << endl << "already created new tracks (age of first track = 0). skipping." << flush;
      return;
    }
  }

  // 1) process all existing tracks in shortTermTracks[nbVStep-1]
  //    and decide about merge/split with shortTermTracks[nbVStep]
  //    which should be equal to / contained in longTermTracks
  decideMerge(frame, reg, true); // TODO (5): switch from output of training data to deciding once it works
  unsigned int merged = 0;
  unsigned int split = 0;
  for (unsigned int idx = 0; idx < stt[nbVStep-1].tracks.size(); ++idx) {
    PointCloudTrack::SPtr track = stt[nbVStep-1].tracks[idx];
    PointCloudTrack::SPtr trackToMergeInto = track->trackMaxMergeScore; // determined in function above
    switch (track->mergeDecision) {
      case PointCloudTrack::MDIgnore:
        // simply ignore/delete this track. but overwrite link, so tracks from the next frame don't loose their potential links
        stt[nbVStep-1].tracks[idx] = trackToMergeInto;
        break;
      case PointCloudTrack::MDMerge:
        // merge this track into linked track (should be long-term-track)
        trackToMergeInto->merge(track.get(), useOverwriteMapping, (useDynamicMapping) || (trackToMergeInto->getUID() == 1), !useSimpleMapping, track->getAge());
        ++merged;
        // after merge/ignore, replace SPtr in track-list by linked track
        stt[nbVStep-1].tracks[idx] = trackToMergeInto;
        break;
      case PointCloudTrack::MDKeep:
        // leave track as it is, but create a copy in longTermTracks.
        if (ltt.size() != 0) // only for moving objects
          track->isMovingObject(true);
        ltt.push_back(track);
        // however, remove corresponding 3D points in linked tracks
        PointCloudTrack::PointIterator ptbegin, ptend;
        track->getPoints(ptbegin, ptend);
        if (stt[nbVStep-1].linksT1.get() != NULL) { // only if links are present at all
          LidarFrame::TrackLinks links = frame->getMergedSTTLinks(nbVStep-1, idx); // get all links of track, without adding world explicitly
          // TODO (8): merging linkage graph could be carried out only once (modify graph)
          BOOST_FOREACH(LidarFrame::TrackLink &l, links) {
            PointCloudTrack::SPtr ltrack = l.track;
            DMatrix R2l;
            DVector t2l,p;
            track->getRt2OtherTrkCS(R2l, t2l, *ltrack, track->getAge(), track->getAge()); // at time when track was created
            cout << endl << "removing in track [" << ltrack->getUID() << "] with R/t=" << R2l << t2l << flush;
            for (PointCloudTrack::PointIterator pti=ptbegin; pti!=ptend; ++pti) {
              p = ublas::prod(R2l,pti->s->position) + t2l;
              //cout << endl << pti->s->position << " -> " << p << flush;
              ltrack->removeNeighborhoodTrkCS(p, 1); // remove in 3x3x3 neighborhood
            } // end: for each point in linked track
          } // end: for each linked track
        } // end: links exist
        ++split;
        break;
    }
  } // end: loop over all tracks in shortTermTracks[STTSize-2]
  cout << "created " << split << " new ltt-tracks, merged " << merged << "..." << flush;

  // 2) tracks in shortTermTracks[nbVStep-1] have now been merged into longTermTracks
  //    with the help of shortTermTracks[nbVStep] which can be emptied and all tracks
  //    moved from trackBuffer[i] into trackBuffer[i+1]
  for (int bi=stt.size()-2; bi>=0; --bi) {
    stt[bi].tracks.swap(stt[bi+1].tracks);
    stt[bi].linksT1.swap(stt[bi+1].linksT1);
    stt[bi].linksT2.swap(stt[bi+1].linksT2);
    stt[bi].lttLinks.swap(stt[bi+1].lttLinks);
  }
  if (frame->projSttIdx < UINT_MAX) frame->projSttIdx++; // shift makes corresponding projection-buffer refer to next buffer index
  // [nbVStep+1] was swapped through to [0] --> reset
  stt[0].tracks.clear();
  stt[0].linksT1.reset();
  stt[0].linksT2.reset();
  stt[0].lttLinks.reset();
  stt[nbVStep].linksT2.reset(); // the track list cannot have valid linkages (which would be in next non-existent buffer)
  stt[nbVStep+1].linksT1.reset(); // the track list cannot have valid linkages (which would be in next non-existent buffer)
  stt[nbVStep+1].linksT2.reset(); // the track list cannot have valid linkages (which would be in next non-existent buffer)

  // 3) create new tracks (from segmentation into shortTermTracks[0])
  bool lttExists = !(ltt.empty() && stt[1].tracks.empty()); // no long-term tracks exist and no tracks in queue
  if ((!lttExists) || (params->disableTracking)) {
    // 3-A) create one world-track only, ignoring segmentation
    ColRowPair cr = frame->distance.firstValidCR(DBL_MAX);
    ColRowPair crEnd = frame->distance.endCR();
    DVector initState = createInitState(frame->point3D.get(cr), params->disableLocalization ? 0.0 : frame->speed);
    PointCloudTrack::SPtr newTrack(new PointCloudTrack(initState, initVariance, predictionVariance, hsize/2, vsize/2));
    if (!lttExists) {
      newTrack->setColor(0.3, 0.3, 0.3); // world-track = gray
      newTrack->isMovingObject(false);
    }
    while (cr != crEnd) {
      newTrack->addPoint(frame->point3D.get(cr),frame->normal3D.get(cr),frame->normalConfidence.get(cr),frame->normalStdDevRAD.get(cr),
          !useSimpleMapping,frame->pointVariance3D.get(cr),frame->normalVariance3D.get(cr));
      cr = frame->distance.nextValidCR(cr, DBL_MAX);
    }
    stt[0].tracks.clear();
    stt[0].tracks.push_back(newTrack);
    cout << "done. created 1 world-track with age " << newTrack->getAge() << "." << flush;
  } else {
    // 3-B-1) create one track from each segment into shortTermTracks[0]
    typedef map< LidarSegment*, LidarFrame::TrackIndex > Seg2TrackIdxMap;
    Seg2TrackIdxMap segTrackMap; // store segment->track correspondence in segTrackMap (used in 3-B-2/3)
    for (int row = 0; row < vsize; ++row) {
      for (int col = 0; col < hsize; ++col) {
        LidarSegment &seg = frame->segments.get(col,row);
        if ((!seg.isValid()) || (seg.isProxy()))
          continue;
        LidarSegment::ColRowIterator cr, crEnd; unsigned int pixCnt;
        float colC, rowC;
        seg.getColRows(cr, crEnd, pixCnt);
        seg.getCenter(colC, rowC);
        DVector initState = createInitState(frame->point3D.get(*cr), frame->speed);
        PointCloudTrack::SPtr newTrack(new PointCloudTrack(initState, initVariance, predictionVariance, colC, rowC));
        newTrack->setColor(seg.getColorR(), seg.getColorG(), seg.getColorB());
        // add points
        while (cr != crEnd) {
          newTrack->addPoint(frame->point3D.get(*cr),frame->normal3D.get(*cr),frame->normalConfidence.get(*cr),frame->normalStdDevRAD.get(*cr),
              !useSimpleMapping,frame->pointVariance3D.get(*cr),frame->normalVariance3D.get(*cr));
          ++cr;
        }
        // insert into buffer
        cout << " c[" << newTrack->getUID() << "]:0+" << seg.getPxCount();
        segTrackMap[seg.getPtr()] = stt[0].tracks.size();
        stt[0].tracks.push_back(newTrack);
      } // loop over col
    } // loop over row
    // 3-B-2) convenient for visualizing track by track: sort tracks by angle:
    std::vector<PointCloudTrack::SPtr> stt0tracks_copy = stt[0].tracks; // copy shared-pointers in temporary buffer
    sort(stt0tracks_copy.begin(), stt0tracks_copy.end(), trackAngleSort); // sort them
    for (Seg2TrackIdxMap::iterator stii=segTrackMap.begin(); stii!=segTrackMap.end(); ++stii) { // adjust mapping
      PointCloudTrack *track = stt[0].tracks[stii->second].get();
      assert(track != NULL && "mergeSplitCreateTracks: track-pointer is NULL");
      for (unsigned int i=0; i<stt0tracks_copy.size(); ++i) { // search new position
        if (stt0tracks_copy[i].get() == track) {
          stii->second = i; // store new position
          break; // stop search, continue with next mapping entry
        }
      }
    }
    stt[0].tracks = stt0tracks_copy; // apply sort to real buffer
    cout << "..." << flush;
    // 3-B-3) associate created tracks with projected tracks from last time step
    if (frame->projSttIdx != 1) {
      cout << endl << "WARNING: segmentsToTracks: frame->projSttIdx == " << frame->projSttIdx << " != 1 " << flush;
    } else {
      stt[0].linksT1 = LidarFrame::TrackLinkageGraphSPtr(new LidarFrame::TrackLinkageGraph());
      stt[0].linksT2 = LidarFrame::TrackLinkageGraphSPtr(new LidarFrame::TrackLinkageGraph());
      stt[0].lttLinks = LidarFrame::LTTrackLinkageGraphSPtr(new LidarFrame::LTTrackLinkageGraph());
      unsigned int addedLinksT1 = 0;
      unsigned int addedLinksT2 = 0;
      unsigned int addedLinksTL = 0;
      unsigned int validSegs = 0;
      unsigned int validProj = 0;
      for (int row = 0; row < vsize; ++row) {
        for (int col = 0; col < hsize; ++col) {
          LidarSegment *seg = frame->segments.get(col,row).getPtr();
          LidarFrame::TrackIndex trkIdx = seg->isValid() ? segTrackMap[seg] : UINT_MAX;
          LidarFrame::TrackIndex trk1Idx = frame->extProjTracks.get(col,row);
          LidarFrame::TrackIndex trk2Idx = frame->projLCTracks.get(col,row);
          LidarFrame::TrackIndex lttTrkIdx = frame->projLttTracks.get(col,row);
          PointCloudTrack::IdT   lttTrkID = 0;
          if (lttTrkIdx != UINT_MAX) lttTrkID = getTrack(frame->longTermTracks, lttTrkIdx)->getUID();

          if (trkIdx != UINT_MAX) ++validSegs;
          if (trk1Idx != UINT_MAX) ++validProj;
          // add to linksT1
          if ((trkIdx != UINT_MAX) && (trk1Idx != UINT_MAX)) {
            stt[0].linksT1->addConnection(trkIdx, trk1Idx);
            ++addedLinksT1;
//          } else {
//            if (trkIdx != UINT_MAX) stt[0].linksT1->addNode1(trkIdx);
//            if (trk1Idx != UINT_MAX) stt[0].linksT1->addNode2(trk1Idx);
          }
          // add to linksT2
          if ((trkIdx != UINT_MAX) && (trk2Idx != UINT_MAX)) {
            stt[0].linksT2->addConnection(trkIdx, trk2Idx);
            ++addedLinksT2;
//          } else {
//            if (trkIdx != UINT_MAX) stt[0].linksT2->addNode1(trkIdx);
//            if (trk2Idx != UINT_MAX) stt[0].linksT2->addNode2(trk2Idx);
          }
          // add to lttLinks
          if ((trkIdx != UINT_MAX) && (lttTrkIdx != UINT_MAX)) {
            stt[0].lttLinks->addConnection(trkIdx, lttTrkID);
            ++addedLinksTL;
//          } else {
//            if (trkIdx != UINT_MAX) stt[0].lttLinks->addNode1(trkIdx);
//            if (lttTrkIdx != UINT_MAX) stt[0].lttLinks->addNode2(lttTrkID);
          }
        }
      }
      cout << "done. created " << stt[0].tracks.size() << " new tracks. linked " << validSegs << " segment-pix with " << validProj << " projection-pix => " << addedLinksT1 << "+" << addedLinksT2 << "+" << addedLinksTL << " links" << flush;
    }
  } // end: create short-term-tracks from segments

  //cout << "done. mapped " << dynMappingCount << " segments, created " << newTrackCount << " new tracks resulting in " << frame->trackList.size() << " total tracks." << flush;

}

void outputSurface(ofstream &outfile, PointCloudTrack::Surface::SPtr s, const mdefs::DMatrix &R, const mdefs::DVector &t, mdefs::DVector &p)
{
  // points, normal, NC
  p = ublas::prod(R, s->position) + t;
  for (size_t i=0;i<3; ++i) outfile << " " << p[i];
  p = ublas::prod(R, s->normal);
  for (size_t i=0;i<3; ++i) outfile << " " << p[i];
  outfile << " " << s->normalConfidence;
}

// output track details in a text-file with same name as current input-image
void outputLTTTrackDetails(FrameSPtr frame, bool skipWorld)
{
  string outsubdir = ParameterHeap::get()->outDir;
  if (outsubdir == "")
    return;

  // determine file name automatically
  string extension = ".txt";
  if (!is_directory(outsubdir))
    create_directory(outsubdir);
  ofstream outfile;

  path imgSource = path(outsubdir) / path("sourcedir");
  outfile.open(imgSource.string().c_str(), ios_base::trunc);
  outfile << path(frame->sourcefile).parent_path() << endl;
  outfile.close();

  path outfilename = path(outsubdir) / path(path(frame->sourcefile).filename().string() + extension);
  if (exists(outfilename))
    remove(outfilename); // delete old file
  cout << endl << "writing tracks into " << outfilename.string() << " " << flush;

  // output tracks
  mdefs::DVector x,t,p;
  mdefs::DMatrix R, covar;
  PointCloudTrack::PointIterator ptI, ptFirst, ptEnd;
  size_t ptCnt, vipCnt;
  outfile.open(outfilename.string().c_str());
  BOOST_FOREACH(PointCloudTrack::SPtr track, frame->longTermTracks) {
    if ((skipWorld) && (track == frame->getLttWorldTrack()))
      continue;
    track->getState(x); // [rx,ry,rz,tx,ty,tz,rxDot,ryDot,rzDot,txDot,tyDot,tzDot]
    track->getCovar(covar);
    track->getRt2EgoCS(R, t); // transform to Ego-CS: p'=R*p+t
    track->getPoints(ptFirst, ptEnd);
    ptCnt=0; vipCnt=0;
    for (ptI = ptFirst; ptI != ptEnd; ++ptI) {
      ptCnt++;
      if (ptI->s->isVIP)
        vipCnt++;
    }
    outfile << "[" << track->getUID() << "]";
    outfile << " " << track->getAge();
    outfile << " " << track->getLastUpdateCounter();
    outfile << " " << track->getColorR();
    outfile << " " << track->getColorG();
    outfile << " " << track->getColorB();
    BOOST_FOREACH(double d, x) // state-vector
      outfile << " " << d;
    for (size_t row=0;row<covar.size1(); ++row) // covariance
      for (size_t col=0;col<covar.size2(); ++col)
        outfile << " " << covar(row,col);
    // output point cloud: nbVIPPts, totalNbPts, surfaces[VIP first]
    outfile << " " << vipCnt;
    outfile << " " << ptCnt;
    for (ptI = ptFirst; ptI != ptEnd; ++ptI) {
      if (ptI->s->isVIP)
        outputSurface(outfile, ptI->s, R, t, p);
    }
    for (ptI = ptFirst; ptI != ptEnd; ++ptI) {
      if (!ptI->s->isVIP)
        outputSurface(outfile, ptI->s, R, t, p);
    }
    outfile << endl;
  }
}
