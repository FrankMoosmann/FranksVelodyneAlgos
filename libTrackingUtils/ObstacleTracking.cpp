#include "ObstacleTracking.hpp"

#include <cmath>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO

#include "ParameterHeap.hpp"
#include "AlgoFrameFeatures.hpp"
#include "AlgoSegmentation.hpp"
#include "AlgoRegistration.hpp"
#include "AlgoTrackGeneration.hpp"

using namespace std;

ObstacleTracking::ObstacleTracking(const unsigned int buffersize, const LidarImageProjector &imageProj)
  : imageProjector(imageProj)
   ,frameBuffer(buffersize, imageProjector.getImgHorizSize(), imageProjector.getImgVertSize())
{
}

ObstacleTracking::~ObstacleTracking()
{
  BOOST_FOREACH(PointCloudTrack::SPtr tsp, frameBuffer.get(0)->longTermTracks)
    tsp->removeFarPoints(0.0); // write points of world-track to file
}

FrameSPtr ObstacleTracking::getFrame(int relPastFrame)
{
  return frameBuffer.get(relPastFrame);
}

void ObstacleTracking::next()
{
  frameBuffer.next();
  registrator.reset();
}

void ObstacleTracking::reset(const unsigned int buffersize)
{
  frameBuffer = LidarFrameBuffer(buffersize, imageProjector.getImgHorizSize(), imageProjector.getImgVertSize());
}

void ObstacleTracking::processFrame(unsigned int relPastFrame)
{
  calcFeatures(relPastFrame);
  calcRegistration(relPastFrame);
//  calcEgoMotion(relPastFrame);
//  calcTrackEgoMotionUpdate(relPastFrame);
  calcSegmentation(relPastFrame);
  calcTracks(relPastFrame);
}

void ObstacleTracking::calcFeatures(unsigned int relPastFrame)
{
  try {
    FrameSPtr f = frameBuffer.get(relPastFrame);
    cerr << endl << "calculating features..." << flush;
    postprocessDistanceData(f, imageProjector);
    calculateDerivative(f);
    calculateConnections(f);
    calculateNormals(f, imageProjector);
    calculateMatchingFeatures(f); // only calculates the features, does not match between frames!
    FrameSPtr lf = frameBuffer.get(relPastFrame+1);
    calculateDistDiff(lf, f);
  } catch (std::range_error e) {
    // do nothing
  }
  cerr << "done" << flush;
}

void ObstacleTracking::calcFeatureMatches(unsigned int relPastFrame)
{
  try {
    FrameSPtr f = frameBuffer.get(relPastFrame);
    FrameSPtr lf = frameBuffer.get(relPastFrame+1);
    cerr << endl << "calculating feature matches..." << flush;
    calculateMatches(lf, f);
  } catch (std::range_error e) {
    // do nothing
  }
  cerr << "done" << flush;
}

void ObstacleTracking::calcRegistration(unsigned int relPastFrame)
{
  try {
    cerr << endl << "calculating registration..." << flush;
    registrationInit(relPastFrame);
    while (registrationNext()) {};
    cerr << "registered " << registrator->nbRegistered << " tracks" << flush;
  } catch (std::range_error e) {
    // do nothing
  }
}

void ObstacleTracking::calcSegmentation(unsigned int relPastFrame)
{
  try {
    FrameSPtr f = frameBuffer.get(relPastFrame);
    if (registrator.get() != NULL) {
      registrator->flushTrackProjListAndExtend();
    }
    segment(f);
    if (registrator.get() != NULL) {
      decideMerge(f, registrator, false); // really decide about merging
    }
    FrameSPtr lf = frameBuffer.get(relPastFrame+1);
    //linkSegments(f, lf); // builds a temporary TwoLayerGraph and saves it to file "segLinks.graphml"
  } catch (std::range_error e) {
    // do nothing
  }
}

void ObstacleTracking::calcTracks(unsigned int relPastFrame)
{
  try {
    FrameSPtr f = frameBuffer.get(relPastFrame);
    cerr << endl << "creating tracks..." << flush;
    f->checkConsistency();
    mergeSplitCreateTracks(f, registrator);
    if (registrator.get() != NULL) {
      BOOST_FOREACH(TrackRegistration::SubTrackSPtr subtrack, registrator->subsampledTracks)
        subtrack->trackBufferIdx++;
    }
    cerr << "done" << flush;
    f->checkConsistency();
    outputLTTTrackDetails(f, false);
  } catch (std::range_error e) {
    // do nothing
  }
}

void ObstacleTracking::registrationInit(unsigned int relPastFrame)
{
  cout << "      RegInit      " << flush;
  try {
    FrameSPtr cf = frameBuffer.get(relPastFrame);
    FrameSPtr lf = frameBuffer.get(relPastFrame+1);
//    visibilityNextFrame(lf, cf);
//    visibilityPreviousFrame(cf, lf);
//    calculateMatches(lf, cf);
    predictTracks(lf, cf); // keep here, so registration can be reinitiated
    registrator = TrackRegistration::SPtr(new TrackRegistration(cf, lf, imageProjector));
  } catch (std::range_error e) {
    // do nothing
  }
  cout << "      RegInitFinit      " << flush;
}

bool ObstacleTracking::registrationNext(bool retFalseOnLayerChange)
{
  cout << "      Next      " << flush;
  bool r = registrator->nextTrack();
  if (retFalseOnLayerChange && r)
    r = ((**registrator->currTrackRegistered).trackBufferIdx == (**registrator->nextTrackToRegister).trackBufferIdx);
  return r;
}


