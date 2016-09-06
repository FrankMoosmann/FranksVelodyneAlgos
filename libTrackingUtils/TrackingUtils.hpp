/*
 * TrackingUtils.h
 *
 *  Created on: Aug 4, 2010
 *      Author: moosmann
 */

#ifndef TRACKINGUTILS_H_
#define TRACKINGUTILS_H_

#include "PngDistanceImage.hpp"

#include "PointCloudTrack.hpp"
#include "LidarFrame.hpp"

void saveToFile(const LidarImage<LidarSegment> &segments, std::string filename);
void loadFromFile(LidarImage<LidarSegment> &segments, std::string filename);

void distanceFromPng(LidarImage<double> &distance, const PngDistanceImage &dstImg);

template <typename T>
class Histogram {
private:
  typedef std::map<T, unsigned int> mapT;
  mapT votes;
public:
  Histogram() {};
  template <class InputIterator>
  Histogram(InputIterator begin, InputIterator end) {
    while (begin != end) vote(*begin++);
  };

  void vote(const T& val) {
    BOOST_AUTO(p, votes.insert( make_pair(val,1) )); // pair<iterator,bool>
    if (!p.second)
      *(p.first)++;
  };

  size_t size() const { return votes.size();};

  T getMaxVotes() const {
    if (votes.empty())
      throw std::out_of_range("");
    typedef typename mapT::const_iterator itT;
    itT i=votes.begin();
    unsigned int maxVote = i->second;
    T t = i->first;
    ++i;
    for (; i!=votes.end(); ++i) {
      if (i->second > maxVote) {
        maxVote = i->second;
        t = i->first;
      }
    }
    return t;
  };
};

// to do: write class for covariance-ellipses (nD)
// currently implemented in Utils/CovarEllipsoidRendering.h
// this class must deliver: angle in 2D case, axis in nD, scales along axis, maybe even whole basis-matrix

// SEGMENT:
//void linkWithTrack(PointCloudTrack *t);
//void removeLowCountLinks(double unlinkThresh); // remove links low count relative to strongest link -->  0<=unlinkThresh<=1
//void replaceTrack(PointCloudTrack *orig, PointCloudTrack *repl);
//PointCloudTrack* getFirstTrack();
//void getTracks(TrackAssocIterator &begin, TrackAssocIterator &end, unsigned int &count);
//void unlinkAllTracks() {tracks.clear();};
//void mergeSegmentsAndCreateTracks(LidarFrame *frame, const double currEgoXVelocity, const mdefs::DSMatrix &initVariance, const mdefs::DSMatrix &predictionVariance);

// TRACK:
// methods concerning temporary frame-based connections
//unsigned int getPixCount() const {return projPixCount;};
//void incPixCount() {projPixCount++;};
//void decPixCount() {if (projPixCount > 0) projPixCount--;};
//void linkWithSegment(LidarSegment* s);
//void unlinkSegment(LidarSegment* s);
//unsigned int getLinkedSegmentCount() {return segments.size();};
//unsigned int getLinkedSegmentLastFrameCount() {return segLastFrame.size();};
//LidarSegment* getFirstSegment() const {return segments.empty() ? NULL : *(segments.begin());};
//LidarSegment* getFirstSegmentLastFrame() const {return segLastFrame.empty() ? NULL : *(segLastFrame.begin());};
//void mergeNonConnectedSegments();
//void splitTrack(PointCloudTrack::SPtr track, LidarFrame *frame, std::vector<PointCloudTrack::SPtr> &tmpList, unsigned int &newIdx);


#endif /* TRACKINGUTILS_H_ */
