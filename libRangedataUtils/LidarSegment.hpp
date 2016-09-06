#ifndef LIDARSEGMENT_H_
#define LIDARSEGMENT_H_

#include <boost/shared_ptr.hpp>
#include <list>
#include <map>

#include <MatrixDefs.hpp>
#include "LidarImage.hpp"

namespace mdefs = matrixTools;

/*! 
 * \class LidarSegment
 * \brief A class that can be used to create a "segmentation forest"
 * 
 * Each pixel of the image has a LidarSegment object. The initial state is "valid" and "leaf/noProxy" meaning that each Pixel represents one segment.
 * When two pixels merge into one segment, one pixel will pass all its data to the other "master pixel" and just store its pointer, turning into a "proxy".
 * This can be done multiple times, thus a segmentation forest will be created where a pixel references another pixel, which again references a different one.
 * If a pixel is proxy, all function calls retrieving/modifying information will be directed to the leaf.
 * A pixel can also be deactivated (turned "invalid"), meaning this pixel has no valid segment assigned.  
 */
class LidarSegment
{
public:
//  typedef boost::shared_ptr< LidarSegment > SPtr;
  typedef std::list< ColRowPair >::const_iterator ColRowIterator;
  typedef std::list< ColRowPair >::const_iterator ColRowConstIterator;
  typedef std::map<LidarSegment*, double>::const_iterator NeighborIterator;
  
  LidarSegment();
  virtual ~LidarSegment();
  
  void reset(); //!< reset segement, activating it and removing the use as proxy
  void update(); //!< trace all pointers to other segments down to leafs, assure non-self-neighborness and unique occurance of each neighbor in list
  void merge(LidarSegment *other); //!< merge two segments: moves all contents of other segment into this, turning the other into a proxy

  LidarSegment* getPtr() {return linkTo ? linkTo->getPtr() : this;}; //!< returns the pointer to the leaf in the segmentation forest
  const LidarSegment* getPtr() const {return linkTo ? linkTo->getPtr() : this;}; //!< returns the pointer to the leaf in the segmentation forest
  bool isProxy() const {return (bool)linkTo;}; //!< returns true if the current object serves as proxy, i.e. directs all calls to another segment
  bool isValid() const {return getPtr()->valid;}; //!< returns true if the segment is valid, false otherwise
  void setValid(bool val); //!< can be used to mark the segment as valid/invalid

  double getColorR() const {return getPtr()->colorR;};
  double getColorG() const {return getPtr()->colorG;};
  double getColorB() const {return getPtr()->colorB;};
  void setColor(double r, double g, double b) {LidarSegment* s=getPtr(); s->colorR=r; s->colorG=g; s->colorB=b;};

  void addPix(ColRowPair cr) {LidarSegment* s=getPtr(); s->colRows.push_back(cr); s->pxCount++;};
  unsigned int getPxCount() const {return getPtr()->pxCount;};
  void getCenter(float &colC, float &rowC) const;
  void getColRows(ColRowIterator &begin, ColRowIterator &end, unsigned int &count);
  void getColRows(ColRowConstIterator &begin, ColRowConstIterator &end, unsigned int &count) const;

  static void makeNeighbors(LidarSegment *s1, LidarSegment *s2, double maxSegScore); // connect 2 segments as neighbors with the specified segmentation score
  bool hasNeighbor(LidarSegment *s) {return (neighbors.find(s) != neighbors.end());};
  double getNeighborSegScore(LidarSegment *s); // returns the maximum segmentation score to the specified segment or a value <0 if there is no connection
  void getNeighbors(NeighborIterator &begin, NeighborIterator &end, unsigned int &count);

private:
  LidarSegment* linkTo; // !=NULL if it serves as proxy and redirects all queries to the Segment pointed to
  bool valid;
  double colorR;
  double colorG;
  double colorB;
  unsigned int pxCount; // number of pixels belonging to this segment
  std::list<ColRowPair> colRows; // pixels of this segment. implicitly referencing the frame this segment is stored in
  std::map<LidarSegment*, double> neighbors; // pointers no neighboring segments in the image together with maximum segmentation score of any link (which should be below the segmentation threshold)
};


//struct SegmentSizeComp { bool operator()(const LidarSegment::SPtr &p1, const LidarSegment::SPtr &p2) const 
//                         {return (p1->ptCount > p2->ptCount);}
//                       };
//struct SegmentVisibleSizeComp { bool operator()(const LidarSegment::SPtr &p1, const LidarSegment::SPtr &p2) const 
//                         {return (p1->ptCountVisibleNextFrame > p2->ptCountVisibleNextFrame);}
//                       };

#endif /*LIDARSEGMENT_H_*/
