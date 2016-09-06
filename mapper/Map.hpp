#ifndef MAP_H_
#define MAP_H_
/*
 * Map.h
 *
 *  Created on: Nov 8, 2010
 *      Author: moosmann
 */
#include <deque>
#include <ANN/ANN.h> // (approximate) nearest neighbor
#include <GridND.hpp>
#include <MatrixDefs.hpp>
#include <ICP.hpp>
#include "Frame.hpp"
#include "Surface.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/function.hpp>


// TODO (3): restructure class (e.g. shift registration to mapper)
class Map {
public:
  typedef GridND< 3, SurfaceProxy, SingleElementCell<SurfaceProxy> >  GridT;
  typedef std::deque<ScanPose>  TrajT;
  typedef GridT::iterator iterator;
  typedef GridT::const_iterator const_iterator;
  typedef TrajT::iterator trajectory_iterator;
  typedef TrajT::const_iterator trajectory_const_iterator;
  // TODO (8): move registration into separate class
  typedef std::list<matrixTools::DVector>::const_iterator const_iterator_samplepoints;
  enum ExportType {Surfaces, Points, PointIntens};
  typedef boost::function<void (double x, double y, double z, unsigned char intensity, unsigned int hitCnt)> CellDeleteNotifier; //!< returns corresponding point and distance (in case returned dist == DBL_MAX correspondence was rejected)

public: // TODO (6): change back to private after debugging finished
//private:
  class SearchResults {
  private:
    const Map &map;
    unsigned int currIndex;
  public:
    SearchResults(const Map &map_) : map(map_), currIndex(0) {};
    ~SearchResults() {};

    bool next() {++currIndex; return (currIndex < map.nbValidResults);}; //!< increments iterator and returns false in case end is reached
    const matrixTools::DVector& getCorrespPoint() {if (currIndex >= map.nbValidResults) throw std::out_of_range("Map::SearchResults: invalid index"); /*std::cout << map.neighbors[currIndex]->position << std::endl;*/ return map.neighbors[currIndex]->position;}; //!< \throws std::exception if invalid correspondence
    const matrixTools::DVector& getCorrespOrigPoint() {if (currIndex >= map.nbValidResults) throw std::out_of_range("Map::SearchResults: invalid index"); /*std::cout << map.neighbors[currIndex]->position << std::endl;*/ return map.neighbors[currIndex]->origPosition;}; //!< \throws std::exception if invalid correspondence
    const matrixTools::DVector& getCorrespNormal() {if (currIndex >= map.nbValidResults) throw std::out_of_range("Map::SearchResults: invalid index"); /*std::cout << map.neighbors[currIndex]->normal << std::endl;*/ return map.neighbors[currIndex]->normal;}; //!< \throws std::exception if invalid correspondence
    double getNormConf() {if (currIndex >= map.nbValidResults) throw std::out_of_range("Map::SearchResults: invalid index"); return map.neighbors[currIndex]->normalConfidence;}; //!< \throws std::exception if invalid correspondence
    //double getNormStdDevRAD() {if (currIndex >= map.nbValidResults) throw std::out_of_range("Map::SearchResults: invalid index"); return map.neighbors[currIndex]->normalStdDevRAD;}; //!< \throws std::exception if invalid correspondence
    double getDist() {if (currIndex >= map.nbValidResults) throw std::out_of_range("Map::SearchResults: invalid index"); return map.dists[currIndex];}; //!< \throws std::exception if invalid correspondence
    double getWeight() {if (currIndex >= map.nbValidResults) throw std::out_of_range("Map::SearchResults: invalid index"); return map.weights[currIndex];}; //!< \throws std::exception if invalid correspondence
  };
  friend class SearchResults; // allow SearchResults access to private variables of ProjectedNeighborSearch

  const double        mapResol; //!< resolution of map in meter
  const double        cellDiag; //!< extent of diagonal of a cell
  const double        maxDist; //!< maximum scanning distance to use points from
  const double         weightZ; //!< weight z-component in kd-tree search space
  const unsigned int  treeDim; //!< dimensionality of kd-tree search space: x,y,z,nx,ny,nz
  const unsigned int  icpNbCorresp; //!< number of correspondences to use per search point within ICP
  const double        icpRegCoeff; //!< regularization coefficient for ICP
  const unsigned int   nbBlocks2Unwarp; //! number of chunks of an image to unwarp
  const unsigned int  filterNCellRad;
  const double         filterMinWSum;
  const double         filterMinNormalConfidence;
  const unsigned int   filterNbNSearched;
  const unsigned int   adeptNbNSearched;
  const unsigned int   nbNMinFound;
  const double         nMaxDist;
  const double         minNormalCoincidence;

  // grid holding the complete map:
  TrajT                trajectory; //!< series of poses of the scans incorporated into the map
  TrajT                insTrajectory; //!< series of poses of the scans incorporated into the map
  GridT                *pointCloudGrid; //!< 3D point index on the point cloud, uses (shared) pointers to the real surface points
  unsigned int        nbAddedScans;
  CellDeleteNotifier  cellDeleteNotifier;

  // kd-Tree storage for neighbor-search:
  unsigned int        lastNbAddedScans; //!< used to decide on re-built of tree
  ANNpointArray       dataPts;
  Surface::SPtr       *idxMap;
  ANNkd_tree           *kdTree;
  // query variables, might even be changed in "const" retrieval-functions:
  mutable ANNpoint      queryPt; // query point, transformed in kd-tree format
  mutable unsigned int  resultArraySize; // holds size of array-index
  mutable unsigned int  nbValidResults; // number of valid results in arrays (<=arraySize)
  mutable ANNidxArray   nnIdx; // query-result: index
  mutable ANNdistArray  dists; // query-result: distance to neighbor (after kd-query: 6D-distance, thus has to be overridden)
  mutable Surface::SPtr *neighbors; // query-result
  mutable double *weights; // query-result

  // sub-sampled frame for registration (buffered for visualization)
  std::list<matrixTools::DVector> pts,nrs;
  std::list<double> nCs;//,nSs;
  matrixTools::DMatrix subPtsHTM;
  ICP::ICPBase *icp;

  //! method to cleanly delete search index
  void deinitKdTree();
  //! resets search and thus accounts for new data in map. Search is limited to points within maxDist around a given center
  void reinitializeSearch(const matrixTools::DVector &center, double maxDistance);
  //! \returns SearchResults object containing nnCount neighbors for given search point + normal
  SearchResults findClosestNeighbors(const matrixTools::DVector &searchPoint, const matrixTools::DVector &searchNormal, /*double normalConfidence,*/ unsigned int nnCount) const;
  //! \returns true if the surface was accepted (even returns true if an existing surface was replaced by s)
  bool addToMap(Surface::SPtr s, bool adept);
  //! modifies the position. this is correctly handled by contingently moving the point into another grid cell
  iterator modifyPos(iterator si, matrixTools::DVector pos);

  //! for storing and loading this class with boost:serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    (void)version;
    unsigned int f = 0; // dummy
    int64_t t = 0;
    matrixTools::DMatrix m = matrixTools::DIdMatrix(4);
    unsigned int n;
    ar & boost::serialization::make_nvp("mapResol",const_cast<double &>(mapResol));
    ar & boost::serialization::make_nvp("maxDist",const_cast<double &>(maxDist));
    ar & boost::serialization::make_nvp("weightZ",const_cast<double &>(weightZ));
    ar & boost::serialization::make_nvp("treeDim",const_cast<unsigned int &>(treeDim));
    ar & boost::serialization::make_nvp("icpNbCorresp",const_cast<unsigned int &>(icpNbCorresp));
    ar & boost::serialization::make_nvp("icpRegCoeff",const_cast<double &>(icpRegCoeff));
    ar & boost::serialization::make_nvp("filterEachXFrames",f);
    ar & boost::serialization::make_nvp("nbBlocks2Unwarp",const_cast<unsigned int &>(nbBlocks2Unwarp));
//    ar & BOOST_SERIALIZATION_NVP(pointCloudGrid);
    ar & boost::serialization::make_nvp("nbAddedScans",n);
    //ar & BOOST_SERIALIZATION_NVP(timeDiffBuffer);
    ar & boost::serialization::make_nvp("avgFrameDuration",t);
    ar & boost::serialization::make_nvp("currFrameMoveHTM",m);
    ar & boost::serialization::make_nvp("currFrameDiffHTM",m);
    ar & boost::serialization::make_nvp("currFrameDiffTime",t);
    //lastNbAddedScans ???
  }
  // TODO: split serialization members and initialize temporary members on load

public:
  Map(double mapResol, double maxDist);
  virtual ~Map();

  void registerNotifier(CellDeleteNotifier notifier) {cellDeleteNotifier = notifier;};

  // TODO (8): remove, introduced for debug purposes only (split algorithm into parts)
  iterator csi;
  double avgC;
  unsigned int nbChanged;
  unsigned int nbDeleted;

  unsigned int nbAddingChangedIdx;
  unsigned int nbAddingDidntChangeIdx;

  iterator begin() {return pointCloudGrid->begin();};
  iterator end() {return pointCloudGrid->end();};
  const_iterator begin() const {return pointCloudGrid->begin();};
  const_iterator end() const {return pointCloudGrid->end();};
  trajectory_iterator beginTrajectory() {return trajectory.begin();};
  trajectory_iterator endTrajectory() {return trajectory.end();};
  trajectory_const_iterator beginTrajectory() const {return trajectory.begin();};
  trajectory_const_iterator endTrajectory() const {return trajectory.end();};
  trajectory_iterator beginInsTrajectory() {return insTrajectory.begin();};
  trajectory_iterator endInsTrajectory() {return insTrajectory.end();};
  trajectory_const_iterator beginInsTrajectory() const {return insTrajectory.begin();};
  trajectory_const_iterator endInsTrajectory() const {return insTrajectory.end();};
  const_iterator_samplepoints beginSubsampledPoints() const {return pts.begin();};
  const_iterator_samplepoints endSubsampledPoints() const {return pts.end();};
  matrixTools::DMatrix getSubPtsHTM() const {return subPtsHTM;};

  void registerScan(LFrameSPtr frame, LFrameSPtr lastFrame, unsigned int sampleCnt = UINT_MAX, bool unwarp = false);
  void reg1(LFrameSPtr frame, LFrameSPtr lastFrame, unsigned int sampleCnt = UINT_MAX, bool unwarp = false);
  void reg2();
  void reg3(LFrameSPtr frame);
  void addToMap(LFrameSPtr frame, bool adept, bool unwarp, bool remOldMap, bool remFarPts);
  void filter();
  void filter_init();
  void filter_next();
  void filter_finish();

  void saveToFile(std::string file) const;
  void loadFromFile(std::string file);
  void exportMap(std::string filename, unsigned int resolutionFactor = 1, ExportType type = Surfaces);
  void importMap(std::string filename);
  void importPointCloud(std::string filename, bool withIntensity);
  void exportTrajectories(std::string basefilename);
  void importTrajectories(std::string basefilename);
};

BOOST_CLASS_VERSION(Map, 1)

#endif /* MAP_H_ */
