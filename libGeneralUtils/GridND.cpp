/*!
    \file   GridND.cpp
    \brief  Implementation of a test routine for the GridND class
    \author  Frank Moosmann (<moosmann@mrt.uka.de>),
    \date    2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#include "GridND.hpp"
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/foreach.hpp>

using namespace std;
using namespace matrixTools;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
///////////// the following lines are a ////////////////
///////////// test implementation for   ////////////////
///////////// GridND template class     ////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
struct GridNDtestDVProxy {
  DVector *v;
  GridNDtestDVProxy (DVector *v_) : v(v_) {};
  inline double operator() (int i) const {return (*v)(i);};
  inline bool operator== (GridNDtestDVProxy const &b) {return (v == b.v);};
  friend std::ostream &operator<< (std::ostream &ostr, const GridNDtestDVProxy &p) {return ostr<<*(p.v);};
};

//struct GridNDtestSurface {
//  DVector position;
//  DVector normal;
//  double variance; // each point as uncertainty relative to current position. assumed that covariance matrix = covariance * I3 (3x3 identity matrix)
//  GridNDtestSurface (DVector &p, DVector &n, double v) : position (p), normal(n), variance(v) {};
//  inline double operator() (int i) const {return position(i);};
//};
//
//struct GridNDtestSurfacePProxy {
//  GridNDtestSurface *s;
//  GridNDtestSurfacePProxy (GridNDtestSurface *s_) : s(s_) {};
//  inline double operator() (int i) const {return (*s)(i);};
//};

bool operator==(const DVector &a, const DVector &b) {
  return (norm_2(a-b) < 0.00001);
};

void GridNDtest(bool verbose) {
  using namespace std;
  bool success = true;

  // setting up 3D grid
  if (verbose) cout << endl << "############# testing 3D grid #############" << endl;
  DVector p1 = ublas::zero_vector<double>(3); p1(0) = 0.01;                 // index 0/0/0
  DVector p2 = ublas::zero_vector<double>(3); p2(1) = 0.02;                 // index 0/0/0
  DVector p3 = ublas::zero_vector<double>(3); p3(1) = 2.0;                   // index 0/19/0 (2x) ??
  DVector p4 = ublas::zero_vector<double>(3); p4(0) = 0.05; p4(2) = 2.0;     // index 0/0/19 ??
  DVector p5 = ublas::zero_vector<double>(3); p5(0) = 2.01;                 // index 20/0/0
  DVector p6 = ublas::zero_vector<double>(3); p6(0) = 2.03;                 // index 20/0/0
  DVector p7 = ublas::zero_vector<double>(3); p7(0) = 2.02;                 // index 20/0/0
  DVector p8 = ublas::zero_vector<double>(3); p8(1) = 2.01;                 // index 0/20/0
  DVector p9 = ublas::zero_vector<double>(3); p9(1) = 2.02;                 // index 0/20/0
  DVector n = ublas::zero_vector<double>(3); n(0) = 1.0;                     // index 9/0/0 ??
  GridND<3,DVector> grid(0.1);
  grid.insert(p1);
  grid.insert(p2);
  grid.insert(p3);
  grid.insert(p3); // double insert
  grid.insert(p4);
  grid.insert(p5);
  grid.insert(p6);
  grid.insert(p7);
  grid.insert(p8);
  grid.insert(p9);
  grid.insert(n);
  DVector pfw = ublas::zero_vector<double>(3); pfw(0) = 5.0; pfw(1) = 10.0; // index 49/9/0 ??
  DVector pS = ublas::zero_vector<double>(3); pS(0) = 0.05;                 // index 0/0/0
  DVector pS2 = ublas::zero_vector<double>(3); pS2(1) = 2.0;                 // index 0/19/0 ??

  if ((success) && (verbose))
      cout << "grid: " << grid << endl;

  // testing 3D grid
  // testing basic information:
  if (success) {
    if (grid.dimensionality() != 3) {
      cerr << endl << "Test failed for member dimensionality()" << flush;
      success = false;
    }
  }
  if (success) {
    if (grid.empty()) {
      cerr << endl << "Test failed for member empty()" << flush;
      success = false;
    }
  }
  if (success) {
    unsigned int nbElem = 0;
    if (verbose) cout << "grid elements:";
    for (GridND<3,DVector>::const_iterator it=grid.begin(); it!=grid.end(); ++it) {
      if (verbose) cout << " " << *it << flush;
      ++nbElem;
    }
    if (verbose) cout << endl;
    if (nbElem != grid.size()) {
      cerr << endl << "Test failed for member size() / iterators" << flush;
      success = false;
    }
  }
  // testing cell access:
  if (success) {
    GridND<3,DVector>::grid_const_iterator itrcend = grid.gridEnd();
    unsigned int nbCells = 0;
    if (verbose)  cout << "cells:" << endl;
    for (GridND<3,DVector>::grid_const_iterator itrc = grid.gridBegin(); itrc!=itrcend; ++itrc) {
      if (verbose) cout << itrc->first[0] << "/" << itrc->first[1] << "/" << itrc->first[2] << " -> " << itrc->second.size() << " points" << endl;
      ++nbCells;
    }
    if (nbCells != grid.getCellCount()) {
      cerr << endl << "Test failed for getCellCount() / gridBegin()/gridEnd()" << flush;
      success = false;
    }
    if (grid.getCellCount() != 6) {
      cerr << endl << "Test failed for getCellCount(): nbCells=" << grid.getCellCount()<<flush;
      success = false;
    }
  }
  if (success) {
    list< DVector > elems;
    list< GridND<3,DVector>::grid_iterator > cells;
    back_insert_iterator< list<DVector> > elemsInserter(elems);
    back_insert_iterator< list< GridND<3,DVector>::grid_iterator > > cellsInserter(cells);
    grid.get_cubic_neighbors(elemsInserter, pS2, 0);
    grid.getCubicCellNeighbors(cellsInserter, pS2, 0);
    if (verbose) cout << "when searching cubic neighborhood of radius 0:" << endl;
    if (verbose) cout << "found " << elems.size() << "elements" << endl;
    if (verbose) cout << "found " << cells.size() << "cells" << endl;
    if (elems.size() != 2) {
      cerr << endl << "Test 1 failed for get_cubic_neighbors()" << flush;
      success = false;
    }
    if (cells.size() != 1) {
      cerr << endl << "Test 1 failed for getCubicCellNeighbors()" << flush;
      success = false;
    }
  }
  if (success) {
    list< DVector > elems;
    list< GridND<3,DVector>::grid_iterator > cells;
    back_insert_iterator< list<DVector> > elemsInserter(elems);
    back_insert_iterator< list< GridND<3,DVector>::grid_iterator > > cellsInserter(cells);
    grid.get_cubic_neighbors(elemsInserter, pS2, 1);
    grid.getCubicCellNeighbors(cellsInserter, pS2, 1);
    if (verbose) cout << "when searching cubic neighborhood of radius 1:" << endl;
    if (verbose) cout << "found " << elems.size() << "elements" << endl;
    if (verbose) cout << "found " << cells.size() << "cells" << endl;
    if (elems.size() != 4) {
      cerr << endl << "Test 2 failed for get_cubic_neighbors()" << flush;
      success = false;
    }
    if (cells.size() != 2) {
      cerr << endl << "Test 2 failed for getCubicCellNeighbors()" << flush;
      success = false;
    }
  }
  // testing neighbor access:
  if (success) {
    if (verbose) cout << "grid elements around point " << n << ": ";
    GridND<3,DVector>::const_iterator itrpend = grid.end(n);
    unsigned int nbElem = 0;
    for (GridND<3,DVector>::const_iterator itrp = grid.begin(n); itrp!=itrpend; ++itrp) {
      if (verbose) cout << *itrp << " " << flush;
      ++nbElem;
    }
    if (verbose) cout << endl;
    if (nbElem != 1) {
      cerr << endl << "Test 1 failed for 1-cell-only iterators" << flush;
      success = false;
    }
  }
  if (success) {
    if (verbose) cout << "grid elements around point " << p1 << ": ";
    GridND<3,DVector>::const_iterator itrpend = grid.end(p1);
    unsigned int nbElem = 0;
    for (GridND<3,DVector>::const_iterator itrp = grid.begin(p1); itrp!=itrpend; ++itrp) {
      if (verbose) cout << *itrp << " " << flush;
      ++nbElem;
    }
    if (verbose) cout << endl;
    if (nbElem != 2) {
      cerr << endl << "Test 2 failed for 1-cell-only iterators (nbEleme = " << nbElem << ")" << flush;
      success = false;
    }
  }
  if (success) {
    DVector pSn;
    if (grid.find_lazy_neighbor(pS,pSn))
      if (verbose) cout << "neighbor of " << pS << " is " << pSn << endl;
    if (!(pSn == p1)) {
      cerr << endl << "Test failed for lazy neighbor search" << flush;
      success = false;
    }
  }
  // testing sampling:
  if (success) {
    if (verbose) {

      list<DVector> samples;
      insert_iterator< list<DVector> > iit(samples, samples.begin());
      if (verbose) cout << grid.uniform_sample(9,iit);
      if (verbose) cout << " samples:";
      BOOST_FOREACH(DVector v, samples) {
        if (verbose) cout << " " << v;
      }
      if (verbose) cout << endl;
    }
  }
  // testing deletion:
// The following does not work because DVector has no operator== defined
//  if (success) {
//    unsigned int sizeBefore = grid.size();
//    unsigned int removed = grid.remove(p3);
//    unsigned int sizeAfter = grid.size();
//    if ((removed != 2) || (sizeBefore != removed+sizeAfter)) {
//      cerr << endl << "Test failed for deleting a double occuring point" << flush;
//      success = false;
//    }
//  }
//  if (success) {
//    unsigned int sizeBefore = grid.size();
//    unsigned int removed = grid.remove(p3); // removing the same point once more
//    unsigned int sizeAfter = grid.size();
//    if ((removed != 0) || (sizeBefore != sizeAfter)) {
//      cerr << endl << "Test failed for deleting a non-existing point" << flush;
//      success = false;
//    }
//  }
  if (success) {
    GridND<3,DVector>::iterator nIt = grid.find_lazy_neighbor(pS);
    if (nIt == grid.end()) {
      cerr << endl << "Test 1  failed for lazy iterator-neighbor search" << flush;
      success = false;
    } else {
      if (verbose) cout << "neighbor of " << pS << " is " << *nIt << endl;
    }
    if (norm_2(*nIt-p1) > 0.00001) {
      cerr << endl << "Test 2 failed for lazy iterator-neighbor search" << flush;
      success = false;
    }
    unsigned int sizeBefore = grid.size();
    if (verbose) cout << "after removing iterator on " << *nIt << flush;
    nIt = grid.erase(nIt);
    if (verbose) cout << " iterator points to " << *nIt << endl;
    unsigned int sizeAfter = grid.size();
    if (sizeBefore-sizeAfter != 1) {
      cerr << endl << "Test 1 failed for erasing an iterator" << flush;
      success = false;
    }
    if (norm_2(*nIt-p1) < 0.00001) { // nIt should now point to a different element
      cerr << endl << "Test 2 failed for erasing an iterator" << flush;
      success = false;
    }
  }



  if (verbose) cout << endl << "############# testing 2D grid #############" << endl;
  n = ublas::zero_vector<double>(3); n(0) = 1.0;
  GridNDtestDVProxy p1p(&p1);
  GridNDtestDVProxy p2p(&p2);
  GridNDtestDVProxy p3p(&p3);
  GridNDtestDVProxy p4p(&p4);
  GridNDtestDVProxy p5p(&p5);
  GridNDtestDVProxy p6p(&p6);
  GridNDtestDVProxy p7p(&p7);
  GridNDtestDVProxy np(&n);
  GridND<2,GridNDtestDVProxy> gridp;
  gridp.insert(p1p);
  gridp.insert(p2p);
  gridp.insert(p3p);
  gridp.insert(p3p); // double insert
  gridp.insert(p4p);
  gridp.insert(p5p);
  gridp.insert(p6p);
  gridp.insert(p7p);
  gridp.insert(np);

  if (verbose) {
    cout << "grid: " << gridp << endl;
    GridNDtestDVProxy pSp(&pS);
    cout << "n:" << n << endl;
    if (gridp.find_lazy_neighbor(pSp,np))
      cout << "neighbor of " << pSp << " is " << np << endl;
    cout << "np:" << np << endl;
    cout << "n:" << n << endl;
    cout << gridp << endl;
    gridp.remove(p3p);
    cout << gridp << endl;

    /*
    GridNDtestSurface s1(p1, n, 1.0);
    GridND<3,GridNDtestSurface> sgrid;
    sgrid.insert(s1);
    if (verbose) cout << sgrid << endl;
    */
  }
};
