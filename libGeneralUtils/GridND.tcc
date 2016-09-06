/*!
    \file   GridND.tcc
    \brief  Template Implementation for the GridND class
    \author	Frank Moosmann (<moosmann@mrt.uka.de>),
    \date		2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#include <iostream>
#include <list>
#include <stdexcept>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO
#include <boost/foreach.hpp>

#include "MatrixDefs.hpp"

template <size_t const N, typename T, typename C>
GridND<N,T,C>::GridND(const double cellExtent_) : cellExtent(cellExtent_), totalPtCount(0) {
	assert(cellExtent > 0);
}

template <size_t const N, typename T, typename C>
GridND<N,T,C>::~GridND() {
}


template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::const_iterator
GridND<N,T,C>::begin(const T &point) const {
	if (cellGrid.empty()) return end();
  return begin(getIndex(point));
};

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::iterator
GridND<N,T,C>::begin(const T &point) {
	if (cellGrid.empty()) return end();
  return begin(getIndex(point));
};

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::const_iterator
GridND<N,T,C>::end(const T &point) const {
	if (cellGrid.empty()) return end();
  return end(getIndex(point));
};

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::iterator
GridND<N,T,C>::end(const T &point) {
	if (cellGrid.empty()) return end();
	return end(getIndex(point));
};


template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::const_iterator
GridND<N,T,C>::begin(const index_t &idx) const {
  if (cellGrid.empty()) return end();
  GridConstIterator idxIt = cellGrid.find(idx);
  if (idxIt == cellGrid.end()) return end();
  GridConstIterator idxItN = idxIt; idxItN++;
  return const_iterator(idxIt, idxItN, idxIt->second.begin()); // ok, as lists must never be empty
};

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::iterator
GridND<N,T,C>::begin(const index_t &idx) {
  if (cellGrid.empty())
    return end();
  GridIterator idxIt = cellGrid.find(idx);
  if (idxIt == cellGrid.end())
    return end();
  GridIterator idxItN = idxIt; idxItN++;
  return iterator(idxIt, idxItN, idxIt->second.begin()); // ok, as lists must never be empty
};

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::const_iterator
GridND<N,T,C>::end(const index_t &idx) const {
  if (cellGrid.empty()) return end();
  GridConstIterator idxIt = cellGrid.find(idx);
  if (idxIt == cellGrid.end()) return end();
  ++idxIt;
  if (idxIt == cellGrid.end()) return end();
  return const_iterator(idxIt, idxIt, CellConstIteratorT(0));
};

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::iterator
GridND<N,T,C>::end(const index_t &idx) {
  if (cellGrid.empty()) return end();
  GridIterator idxIt = cellGrid.find(idx);
  if (idxIt == cellGrid.end())
    return end();
  ++idxIt;
  if (idxIt == cellGrid.end())
    return end();
  return iterator(idxIt, idxIt, CellIteratorT(0));
};


template <size_t const N, typename T, typename C>
void
GridND<N,T,C>::getIndex(const T &point, index_t &idx) const {
	idx.resize(N);
	int ii;
	for (unsigned int i=0; i<N; ++i) {
		ii = (int)(point(i)/cellExtent);
		if (point(i) < 0) ii-=1;	// correct neg. indices due to rounding
		idx[i] = ii;
	}
}

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::index_t
GridND<N,T,C>::getIndex(const T &point) const {
	index_t i;
	getIndex(point,i);
	return i;
}

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::coordinate_t
GridND<N,T,C>::getCenter(const index_t &idx) const {
	coordinate_t c(N);
	for (unsigned int i=0; i<N; ++i) {
		double d = ((double)(idx[i]) + 0.5)*cellExtent;
		c[i] = d;
	}
	return c;
}

template <size_t const N, typename T, typename C>
void
GridND<N,T,C>::insert(const T &point) {
	index_t idx =	getIndex(point);
	GridIterator it = cellGrid.find(idx);
	if (it == cellGrid.end()) {
		BOOST_AUTO(itinspair,cellGrid.insert(GridElementType(idx,CellT())));
		it = itinspair.first;
	}
  it->second.push_back(point); // inserts a copy of point
  ++totalPtCount;
}

template <size_t const N, typename T, typename C>
size_t
GridND<N,T,C>::remove(const T &point) {
	index_t idx =	getIndex(point);
	GridIterator git = cellGrid.find(idx);
	size_t nbElmRemoved = 0;
	if (git != cellGrid.end()) {
		nbElmRemoved = git->second.size();
		git->second.remove(point);
		size_t nbpNew = git->second.size();
		nbElmRemoved -= nbpNew;
		totalPtCount -= nbElmRemoved;
		if (nbpNew == 0) {// cell became empty -> delete
			cellGrid.erase(git);
		}
	}
	return nbElmRemoved;
}

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::iterator
GridND<N,T,C>::erase(iterator it) {
	iterator next = it; ++next;
	if (it.ceIt == cellGrid.end())
		throw std::out_of_range("GridND::erase(iterator): iterator is invalid");
	it.ceIt->second.erase(it.ptIt); // erase element from cell
	totalPtCount--;
	if (it.ceIt->second.size() == 0) { // cell became empty -> delete
		cellGrid.erase(it.ceIt);
	}
	return next;
}

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::grid_iterator
GridND<N,T,C>::erase(GridND<N,T,C>::grid_iterator it) {
	grid_iterator next = it; ++next;
	if (it == cellGrid.end())
		throw std::out_of_range("GridND::erase(iterator): iterator is invalid");
	totalPtCount -= it->second.size();
	cellGrid.erase(it);
	return next;
}

template <size_t const N, typename T, typename C>
void GridND<N,T,C>::clear() {
	cellGrid.clear();
	totalPtCount = 0;
}

template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::iterator
GridND<N,T,C>::find(const T &point) {
  // first, find corresponding cell
  index_t idx = getIndex(point);
  GridIterator cellIt = cellGrid.find(idx);
  GridIterator cellEnd = cellGrid.end();
  if (cellIt == cellEnd)
    return end();
  CellT &cl = cellIt->second;
  // if cell exists, find element
  CellIteratorT pIt = std::find(cl.begin(), cl.end(), point); // TODO (9): for some containers, an optimized find() could be used
  if (pIt == cl.end())
    return end();
  return iterator(cellIt, cellEnd, pIt);
}

//template <size_t const N, typename T, typename C>
//typename GridND<N,T,C>::iterator
//GridND<N,T,C>::find_lazy_neighbor(const T &point, double *dist) {
//	const_iterator nIt = ((const GridND<N,T,C>*)this)->find_lazy_neighbor(point, dist);
//	return iterator(nIt);
////	return iterator(GridIterator(nIt.ceIt), GridIterator(nIt.ceEnd), CellListIterator(nIt.ptIt));
//}
// TODO(9): copied code to get non-const function. How can you do this nicely?
template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::iterator
GridND<N,T,C>::find_lazy_neighbor(const T &point, double *dist) {
	index_t idx =	getIndex(point);
	GridIterator cellIt = cellGrid.find(idx);
	GridIterator cellEnd = cellGrid.end();
	if (cellIt == cellEnd)
		return end();
	CellT &cl = cellIt->second;
	// search over all elements in the cell to find closest neighbor:
	CellIteratorT pIt = cl.begin();
	CellIteratorT pItNN = pIt;
	double minDSq = DBL_MAX;
	for (; pIt!=cl.end(); ++pIt) {
		const T &n = *pIt; // retrieve point
		double dsq = 0; // calculate its distance
		for (unsigned int i=0; i<N; ++i) {
			dsq += pow(point(i)-n(i),2);
		}
		if (dsq < minDSq) { // accept if so far closest point
			minDSq = dsq;
			pItNN = pIt;
		}
	}
	assert(minDSq != DBL_MAX);
	if (dist != NULL)
		*dist = sqrt(minDSq);
	return iterator(cellIt, cellEnd, pItNN);
}
template <size_t const N, typename T, typename C>
typename GridND<N,T,C>::const_iterator
GridND<N,T,C>::find_lazy_neighbor(const T &point, double *dist) const {
	index_t idx =	getIndex(point);
	GridConstIterator cellIt = cellGrid.find(idx);
	GridConstIterator cellEnd = cellGrid.end();
	if (cellIt == cellEnd)
		return end();
	const CellT &cl = cellIt->second;
	// search over all elements in the cell to find closest neighbor:
	CellConstIteratorT pIt = cl.begin();
	CellConstIteratorT pItNN = pIt;
	double minDSq = DBL_MAX;
	for (; pIt!=cl.end(); ++pIt) {
		const T &n = *pIt; // retrieve point
		double dsq = 0; // calculate its distance
		for (unsigned int i=0; i<N; ++i) {
			dsq += pow(point(i)-n(i),2);
		}
		if (dsq < minDSq) { // accept if so far closest point
			minDSq = dsq;
			pItNN = pIt;
		}
	}
	assert(minDSq != DBL_MAX);
	if (dist != NULL)
		*dist = sqrt(minDSq);
	return const_iterator(cellIt, cellEnd, pItNN);
}

// TODO(9): copied code to get const function. How can you do this nicely?
template <size_t const N, typename T, typename C>
bool GridND<N,T,C>::find_lazy_neighbor(const T &point, T &neighbor) const {
	const_iterator nIt = find_lazy_neighbor(point, NULL);
	if (nIt == end())
		return false;
	neighbor = *nIt;
	return true;
}

template <size_t const N, typename T, typename C>
bool GridND<N,T,C>::find_lazy_neighbor(const T &point, T &neighbor, double &dist) const {
	const_iterator nIt = find_lazy_neighbor(point, &dist);
	if (nIt == end())
		return false;
	neighbor = *nIt;
	return true;
}

template <size_t const N, typename T, typename C>
template <typename InsertIterator>
void GridND<N,T,C>::cubicNeighborRecursiveHelper(InsertIterator ii, index_t idx, unsigned int dim, index_t center_idx, int radius) {
	//std::cout << "entering at idx " << idx[0] << "/" << idx[1] << "/" << idx[2] << " permut dim " << dim << " from " << center_idx[dim]-radius << " to " << center_idx[dim]+radius << std::endl;
	for (int val = center_idx[dim]-radius; val <= center_idx[dim]+radius; ++val) {
		idx[dim] = val;
		if (dim+1 < N) {
			cubicNeighborRecursiveHelper(ii, idx, dim+1, center_idx, radius);
		} else {
			// if cell exists -> add
			GridIterator cellIt = cellGrid.find(idx);
			if (cellIt != cellGrid.end()) {
				GridIterator cellItN = cellIt; cellItN++;
				CellIteratorT elem = cellIt->second.begin();
				CellIteratorT elemEND = cellIt->second.end();
				while (elem != elemEND) {
					//++ii = iterator(cellIt, cellItN, elem);
					*ii++ = *elem;
					++elem;
				}
			}
		}
	}
}

template <size_t const N, typename T, typename C>
template <typename InsertIterator>
void GridND<N,T,C>::cubicCellNeighborRecursiveHelper(InsertIterator ii, index_t idx, unsigned int dim, index_t center_idx, int radius) {
	if (N == 3) {
//		int nbFound = 0;
		for (int val1 = center_idx[1]-radius; val1 <= center_idx[1]+radius; ++val1) {
			for (int val2 = center_idx[2]-radius; val2 <= center_idx[2]+radius; ++val2) {
				for (int val3 = center_idx[3]-radius; val3 <= center_idx[3]+radius; ++val3) {
					idx[0] = val1; idx[1] = val2; idx[2] = val3;
					GridIterator cellIt = cellGrid.find(idx);
					if (cellIt != cellGrid.end())
//						++nbFound;
						*ii++ = cellIt;
				}
			}
		}
//		if (nbFound >= 0) return;
		return;
	}
	std::cout << "oops" << std::endl;
	for (int val = center_idx[dim]-radius; val <= center_idx[dim]+radius; ++val) {
		idx[dim] = val;
		if (dim+1 < N) {
			cubicCellNeighborRecursiveHelper(ii, idx, dim+1, center_idx, radius);
		} else {
			// if cell exists -> add
			GridIterator cellIt = cellGrid.find(idx);
			if (cellIt != cellGrid.end())
				*ii++ = cellIt;
		}
	}
}

template <size_t const N, typename T, typename C>
template <typename InsertIterator>
unsigned int GridND<N,T,C>::SampleQueue::sample(const unsigned int desiredCount, InsertIterator container)
{
  assert(desiredCount <= elements.size());
  if (desiredCount == elements.size()) {
    for (size_t i=0; i<elements.size(); ++i) {
      *container++ = elements[i];
    }
    return elements.size();
  }
  // TODO (8): use random sampling that guarantees that desiredCount is met
  double sampFac = (double)desiredCount / cumWeight;
  double accSampling = 1.0f; // always select the first element
  unsigned int remaining = desiredCount;
  for (size_t i=0; i<elements.size(); ++i) {
    accSampling += weights[i];
    if (accSampling >= sampFac) {
      accSampling -= sampFac;
      *container++ = elements[i];
      remaining--;
    }
    if (remaining == 0) break;
  }
  return desiredCount - remaining;
}

template <size_t const N, typename T, typename C>
template <typename InsertIterator>
unsigned int GridND<N,T,C>::uniform_sample(const unsigned int desiredCount, InsertIterator container, weight_func weight) const
{
	using namespace std;
	bool useWeighting = !weight.empty();

	unsigned int remaining = desiredCount;
	if (desiredCount > getCellCount()) {
		// we have to pick at least one out of each cell
		// so the challenge is how many we pick out of which cell
		// as there might be cells containing less elements than desiredCount/cellCount
		typedef pair<unsigned int, const CellT*> sortedCellsT;
		multimap<unsigned int, const CellT*> sizeSortedCells;
		for (GridConstIterator it = cellGrid.begin(); it != cellGrid.end(); ++it) {
			const CellT *c = &(it->second);
			sizeSortedCells.insert(sortedCellsT(c->size(),c));
		}
		// now sample into container, starting with the sparsest cell
		unsigned int currentSampleCount;
	  unsigned int remainingCells = getCellCount();
		BOOST_FOREACH(sortedCellsT c, sizeSortedCells) {
			currentSampleCount = ceil((double)remaining/(double)remainingCells); // optimal sample count due to uniform sampling
			currentSampleCount = min(currentSampleCount,c.first); // but: limited by nb samples in cell
			SampleQueue queue;
			for (BOOST_AUTO(it,c.second->begin()); it!=c.second->end(); ++it) {
			  const T &p = *it;
	      double w = useWeighting ? weight(p) : 1.0;
	      queue.push_back(p, w);
			}
			remaining -= queue.sample(currentSampleCount, container);
			if (remaining == 0) break;
			--remainingCells;
		}
	} else {
		// there are more cells than we have to pick
		// thus from each selected cell we pick only the first element
		// so the challenge is to select the cells to pick from
		SampleQueue queue;
		for (GridConstIterator cit = cellGrid.begin(); cit != cellGrid.end(); ++cit) {
      BOOST_AUTO(pit,cit->second.begin()); // simply pick first element in cell
		  if (useWeighting) {
		    // TODO(5): don't pick first element but the one with highest weight
		    // hence, loop through all elements and store keep iterator of the highest weight
		  }
      const T &p = *pit;
      double w = useWeighting ? weight(p) : 1.0;
		  queue.push_back(p, w);
		}
		remaining -= queue.sample(desiredCount, container);
	}
//	cout << endl;
	return desiredCount-remaining;
}

