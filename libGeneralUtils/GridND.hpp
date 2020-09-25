/*!
    \file   GridND.h
    \brief  Provides a template container for storing data in spatial grid-cells
    \author  Frank Moosmann (<moosmann@mrt.uka.de>),
    \date    2009

    Copyright: Karlsruhe Institute of Technology (KIT)
               Institute of Measurement and Control Systems
               All rights reserved
               http://www.mrt.uni-karlsruhe.de
*/
#ifndef GRIDND_H_
#define GRIDND_H_

#include <iostream>
#include <cfloat>
#include <cstddef>
#include <list>
#include <deque>
#include <vector>
#include <map>
#include <unordered_map>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>

/*!
  \class GridND
  \brief Class for efficiently managing objects within a N-dimensional grid-like structure

  This class distributes elements of type T into a ND-grid of cells of type C using hashing.

  The template type T only needs to implement the operator(int) in order
  to access the elements coordinates via the indices (0),(1)..(N-1)
  If elements shall be erased, the operator== must additionally be defined on T.
  If elements shall be searched for, the operator= must additionally be defined on T.

  The cell type C has to implement insertion of elements via push_back(),
  iteration over their elements via begin()/end() delivering both
  iterator and const_iterator and removal of elements via erase(iterator) and remove(T).
  A custom cell implementation is given below (SingleElementCell)

  Only cells that are occupied are stored in memory.
  Insertion does not invalidate any iterator if this is the case for C
*/
template < size_t const N, typename T, typename C = std::list<T> >
class GridND {
////////////////////////////////////////////////////////////////////////////////
////////////////   public typedefs for the start  //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
public:
  // TODO(9): use boost::array for fixed-size containers
  typedef T                                               value_type;
  typedef size_t                                          size_type;
  typedef std::vector<double>                              coordinate_t;
  typedef std::vector<int>                                 index_t;
//  typedef int[N]                                          index_t;
//  typedef boost::array<int>                               index_t;
  typedef boost::function<double (const T &t)>            weight_func;

////////////////////////////////////////////////////////////////////////////////
//////////////// some typedefs and private classes /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
private:
  // Definition of the cell-types:
  // Using std::list has the nice property that insertion does not invalidate iterators
  typedef C                                               CellT;
  typedef typename CellT::iterator                         CellIteratorT;
  typedef typename CellT::const_iterator                   CellConstIteratorT;

  // Definition of the grid-types:
  // Using std::unordered_map has the nice property that insertion does not invalidate iterators
  struct                                                   index_hash
  : std::unary_function<index_t, std::size_t> {
     std::size_t operator() (index_t const& x) const {
       std::size_t seed = 0;
       for (unsigned int i=0; i<N; ++i) {
         boost::hash_combine(seed, x[i]);
       }
       return seed;
     }
  };
  typedef std::unordered_map<index_t, CellT, index_hash>   GridType;
  typedef typename std::pair< index_t, CellT >            GridElementType;
  typedef typename GridType::iterator                     GridIterator;
  typedef typename GridType::const_iterator               GridConstIterator;

  //! \class "iter"  \brief serves as template for "iterator" and "const_iterator"
  template <typename ElemT, typename CellIterator, typename ElemIterator> // all three typenames can switch between const and non-const
  class iter : public boost::iterator_facade<iter<ElemT, CellIterator, ElemIterator>, ElemT, boost::forward_traversal_tag> {
    public:

      iter() {};
      iter(CellIterator ceIt_, CellIterator ceEnd_, ElemIterator ptIt_): ceIt(ceIt_), ceEnd(ceEnd_), ptIt(ptIt_) {};
      template <typename OtherElemT, typename OtherCellIterator, typename OtherElemIterator> // constructor as template to allow construction from const and non-const
      iter(iter<OtherElemT, OtherCellIterator, OtherElemIterator> const& other) : ceIt(other.ceIt), ceEnd(other.ceEnd), ptIt(other.ptIt) {};
//      template <typename OtherElemT, typename OtherCellIterator, typename OtherElemIterator>
//      iter<ElemT, CellIterator, ElemIterator>& operator=(iter<OtherElemT, CellIterator, ElemIterator> const& other) {ceIt = other.ceIt; ceEnd = other.ceEnd; ptIt = other.ptIt; return (*this);};

    private:
      friend class boost::iterator_core_access;
      friend class GridND;
//      template <class> friend class iter;
      CellIterator ceIt;
      CellIterator ceEnd;
      ElemIterator ptIt;

      // template-variant allows comparison between const and non-const, too
      template <typename OtherElemT, typename OtherCellIterator, typename OtherElemIterator>
      bool equal(iter<OtherElemT, OtherCellIterator, OtherElemIterator> const& other) const {
        return (ptIt == other.ptIt);
      };
      ElemT& dereference() const {
        return (*ptIt);
      };
      iter<ElemT, CellIterator, ElemIterator>& increment() {
        ptIt++;
        if (ptIt == ceIt->second.end()) {
          ceIt++;
          if (ceIt == ceEnd)
            ptIt = ElemIterator(0);
          else
            ptIt = ceIt->second.begin();
        }
        return(*this);
      };
  };

  class SampleQueue {
  public:
    SampleQueue() : cumWeight(0) {};
    ~SampleQueue() {};
    void push_back(const T& t, double w) {elements.push_back(t); weights.push_back(w); cumWeight+=w;};
    template <typename InsertIterator>
    unsigned int sample(const unsigned int desiredCount, InsertIterator container);
  private:
    std::deque<T> elements;
    std::deque<double> weights;
    double cumWeight;
  };
////////////////////////////////////////////////////////////////////////////////
////////////////   private methods and variables  //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
private:
  const double     cellExtent;
  GridType         cellGrid;
  unsigned int     totalPtCount;

  template <typename InsertIterator>
  void cubicNeighborRecursiveHelper(InsertIterator ii, index_t idx, unsigned int dim, index_t center_idx, int radius);
  template <typename InsertIterator>
  void cubicCellNeighborRecursiveHelper(InsertIterator ii, index_t idx, unsigned int dim, index_t center_idx, int radius);

  //! for storing and loading this class with boost:serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    (void)version;
    ar & boost::serialization::make_nvp("cellExtent",const_cast<double &>(cellExtent));
    ar & BOOST_SERIALIZATION_NVP(cellGrid);
    ar & BOOST_SERIALIZATION_NVP(totalPtCount);
  };

////////////////////////////////////////////////////////////////////////////////
////////////////   public methods and variables   //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
public:

  typedef iter<T, GridIterator, CellIteratorT>                   iterator; //!< iterates on all elements sorted by cell
  typedef iter<T const, GridConstIterator, CellConstIteratorT>   const_iterator; //!< iterates on all elements sorted by cell
  typedef GridIterator                                          grid_iterator; //!< iterates on the cells
  typedef GridConstIterator                                      grid_const_iterator; //!< iterates on the cells

  GridND(const double cellExtent = 0.1); //!< \param cellExtent specifies the size of the spatial bins in each dimension
  virtual ~GridND();

  void            getIndex(const T &point, index_t &i) const;
  index_t         getIndex(const T &point) const;
  coordinate_t    getCenter(const index_t &i) const;

  size_t           dimensionality() const {return N;}; //!< get the number of dimensions of the grid-structure
  const_iterator   begin() const {return cellGrid.empty() ? end() : const_iterator(cellGrid.begin(), cellGrid.end(), cellGrid.begin()->second.begin());}; //!< returns an iterator on all the points stored within the grid.
  iterator         begin()       {return cellGrid.empty() ? end() :       iterator(cellGrid.begin(), cellGrid.end(), cellGrid.begin()->second.begin());}; //!< returns an iterator on all the points stored within the grid.
  const_iterator   end() const {return const_iterator(cellGrid.end(), cellGrid.end(), CellConstIteratorT(0));}; //!< returns an iterator on all the points stored within the grid.
  iterator         end()       {return       iterator(cellGrid.end(), cellGrid.end(), CellIteratorT(0));}; //!< returns an iterator on all the points stored within the grid.
  const_iterator   begin(const T &point) const;   //!< returns an iterator on all the points within the same grid-cell as the specified point
  iterator         begin(const T &point);         //!< returns an iterator on all the points within the same grid-cell as the specified point
  const_iterator   end(const T &point) const;   //!< returns an iterator on all the points within the same grid-cell as the specified point
  iterator         end(const T &point);         //!< returns an iterator on all the points within the same grid-cell as the specified point
  const_iterator  begin(const index_t &idx) const; //!< returns an iterator on all the points within the same grid-cell as the specified point
  iterator        begin(const index_t &idx);             //!< returns an iterator on all the points within the same grid-cell as the specified point
  const_iterator  end(const index_t &idx) const; //!< returns an iterator on all the points within the grid-cell specified
  iterator        end(const index_t &idx);             //!< returns an iterator on all the points within the grid-cell specified
//  const_iterator   begin(grid_const_iterator git) const; //!< returns an iterator on all the points within the same grid-cell as the specified point
//  iterator         begin(grid_iterator git);             //!< returns an iterator on all the points within the same grid-cell as the specified point
//  const_iterator   end(grid_const_iterator git) const; //!< returns an iterator on all the points within the grid-cell specified
//  iterator         end(grid_iterator git);             //!< returns an iterator on all the points within the grid-cell specified
  grid_const_iterator   gridBegin() const {return cellGrid.begin();}; //!< returns an iterator on the cells of the grid
  grid_iterator         gridBegin()       {return cellGrid.begin();}; //!< returns an iterator on the cells of the grid
  grid_const_iterator   gridEnd() const {return cellGrid.end();}; //!< returns an iterator on the cells of the grid
  grid_iterator         gridEnd()       {return cellGrid.end();}; //!< returns an iterator on the cells of the grid
  grid_const_iterator   fromIterator(const_iterator it) const  {return it.ceIt;}; //!< converts an iterator into a grid iterator
  grid_iterator         fromIterator(iterator it)              {return it.ceIt;}; //!< converts an iterator into a grid iterator
  grid_const_iterator   find(const index_t &idx) const {return cellGrid.find(idx);}; //!< returns an iterator on the cells of the grid
  grid_iterator         find(const index_t &idx)       {return cellGrid.find(idx);}; //!< returns an iterator on the cells of the grid


  bool            empty() const {return totalPtCount == 0;}; //!< returns true if there is nothing stored in the grid
  size_t           size() const {return totalPtCount;}; //!< get the number of elements stored in the grid-structure
  void             insert(const T &point); //!< add an element to the ND grid. Handling of multiple occurrences depends on the "push_back" policy of C. Does not invalidate existing iterators in case C has this property
  size_t           remove(const T &point); //!< erases all elements from the ND grid corresponding to point, returns the number of elements erased. Does not invalidate existing iterators in case C has this property
  iterator        find(const T &point); //!< find a specific element, if it does not exist, return end()
  iterator         erase(iterator it); //!< erases an element from the ND grid. Does not invalidate existing iterators if C has this property
  grid_iterator    erase(grid_iterator it); //!< erases a cell from the ND grid. Does not invalidate existing iterators on other cells
  void             clear();
//  iterator        modified(iterator si); //!< should be called when an element changed, so it can be moved into another cell

  iterator        find_lazy_neighbor(const T &point, double *dist = NULL); //!< find and return neighbor (plus its distance) in the grid-cell of the point only. a closer neighbor might still exist in a neighboring cell
  const_iterator  find_lazy_neighbor(const T &point, double *dist = NULL) const; //!< find and return neighbor (plus its distance) in the grid-cell of the point only. a closer neighbor might still exist in a neighboring cell
  bool             find_lazy_neighbor(const T &point, T &neighbor) const; //!< find and return neighbor in the grid-cell of the point only. a closer neighbor might still exist in a neighboring cell
  bool             find_lazy_neighbor(const T &point, T &neighbor, double &dist) const;  //!< find and return neighbor plus its distance in the grid-cell of the point only. a closer neighbor might still exist in a neighboring cell
//  bool            find_exact_neighbor(const T &point, const double maxDist, T &neighbor);
  template <typename InsertIterator>
  void             get_cubic_neighbors(InsertIterator ii, const T &center, unsigned int radius) { //!< inserts all elements within the cubic region into the insert-iterator
    index_t ci = getIndex(center);
    cubicNeighborRecursiveHelper(ii, ci, 0, ci, radius);
  };

  size_t           getCellCount() const {return cellGrid.size();}; //!< get the number of cells that are occupied
  double           getCellExtent() const {return cellExtent;}; //!< get the edge length of each cell (equal for all dimensions)
  coordinate_t    getCellCenter(const_iterator i) const {return getCenter(i.ceIt->first);}; //!< get the coordinate of the current cell center
  template <typename InsertIterator>
  void             getCubicCellNeighbors(InsertIterator ii, const T &center, unsigned int radius) {//!< inserts grid-iterators to all cells within the cubic region into the insert-iterator
    index_t ci = getIndex(center);
    cubicCellNeighborRecursiveHelper(ii, ci, 0, ci, radius);
  };


  template <typename InsertIterator>
  unsigned int     uniform_sample(const unsigned int desiredCount, InsertIterator container, weight_func w = weight_func()) const; //!< sample uniformly over cells. returns effective count sampled (<=desiredCount)

//  friend std::ostream &operator<< (std::ostream &ostr, const index_t &idx) {
//    ostr << "( ";
//    for (unsigned int i=0; i<N; ++i) ostr << idx[i] << " ";
//    ostr << ")";
//    return ostr;
//  };
  friend std::ostream &operator<< (std::ostream &ostr, const GridND<N,T> &grid) {
    ostr << grid.cellGrid.size() << "[" << grid.totalPtCount << "]";
    for (GridConstIterator it = grid.cellGrid.begin(); it != grid.cellGrid.end(); ++it) {
      ostr << " ( ";
      for (unsigned int i=0; i<N; ++i) ostr << it->first[i] << " ";
      ostr << ")[" << it->second.size() << "]";
    }
    return ostr;
  };
};

void GridNDtest(bool verbose = true);

/*!
  \class SingleElementCell
  \brief Class for efficiently managing only 1 element per cell within the ND-Grid

  This class can be used as container for the ND-Grid and stores only 1 element
  It implements a replace-policy on inserting new elements.
  Maybe this policy could be shifted into a template parameter?
*/
template <class T>
class SingleElementCell {
public:
  typedef T* iterator; // this also allows instantiation of invalid iterator(0)
  typedef const T* const_iterator;
  SingleElementCell() : valid(false), sp() {};
  SingleElementCell(const T &sp_) : valid(true), sp(sp_) {};
  ~SingleElementCell() {};
  void operator=(SingleElementCell const &other) {valid=other.valid; sp=other.sp;};
  //bool operator==(SingleElementCell const &other) const {return ((valid=other.valid) && (sp=other.sp));};
  T& get() {if (!valid) throw std::runtime_error("SingleElementCell::get(): not valid"); return sp;};
  const T& get() const {if (!valid) throw std::runtime_error("SingleElementCell::get(): not valid"); return sp;};
  iterator begin() {if (valid) return &sp; return end();};
  iterator end() {iterator i=&sp; ++i; return i;};
  const_iterator begin() const {if (valid) return &sp; return end();};
  const_iterator end() const {const_iterator i=&sp; ++i; return i;};
  void push_back(const T &sp_) {sp=sp_; valid=true;};
  void remove(const T &sp_) { if (!valid) return; if (sp_ == sp) valid = false;};
  iterator erase(iterator it) {(void)it; valid=false; return end();}; // iterator must be valid
  void clear() {valid=false;};
  iterator find(const T &target) {if ((valid) && (sp == target)) return begin(); return end();};
  size_t size() const {return (valid) ? 1 : 0;};
private:
  bool valid;
  T sp;

  //! for storing and loading this class with boost:serialization
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    (void)version;
    ar & BOOST_SERIALIZATION_NVP(valid);
    ar & BOOST_SERIALIZATION_NVP(sp);
  }
};


#include "GridND.tcc" // include template implementations

#endif /* GRIDND_H_ */
