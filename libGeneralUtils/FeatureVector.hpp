#ifndef FEATURE_VECTOR_H
#define FEATURE_VECTOR_H

/*!
    \file   FeatureVector.h
    \brief  Provides classes representing a FeatureVector (e.g. for classification)

    The advantage here is that a base class features reference counting and
    implicit shallow copies for efficient and safe usage. The storage implementations
    are put into subclasses so that sparse, dense, or memory-chunk-pointer can be
    implemented but used with the same interface
    Additionally, tags can be associated dynamically with each feature vector

    \author  Frank Moosmann (<frank.moosmann@kit.edu>)
    \date    07.11.2005
*/
 
#include <iostream>
#include <ostream>
#include <istream>
#include <stdexcept>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <boost/typeof/typeof.hpp>

#include "Tag.hpp"

typedef float FeatureData;

/*!
 * \class FeatureVectorDataBase
 * \brief This class serves as abstract base class for all storage implementation that can be used together with FeatureVector
 * \note This class should NEVER be used directly, but always within FeatureVector. The only way a storage container
 *        should be used is within the constructor of FeatureVector, e.g.: FeatureVector fv(new FeatureVectorDataBaseDerivedClass());
 */
class FeatureVectorDataBase {
  public:
    FeatureVectorDataBase();
    virtual ~FeatureVectorDataBase();
    virtual FeatureVectorDataBase* clone() const; //!< make a deep copy of everything
    void cloneTags(FeatureVectorDataBase* target) const; //!< deep-copies the tags into the target

    // the following functions must be implemented by derived classes:
    virtual FeatureVectorDataBase* cloneData() const = 0; //!< make a deep copy of only the data, not the tags
    virtual void          clear() = 0; //!< Clears all stored data. If not possible a std::logic_error must be thrown
    virtual unsigned int  size() const = 0; //!< Returns the size of the data
    virtual FeatureData&  getElement(unsigned int idx) = 0; //!< Get one specific data element. Index must be in range 0..size()-1
    virtual FeatureData   getElement(unsigned int idx) const = 0; //!< Get one specific data element. Index must be in range 0..size()-1
    virtual void          setElement(const unsigned int index, const FeatureData elem) = 0; //!< Set one specific data element. if index>=size(), size must be automatically increased. If not possible a std::logic_error must be thrown
    virtual void          push_back(const FeatureData elem) = 0; //!< Appends the data with the given element. If not possible a std::logic_error must be thrown

//    typedef std::vector<FeatureData>::iterator iterator; //!< Provides an iterator to stored data
//    iterator                  begin(); //!< Returns an iterator pointing to the beginning of the data
//    iterator                  end(); //!< Returns an iterator pointing to the end of the data

  private:
    friend class FeatureVector;
    typedef std::map<TagType, Tag*> TagMap;
    TagMap                tags; //!< associated tags. access handled/implemented in FeatureVector
    unsigned int          refCount; //!< reference counting. access handled/implemented in FeatureVector
};


/*!
 * \class FeatureVector
 * \brief This class can be used for classification etc. to store vector data
 *
 * It features reference counting, i.e. assigning or copying a FeatureVector does
 * not copy the complete data!
 */
class FeatureVector {
  public:
    FeatureVector(); //!< default constructor. will implicitly create a FeatureVectorDenseSTLData as storage
    FeatureVector(const FeatureVector &other); //!< copy constructor, creates a shallow copy
    FeatureVector(FeatureVectorDataBase *data); //!< constructor using a new explicit storage. will overtake memory ownership of data. Usage: FeatureVector fv(new myFeatureVectorData());
    template <class FeatureVectorDataBaseT>
    FeatureVector() {fvdata = new FeatureVectorDataBaseT(); fvdata->refCount = 1;}; //!< constructor that implicitly creates a new storage. Usage: FeatureVector fv<myFeatureVectorData>();
    ~FeatureVector();
    
    FeatureVector clone() const {return FeatureVector(fvdata->clone());}; //!< makes a deep copy of the current FeatureVector. keeps the storage type
    template <class DimIndexIterator>
    FeatureVector clone(DimIndexIterator begin, DimIndexIterator end) const; //!< makes a deep copy of the current FeatureVector, but only the specified dimensions (must be in ascending order) if begin!=end. changes storage type to DenseSTL
    FeatureVector cloneData() const {return FeatureVector(fvdata->cloneData());}; //!< makes a deep copy of the current FeatureVector but does not copy tags. keeps the storage type
    FeatureVector& operator= (const FeatureVector &other);
    bool operator== (const FeatureVector &other) const;
    bool operator!= (const FeatureVector &other) const {return !(*this == other);};
    FeatureVector& operator+= (const FeatureVector &other);
    FeatureVector& operator/ (const float value);

//    friend std::ostream &operator<<(std::ostream &ostr, const FeatureVector &f);

    // TODO: change all non-const functions below to create a copy before changing the vector in case reference count > 1
    // (ImageMagick++ behavior) This avoids unexpected changes in one vector when another vector is changed
    void               clear() {fvdata->clear();};
    unsigned int       size() const {return getSize();}; //!< Returns the size of the data
    unsigned int       getSize() const {return fvdata->size();}; //!< Returns the size of the data
    FeatureData&       operator()(unsigned int idx) {return fvdata->getElement(idx);}; //!< same as getElement
    FeatureData        operator()(unsigned int idx) const {return fvdata->getElement(idx);}; //!< same as getElement
    FeatureData&       operator[](unsigned int idx) {return fvdata->getElement(idx);}; //!< same as getElement
    FeatureData        operator[](unsigned int idx) const {return fvdata->getElement(idx);}; //!< same as getElement
    FeatureData        getElement(unsigned int idx) const {return fvdata->getElement(idx);}; //!< Get one specific data element. Index must be in range 0..getSize()-1
    void               push_back(const FeatureData elem) {fvdata->push_back(elem);}; //!< Appends the data with the given element
    void               addElement(const FeatureData elem) {push_back(elem);}; //!< Appends the data with the given element
    void               setElement(const int index, const FeatureData elem) {fvdata->setElement(index,elem);}; //!< Set one specific data element.
//    typedef std::vector<FeatureData>::iterator iterator; //!< Provides an iterator to stored data
//    iterator           begin(); //!< Returns an iterator pointing to the beginning of the data
//    iterator           end(); //!< Returns an iterator pointing to the end of the data

    bool                                   hasTag(const TagType type) const; //!< Returns true if a Tag with the given Type is stored along with the data.
    template< class CustomTag > bool       hasTag() const; //!< Returns true if a Tag with the given Type is stored along with the data.
    template< class CustomTag > void        setTag(const CustomTag &tag); //!< Provides the possibility to store Tags along with the data. Only one Tag per TagType can be stored. Tag will be copied
    template< class CustomTag > CustomTag& getTag() const; //!< Returns the Tag with the given Type stored along with the data. \throws runtime_error if this TagType doesn't exist
    Tag*                                    getTag(const TagType type) const; //!< Returns the Tag with the given Type. Does not throw but return a NULL-pointer if TagType doesn't exist

    static float         cityblockDistance(const FeatureVector &f1, const FeatureVector &f2);
    static float         euclideanDistance(const FeatureVector &f1, const FeatureVector &f2);
    template <class InputIterator, class InsertIterator>
    static void          extractDimension(InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, InsertIterator ins, const unsigned int dimNb, bool checkEntryValidity = false); //!< extracts the specified dimension and inserts all elements via the insert iterator
    template <class InputIterator>
    static FeatureData   compEmpMean (InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, const unsigned int dimNb, bool checkEntryValidity = false); //!< returns the empirical mean of the FeatureVectors at dimension dimNb
    template <class InputIterator>
    static FeatureVector compEmpMean (InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, bool checkEntryValidity = false); //!< returns the mean of the FeatureVectors
    template <class InputIterator>
    static FeatureData   compEmpStdDeviation (InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, const unsigned int dimNb, const FeatureData mean, bool checkEntryValidity = false); //!< returns the standard deviation of the FeatureVectors at dimension dimNb

  private:
    void deinit();
    FeatureVectorDataBase *fvdata; //!< pointer to data. functions as shared pointer, if FeatureVector is copied
};


/*!
 * \class FeatureVectorDenseSTLData
 * \brief This is one specific implementation of storage for dense data
 */
class FeatureVectorDenseSTLData : public FeatureVectorDataBase {
  public:
    FeatureVectorDenseSTLData();
    FeatureVectorDenseSTLData(const FeatureVector &other);
    FeatureVectorDenseSTLData(unsigned int size, FeatureData initValue);
    template<typename floatT>
    FeatureVectorDenseSTLData(const floatT *arr, unsigned int size); //!< copies the data given at some memory chunk
    virtual ~FeatureVectorDenseSTLData();

    virtual FeatureVectorDataBase* cloneData() const;
    virtual void          clear() {data.clear();};
    virtual unsigned int  size() const {return data.size();};
    virtual FeatureData&  getElement(unsigned int idx) {return data[idx];};
    virtual FeatureData   getElement(unsigned int idx) const {return data[idx];};
    virtual void          setElement(const unsigned int index, const FeatureData elem) {if (index>=size()) data.resize(index+1); data[index] = elem;};
    virtual void          push_back(const FeatureData elem) {data.push_back(elem);};

//    typedef std::vector<FeatureData>::iterator iterator; //!< Provides an iterator to stored data
//    iterator           begin(); //!< Returns an iterator pointing to the beginning of the data
//    iterator           end(); //!< Returns an iterator pointing to the end of the data

  private:
    std::vector<FeatureData> data; //!< data storage
};


/*!
 * \class FeatureVectorPointerData
 * \brief This is an implementation of storage that refers to data in some (external) memory block
 */
template <typename dataT>
class FeatureVectorPointerData : public FeatureVectorDataBase {
  public:
    FeatureVectorPointerData(dataT *arr, unsigned int size); //!< does not copy data but only store the pointer
    virtual ~FeatureVectorPointerData();

    virtual FeatureVectorDataBase* cloneData() const;
    virtual void          clear() {throw std::logic_error("FeatureVectorPointerData does not allow to clear data");};
    virtual unsigned int  size() const {return datasize;};
    virtual FeatureData&  getElement(unsigned int idx) {return data[idx];};
    virtual FeatureData   getElement(unsigned int idx) const {return data[idx];};
    virtual void          setElement(const unsigned int index, const FeatureData elem) {if (index>=datasize) throw std::out_of_range("FeatureVectorPointerData::setElement"); data[index] = elem;};
    virtual void          push_back(const FeatureData elem) {(void)elem;throw std::logic_error("FeatureVectorPointerData does not allow to push_back");};

//    typedef std::vector<FeatureData>::iterator iterator; //!< Provides an iterator to stored data
//    iterator           begin(); //!< Returns an iterator pointing to the beginning of the data
//    iterator           end(); //!< Returns an iterator pointing to the end of the data

  private:
    unsigned int datasize;
    dataT *data; //!< external data storage
};

/*!
 * \class FeatureVectorStepPointerData
 * \brief This is an implementation of storage that refers to regular data in some
 *         (external) memory block (it can e.g. be used to refer to a column of a row-wise stored matrix)
 */
template <typename dataT>
class FeatureVectorStepPointerData : public FeatureVectorDataBase {
  public:
  FeatureVectorStepPointerData(dataT *arr, unsigned int stepSize, unsigned int size); //!< does not copy data but only store the pointer
    virtual ~FeatureVectorStepPointerData();

    virtual FeatureVectorDataBase* cloneData() const;
    virtual void          clear() {throw std::logic_error("FeatureVectorStepPointerData does not allow to clear data");};
    virtual unsigned int  size() const {return datasize;};
    virtual FeatureData&  getElement(unsigned int idx) {return data[idx*stepsize];};
    virtual FeatureData   getElement(unsigned int idx) const {return data[idx*stepsize];};
    virtual void          setElement(const unsigned int index, const FeatureData elem) {if (index>=datasize) throw std::out_of_range("FeatureVectorPointerData::setElement"); data[index*stepsize] = elem;};
    virtual void          push_back(const FeatureData elem) {(void)elem;throw std::logic_error("FeatureVectorStepPointerData does not allow to push_back");};

//    typedef std::vector<FeatureData>::iterator iterator; //!< Provides an iterator to stored data
//    iterator           begin(); //!< Returns an iterator pointing to the beginning of the data
//    iterator           end(); //!< Returns an iterator pointing to the end of the data

  private:
    unsigned int datasize;
    unsigned int stepsize;
    dataT *data; //!< external data storage
};



#include "FeatureVector.tcc"

#endif /*FEATURE_VECTOR_H*/
