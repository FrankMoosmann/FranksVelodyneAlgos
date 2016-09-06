/*
 * FeatureVector.tcc
 *
 * Template implementation for FeatureVector.h
 */
#include "TagInvalidEntries.hpp"
#include <exception>
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////        Class  FeatureVector           /////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

template< class CustomTag >
void FeatureVector::setTag(const CustomTag &tag) {
	CustomTag *tagcopyp = new CustomTag(tag); // we manage memory! also need pointer to add to list
	TagType type = tagcopyp->getType();
	BOOST_AUTO(i,fvdata->tags.find(type));
	if (i != fvdata->tags.end()) {
		delete i->second; //overwriting stored tag. so free it first
	}
	fvdata->tags[type] = tagcopyp;
}

template< class CustomTag >
CustomTag& FeatureVector::getTag() const {
	TagType type = CustomTag::Type();
	BOOST_AUTO(i,fvdata->tags.find(type));
	if (i != fvdata->tags.end()) {
		CustomTag *ct = dynamic_cast<CustomTag*>(i->second);
		return *ct;
	} else {
		std::string msg = "FeatureVector::getTag(): Type "+CustomTag::TypeStr()+" not available";
		throw std::runtime_error(msg);
	}
}

template< class CustomTag >
bool FeatureVector::hasTag() const {
	TagType type = CustomTag::Type();
	BOOST_AUTO(i,fvdata->tags.find(type));
	return (i!=fvdata->tags.end());
}

template <class DimIndexIterator>
FeatureVector FeatureVector::clone(DimIndexIterator diBegin, DimIndexIterator diEnd) const
{
  if (diBegin == diEnd)
    return clone(); // more efficient, will keep storage type
  FeatureVectorDataBase *data = new FeatureVectorDenseSTLData();
  fvdata->cloneTags(data); // copy tags
  for (unsigned int i=0; i<fvdata->size(); ++i) {
    if ((diBegin != diEnd) && (i == *diBegin)) { // copy dimension
      data->push_back(fvdata->getElement(i));
      ++diBegin;
    }
  }

  return FeatureVector(data);
}

template <class InputIterator, class InsertIterator>
void FeatureVector::extractDimension(InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, InsertIterator ins, const unsigned int dimNb, bool checkEntryValidity)
{
  for (InputIterator fVec = FeatureVectorBegin; fVec != FeatureVectorEnd; ++fVec) {
    if (checkEntryValidity) {
      try{
        TagInvalidEntries tag = fVec->getTag<TagInvalidEntries>();
        if (tag.isValid(dimNb))
          *ins++ = fVec->getElement(dimNb);
      } catch(std::exception& e) {
      }
    } else {
      *ins++ = fVec->getElement(dimNb);
    }
  }
}

template <class InputIterator>
FeatureData FeatureVector::compEmpMean (InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, const unsigned int dimNb, bool checkEntryValidity)
{
//	if(dimNb == 16 ||dimNb == 15 || dimNb == 17)
//		std::cout << std::endl << " compEmpMean dimNb: " << dimNb << std::endl;
  FeatureData sum = 0;
  FeatureData count = 0;
  while (FeatureVectorBegin != FeatureVectorEnd) {
  	bool use = true;
  	if (checkEntryValidity){
  		try{
  			use = false;
  			TagInvalidEntries tag = FeatureVectorBegin->getTag<TagInvalidEntries>();
  			if(tag.isValid(dimNb)){
//  				if(dimNb == 16 ||dimNb == 15 || dimNb == 17)
//  					std::cout << FeatureVectorBegin->getElement(dimNb)<< ", " << std::flush;
  				sum += FeatureVectorBegin->getElement(dimNb);
  				count += 1;
  			}
  		}catch(std::exception& e){
  			use = true;
  		}
  	}
  	if (use) {
			sum += FeatureVectorBegin->getElement(dimNb);
			count += 1;
  	}
    ++FeatureVectorBegin;
  }
  //cout << "  attrNo: " << attrNo << "  sum: " << sum << "  count: " << count;
  if (count > 0)
    return sum/count;
  else
    return 0.0;
}

template <class InputIterator>
FeatureVector FeatureVector::compEmpMean (InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, bool checkEntryValidity)
{
	unsigned int nbDim = FeatureVectorBegin->size();
  FeatureData *sum = new FeatureData[nbDim];
  for (unsigned int i=0; i<nbDim; ++i)
  	sum[i] = 0;
  std::vector<FeatureData> count(nbDim,0.0);
  while (FeatureVectorBegin != FeatureVectorEnd) {
  	bool use = true;
  	if (checkEntryValidity) {
  		try{
  			TagInvalidEntries tag = FeatureVectorBegin->getTag<TagInvalidEntries>();
  			for (unsigned int i=0; i<nbDim; ++i){
  				if(tag.isValid(i)){
  				sum[i] += FeatureVectorBegin->getElement(i);
  				count.at(i) += 1;
  				}
  			}
  			use = false;
  		} catch(std::exception& e) {
  			use = true;
  		}
  	}
  	if( use ){
  		for (unsigned int i=0; i<nbDim; ++i){
  			sum[i] += FeatureVectorBegin->getElement(i);
  			count.at(i) += 1;
  		}
  	}
    ++FeatureVectorBegin;
  }
  //cout << "  attrNo: " << attrNo << "  sum: " << sum << "  count: " << count;
  for (unsigned int i=0; i<nbDim; ++i){
		if (count.at(i) > 0) {
			sum[i] /= count.at(i);
		}else{
			sum[i] = 0.0;
		}
  }
  
  delete[] sum;
  
  return FeatureVector(new FeatureVectorDenseSTLData(sum,nbDim));
}

template <class InputIterator>
FeatureData FeatureVector::compEmpStdDeviation (InputIterator FeatureVectorBegin, InputIterator FeatureVectorEnd, const unsigned int dimNb, const FeatureData mean, bool checkEntryValidity)
{
  FeatureData sum = 0.0;
  FeatureData count = 0;
  while (FeatureVectorBegin != FeatureVectorEnd) {
  	bool use = true;
  	if (checkEntryValidity){
  		try{
  			TagInvalidEntries tag = FeatureVectorBegin->getTag<TagInvalidEntries>();
  			if(tag.isValid(dimNb)){
  				FeatureData fVal = FeatureVectorBegin->getElement(dimNb);
  				sum += pow(fVal - mean,2);
  				count += 1.0;
  			}
  		}catch(std::exception& e){
  			use = true;
  		}
  	}
  	if( use ) {
  		FeatureData fVal = FeatureVectorBegin->getElement(dimNb);
  		sum += pow(fVal - mean,2);
			count += 1.0;
  	}
    ++FeatureVectorBegin;
  }
  if (count > 1)
    return sqrt(sum/(count-1));
  else
    return 0.0;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////   Class  FeatureVectorDenseSTLData    /////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

template<typename floatT>
FeatureVectorDenseSTLData::FeatureVectorDenseSTLData(const floatT *arr, unsigned int size)
	: FeatureVectorDataBase()
{
	data.reserve(size);
  for (unsigned int i=0; i<size; ++i)
		data.push_back(arr[i]);
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////   Class  FeatureVectorPointerData     /////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

template <typename dataT>
FeatureVectorPointerData<dataT>::FeatureVectorPointerData(dataT *arr, unsigned int size)
	: FeatureVectorDataBase()
{
	datasize = size;
	data = arr;
}

template <typename dataT>
FeatureVectorPointerData<dataT>::~FeatureVectorPointerData()
{
}

template <typename dataT>
FeatureVectorDataBase* FeatureVectorPointerData<dataT>::cloneData() const
{
	throw std::logic_error("FeatureVectorPointerData does not allow to clone");
	return NULL;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////  Class  FeatureVectorStepPointerData  /////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

template <typename dataT>
FeatureVectorStepPointerData<dataT>::FeatureVectorStepPointerData(dataT *arr, unsigned int step, unsigned int size)
	: FeatureVectorDataBase()
{
	datasize = size;
	stepsize = step;
	data = arr;
}

template <typename dataT>
FeatureVectorStepPointerData<dataT>::~FeatureVectorStepPointerData()
{
}

template <typename dataT>
FeatureVectorDataBase* FeatureVectorStepPointerData<dataT>::cloneData() const
{
	throw std::logic_error("FeatureVectorStepPointerData does not allow to clone");
	return NULL;
}

