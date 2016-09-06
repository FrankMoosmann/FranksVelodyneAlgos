
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////        template implementation       //////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include <cfloat>
#include <stdexcept>
#include <boost/bind.hpp>
#include <boost/static_assert.hpp>

#include "HomogeneousTransformationMatrix.hpp"


//template <typename DVectorInputIterator, typename SearchResultsT>
//ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::
//ICPBesl(ICP::cpSearch corrFunc, DVectorInputIterator ptBegin_, DVectorInputIterator ptEnd_, DVectorInputIterator nBegin_, DVectorInputIterator nEnd_)
//	:searchCorrespondence(boost::bind(&cpSearch2cpwSearch, _1, _2, _3, _4, _5, corrFunc))
//	,useWeights(false)
//	,ptBegin(ptBegin_)
//	,ptEnd(ptEnd_)
//	,nBegin(nBegin_)
//	,nEnd(nEnd_)
//{
//  init();
//}

template <typename DVectorInputIterator, typename SearchResultsT>
ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::
ICPBesl(DVectorInputIterator ptBegin_, DVectorInputIterator ptEnd_, DVectorInputIterator nBegin_, DVectorInputIterator nEnd_, corrFuncT corrFunc, unsigned int nCorr)
	: ICPBase()
	,searchCorrespondence(corrFunc)
	,nbCorresp(nCorr)
	,useWeights(true)
	,ptBegin(ptBegin_)
	,ptEnd(ptEnd_)
	,nBegin(nBegin_)
	,nEnd(nEnd_)
{
  init();
}


template <typename DVectorInputIterator, typename SearchResultsT>
ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::
~ICPBesl()
{
}

template <typename DVectorInputIterator, typename SearchResultsT>
void ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::init()
{
	ptCount = 0;
	for (DVectorInputIterator tmp=ptBegin; tmp!=ptEnd; ++tmp)
		++ptCount;
	transformedModel.resize(3,ptCount,false);
	worldCorresp.resize(3,ptCount*nbCorresp,false);
  modelCenter.resize(3,false);
  worldCenter.resize(3,false);
  centeredModel.resize(3,ptCount*nbCorresp,false);
  centeredWorld.resize(3,ptCount*nbCorresp,false);
  A.resize(3,false);
  U.resize(3,3,false);
  Vt.resize(3,3,false);
  if (useWeights)
  	weights = DVector(ptCount*nbCorresp, 1.0); // weights
	ICPBase::reset();
}

//template <typename DVectorInputIterator, typename SearchResultsT>
//void ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::reset(const DMatrix &rot, const DVector &trans)
//{
//	ICPBase::reset(rot, trans);
//}

template <typename DVectorInputIterator, typename SearchResultsT>
bool ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::iterate()
{ 
	return (FindCorrespondences() && ComputeTransformation());
}


template <typename DVectorInputIterator, typename SearchResultsT>
bool ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::FindCorrespondences()
{
  /* reads following member variables:
   * - ptBegin, ptEnd  	iterators over the template points
   * - nBegin, nEnd   	iterators over the template normals
   * - R    						current rotation matrix
   * - t    						current translation vector
   * 
   * writes following member variables:
   * - worldCorresp    	matched scene points (correspondences of template points)
   * - weights    			weight matrix
   * - worldCenter   		weighted center of matched scene points
   * - modelCenter   		weighted center of template
   * - cumWeight 				cumulative weight
   * - cumWeightDist  	cumulative weighted distance of correspondences
   *
   * writes following inherited member variables:
   * - avgDist          average distance of correspondences
   * - wAvgDist         weighted average distance of correspondences
   * - avgErr           average error before minimization
   * - wAvgErr          weighted average error before minimization
   */
	namespace ublas=boost::numeric::ublas;
	using namespace std;
	modelCenter = ublas::zero_vector<double>(3); // Weighted center of template
	worldCenter = ublas::zero_vector<double>(3); // Weighted center of matched scene points
	cumWeight = 0.0;
  avgDist = 0.0;
  wAvgDist = 0.0;
  avgErr = 0.0;
  wAvgErr = 0.0;
  unsigned int validCorrespCnt = 0;

  // for each point in the template
  DVectorInputIterator pi = ptBegin;
  DVectorInputIterator ni = nBegin;
  unsigned int i = 0;
  bool oneValid = false;
	DVector SCorr(3); // temporary: next neighbor (correspondence), crossproduct
	double weight, dist;
//	cout << endl << "searching " << ptCount << " points";
	while ((pi != ptEnd) && (ni != nEnd))
	{
		// apply current R+t
	  DVector pTmp = prod(R,*pi) + t;
    DVector nTmp = prod(R,*ni);
    transformedModel(0,i)=pTmp(0);
    transformedModel(1,i)=pTmp(1);
    transformedModel(2,i)=pTmp(2);
		
    // search nearest neighbor
		SearchResultsT sr = searchCorrespondence(pTmp, nTmp, 1.0, nbCorresp); // will not throw
		for (unsigned int c=0; c<nbCorresp; ++c) {
			unsigned int i2 = i*nbCorresp + c;
			try {
				SCorr = sr.getCorrespPoint();
//				nNN = sr.getCorrespNormal();
				dist = sr.getDist();
				weight = useWeights ? sr.getWeight() : 1.0;
				sr.next();
				cout << "." << flush;
			} catch (exception &e) {
				cout << "E" << flush;
//	      cerr << endl << e.what();
	      SCorr = pTmp;
				dist = 0.0;
				weight = 0.0;
			}
	    // store resulting correspondence in matrix "Matched Scene Points"
			worldCorresp(0,i2) = SCorr(0);
			worldCorresp(1,i2) = SCorr(1);
			worldCorresp(2,i2) = SCorr(2);
	    if (useWeights) {
				weights(i2) = weight;
	    }
			if (weight != 0.0) {
				++validCorrespCnt;
				cumWeight += weight;
				avgDist += dist; // will be normalized after for-loop
				avgErr += dist*dist;
				wAvgDist += weight*dist;
				wAvgErr += weight*dist*dist;

				// calculate weighted center of template and matched scene points
				worldCenter += SCorr*weight;
				modelCenter += pTmp*weight;
				oneValid = true;
			}
		} // end for each correspondence
    ++pi; // next point
    ++ni; // next normal
    ++i; // increase matrix index
	} // end: for each model point
  if (cumWeight > 0.01) {
  	worldCenter = worldCenter/cumWeight;
  	modelCenter = modelCenter/cumWeight;	
    avgDist /= (double)validCorrespCnt;
		wAvgDist /= cumWeight;
    avgErr /= (double)validCorrespCnt;
		wAvgErr /= cumWeight;
  } else {
    avgDist = DBL_MAX;
    wAvgDist = DBL_MAX;
    avgErr = DBL_MAX;
    wAvgErr = DBL_MAX;
  }
//	cout << endl << "world center = " << worldCenter;
//	cout << endl << "model center = " << modelCenter;
//	cout << endl << "cumWeight = " << cumWeight;
  return oneValid;
}

template <typename DVectorInputIterator, typename SearchResultsT>
bool ICP::ICPBesl<DVectorInputIterator, SearchResultsT>::ComputeTransformation()
{	
  /* reads following member variables:
   * - worldCorresp		matched scene points (correpondences of template points)
   * - weights    		weight matrix (update)
   * - worldCenter   	weighted center of matched scene points
   * - modelCenter   	weighted center of template
   * - cumWeight   		cumulative weight
   * - cumWeightDist  cumulative weighted distance of correspondences
   * 
   * writes following member variables:
   * - R    					current rotation matrix
   * - t    					current translation vector
   */
	namespace ublas=boost::numeric::ublas;
	using namespace std;
   
  // directly return if not enough correspondences are available
  if (cumWeight < 0.01) {
    return false;
  }

  // center template and matched scene points
	for (unsigned int i=0; i<ptCount;++i) {
		for (unsigned int c=0; c<nbCorresp;++c) {
			for (unsigned int d=0; d<3;++d) {
				centeredModel(d,i*nbCorresp+c) = transformedModel(d,i) - modelCenter(d);
				centeredWorld(d,i*nbCorresp+c) = worldCorresp(d,i*nbCorresp+c) - worldCenter(d);
			}
		}
	}

  // calculate weighted covariance matrix K
	tmp=trans(centeredModel);
	if (useWeights) {
		//tmp=prod(weights,tmp);
		assert(weights.size() == tmp.size1());
		for (unsigned int i=0; i<tmp.size1(); ++i) {
			for (unsigned int j=0; j<tmp.size2(); ++j) {
				tmp(i,j) = tmp(i,j) * weights(i);
			}
		}
	}
	K=prod(centeredWorld,tmp);
	//K/=cumWeight;
  // do Singular Value Decomposition
	boost::numeric::bindings::lapack::gesvd ('A', 'A', K, A, U, Vt);

  // calculate R+t
	DMatrix Ri=prod(U,Vt);
	double det = Ri(0,0)*Ri(1,1)*Ri(2,2) + Ri(0,1)*Ri(1,2)*Ri(2,0) + Ri(0,2)*Ri(1,0)*Ri(2,1)
			       - Ri(0,2)*Ri(1,1)*Ri(2,0) - Ri(0,1)*Ri(1,0)*Ri(2,2) - Ri(0,0)*Ri(1,2)*Ri(2,1);
	if (det < 0) {
	  // calculate Ri=prod(U, B, Vt) with B = determinant-modified Id3
    cout << " X " << flush;
		DMatrix B(ublas::identity_matrix<double>(3));
		B(2,2)=det;
		Ri=prod(B,Vt);
		Ri=prod(U,Ri);
	}

//  Ri = trans(Ri); // transpose (i.e. inverse)
  DVector ti = worldCenter-prod(Ri,modelCenter);

//Activate this code to disable rotation estimation:
//  DMatrix Ri = ublas::identity_matrix<double>(3);
//  DVector ti = worldCenter-modelCenter;

//  cout << endl << "    ti=" << ti;
//  cout << endl << "    Ri=" << Ri;
	R = ublas::prod(Ri,R);
	t = ublas::prod(Ri,t) + ti;
//  cout << endl << "    t=" << t;
//  cout << endl << "    R=" << R;
	return true;
}


