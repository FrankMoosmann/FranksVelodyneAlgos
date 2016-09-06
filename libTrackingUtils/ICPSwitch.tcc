
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////        template implementation       //////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>

#include "HomogeneousTransformationMatrix.hpp"
#include "TrackingUtils.hpp"


template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPSwitch<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
ICPSwitch(DVectorInputIterator ptBegin_, DVectorInputIterator ptEnd_,
							DVectorInputIterator nBegin_, DVectorInputIterator nEnd_,
							DInputIterator wBegin_, DInputIterator wEnd_,
							DInputIterator nConfBegin_, DInputIterator nConfEnd_,
							DInputIterator nStdDevBegin_, DInputIterator nStdDevEnd_,
              DMatrixInputIterator pCovarBegin, DMatrixInputIterator pCovarEnd,
              DMatrixInputIterator nCovarBegin, DMatrixInputIterator nCovarEnd,
							corrFuncT corrFunc, EnergyFunction &energy1_, EnergyFunction &energy2_,
							unsigned int nbIterSwitch_, switchFunc sFunc_, bool sFuncInitVal_,
							unsigned int nCorr, double regCoeff_, EstimationMode estMode_)
	:ICPLin(ptBegin_, ptEnd_, nBegin_, nEnd_, wBegin_, wEnd_, nConfBegin_, nConfEnd_, nStdDevBegin_, nStdDevEnd_,
	    pCovarBegin, pCovarEnd, nCovarBegin, nCovarEnd,
	    corrFunc, energy1_, nCorr, regCoeff_, estMode_)
	,energy1(&energy1_)
	,energy2(&energy2_)
	,nbIterSwitch(nbIterSwitch_)
	,sFunc(sFunc_)
	,sFuncInitVal(sFuncInitVal_)
	,currIter(0)
{
	if (!sFunc.empty()) sFunc(sFuncInitVal);
	ICPBase::reset();
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPSwitch<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
~ICPSwitch()
{
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
void ICP::ICPSwitch<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
reset(const DMatrix &rot, const DVector &trans)
{
	currIter = 0;
	ICPLin::energy = energy1;
	if (!sFunc.empty()) sFunc(sFuncInitVal);
	ICPLin::reset(rot, trans);
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
bool ICP::ICPSwitch<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
iterate()
{ 
	if (++currIter <= nbIterSwitch) {
		ICPLin::energy = energy1;
		if (!sFunc.empty()) sFunc(sFuncInitVal);
	} else {
	  if (currIter == nbIterSwitch+1) std::cout << " SWITCH" << std::flush;
		ICPLin::energy = energy2;
		if (!sFunc.empty()) sFunc(!sFuncInitVal);
	}
	bool success = ICPLin::iterate();
	ICPLin::energy = energy1; // for manual covariance/error calculation always use first energy

	return success;
}

