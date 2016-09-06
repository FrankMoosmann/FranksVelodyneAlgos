
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////        template implementation       //////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include <cfloat>
#include <cmath>
#include <stdexcept>
#include <boost/bind.hpp>
#include <boost/static_assert.hpp>
#include <cvm.h>

#include "HomogeneousTransformationMatrix.hpp"


template <typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPPlanePlane<DVectorInputIterator, DInputIterator, SearchResultsT>::
ICPPlanePlane(DVectorInputIterator ptBegin_, DVectorInputIterator ptEnd_,
							DVectorInputIterator nBegin_, DVectorInputIterator nEnd_,
							DInputIterator nConfBegin_, DInputIterator nConfEnd_,
							DInputIterator nStdDevBegin_, DInputIterator nStdDevEnd_,
							corrFuncT corrFunc, unsigned int nCorr, double regCoeff_)
	:ICPBase()
	,ptBegin(ptBegin_)
	,ptEnd(ptEnd_)
	,nBegin(nBegin_)
	,nEnd(nEnd_)
	,nConfBegin(nConfBegin_)
	,nConfEnd(nConfEnd_)
	,nStdDevBegin(nStdDevBegin_)
	,nStdDevEnd(nStdDevEnd_)
	,searchCorrespondence(corrFunc)
	,nbCorresp(nCorr)
	,regCoeff(regCoeff_)
	,useWeights(true)
  ,U(3,3)
  ,Vt(3,3)
  ,Dv(3)
	,AtA(3,3)
	,Atb(3)
	,Ri(3,3)
	,ti(3)
{
	nbModelPts = 0;
	DVectorInputIterator pTmp=ptBegin;
	for (DVectorInputIterator nTmp=nBegin; nTmp!=nEnd; ++nTmp,++pTmp)
		++nbModelPts;
	if (pTmp != ptEnd)
		throw std::out_of_range("ICP-Constructor: number of points != number of normals"); // loop stopped on end of Normal-Iterator, so Point-Iterator should be "End" as well
	transModelPts.resize(3,nbModelPts,false);
	transModelNrm.resize(3,nbModelPts,false);
	CorrPts.resize(3,nbModelPts*nbCorresp,false);
	CorrNrm.resize(3,nbModelPts*nbCorresp,false);
	CorrWeights = DVector(nbModelPts*nbCorresp, 1.0); // assume equal full weight
	CorrNConf = DVector(nbModelPts*nbCorresp, 1.0); // assume equal full confidence
	CorrNStdDevRAD = DVector(nbModelPts*nbCorresp, 0.0); // assume equal zero uncertainty
	A.resize(nbModelPts*nbCorresp,3,false);
	At.resize(3,nbModelPts*nbCorresp,false);
	b.resize(nbModelPts*nbCorresp);
  ICPBase::reset();
}

template <typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPPlanePlane<DVectorInputIterator, DInputIterator, SearchResultsT>::
~ICPPlanePlane()
{
}


template <typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
bool ICP::ICPPlanePlane<DVectorInputIterator, DInputIterator, SearchResultsT>::
iterate()
{ 
	using namespace matrixTools;
	using namespace std;

	// matrix proxies are references, so no memory is copied
	typedef ublas::matrix_column< DCMatrix > DCMatrixCol;
	typedef ublas::matrix_range< DCMatrix > DCMatrixRange;
	typedef ublas::matrix_column< DCMatrixRange > DCMatrixRangeCol;
	typedef ublas::matrix_vector_range< DCMatrix > DCMatrixVRange;

	if (nbModelPts < 3) return false;

	/////////////////////////////////////////////////////////////
  //               search correspondences                    //
  /////////////////////////////////////////////////////////////
  avgDist=0.0;
	wAvgDist=0.0;
  avgErr=0.0;
	wAvgErr=0.0;
	cumWeight = 0.0;
  DVectorInputIterator pi = ptBegin;
  DVectorInputIterator ni = nBegin;
  unsigned int i = 0;
  unsigned int i2;
  DVector pt(3), nr(3);
	double dist,err,weight;
	while ((pi != ptEnd) && (ni != nEnd)) // for each model point
	{
		pt = prod(R,*pi)+t;
		nr = prod(R,*ni);
		SearchResultsT sr = searchCorrespondence(pt, nr, nbCorresp);
		DCMatrixCol(transModelPts,i) = pt;
		DCMatrixCol(transModelNrm,i) = nr;
		for (unsigned int c=0; c<nbCorresp; ++c) { // for each correspondence of the current model point
			i2 = i*nbCorresp+c;
			try {
				DCMatrixCol(CorrPts,i2) = sr.getCorrespPoint();
				DCMatrixCol(CorrNrm,i2) = sr.getCorrespNormal();
				weight = useWeights ? sr.getWeight() : 1.0;
				CorrNConf(i2) = sr.getNormConf();
				CorrNStdDevRAD(i2) = sr.getNormStdDevRAD();
				dist = sr.getDist();
				err = (1.0-fabs(ublas::inner_prod(nr,sr.getCorrespNormal()))) + pow(ublas::inner_prod(pt-sr.getCorrespPoint(),sr.getCorrespNormal()),2);
				sr.next();
			} catch (const std::exception &e) {
				DCMatrixCol(CorrPts,i2) = pt;
				DCMatrixCol(CorrNrm,i2) = nr;
				CorrNConf(i2) = 0.0;
				CorrNStdDevRAD(i2) = M_PI/4;
				weight = 0.0;
				dist = 0.0;
				err = 0.0;
			}
			if (useWeights)
				CorrWeights(i2) = weight;
			avgDist += dist;
			avgErr += err;
			wAvgDist += weight*dist;
			wAvgErr += weight*err;
			cumWeight += weight;
		} // end: for each nearest neighbor to the current model point
    ++pi;
    ++ni;
    ++i;
	} // end: for each model point
	if (i != nbModelPts)
		throw std::logic_error("ICPLinearized::establishCorrespondences: number of points has changed!");
	avgDist /= (double)(nbModelPts*nbCorresp);
	avgErr /= (double)(nbModelPts*nbCorresp);
	wAvgErr /= cumWeight;
	wAvgDist /= cumWeight;

	/////////////////////////////////////////////////////////////
  //                   compute rotation                      //
  /////////////////////////////////////////////////////////////
	H = DCMatrix(3,3,0.0);
	for (unsigned int i=0; i<nbModelPts; ++i) {
		for (unsigned int c=0; c<nbCorresp; ++c) {
			i2 = i*nbCorresp+c;
			DCMatrixCol n(transModelNrm,i);
			DCMatrixCol m(CorrNrm,i2);
			H += CorrWeights[i2]*ublas::outer_prod(n, m);
		}
	}
  // SVD with H = U D Vt

// Boost-Version does not work:
//	lapack::gesvd(H, Dv, U, Vt); // returns H = U*D*Vt
//  Ri = trans(prod(U, Vt));
// So copy everything into a CVM-Matrix and do calculations there:
  cvm::srmatrix cH(3);
  for (unsigned int row=0; row<3; ++row)
  	for (unsigned int col=0; col<3; ++col)
  		cH(row+1,col+1) = H(row,col);
  cvm::srmatrix cU(3), cVt(3);
  cvm::rvector cv;
  cv << cH.svd(cU, cVt);
  cvm::srmatrix cRi = ~(cU * cVt);
  //srmatrix estmated_Rotation = ~(cVt * cU);
  for (unsigned int row=0; row<3; ++row)
  	for (unsigned int col=0; col<3; ++col)
  		Ri(row,col) = cRi(row+1,col+1);
  cout << endl << "Rot:" << Ri << flush;


  // Double check H
//  DMatrix D(3,3,0.0);
//  D(0,0) = Dv(0);
//  D(1,1) = Dv(1);
//  D(2,2) = Dv(2);
//  DMatrix H_ = prod(U, D);
//  H_ = prod(H_, Vt); // H_ now should be = H
//  cout << endl << "H-H_: \n"<< H-H_ << endl;

//  cvm::srmatrix cD(3);
//  cD.diag(0) = cv;
////  cvm::rvector vec(3);
////  vec(1) = 1; vec(2) = 1; vec(3) = (cU * cVt).det();
////  cvm::rmatrix diag(3,3);
////  diag.diag(0) = vec;
//  cvm::srmatrix cH_ = (cU * cD) * cVt;
//	cout << "H-H_: \n"<< cH-cH_ << endl;



	/////////////////////////////////////////////////////////////
  //                  compute translation                    //
  /////////////////////////////////////////////////////////////
  if (false) {
		// The following minimizes point-plane distances:
		for (unsigned int i=0; i<nbModelPts; ++i) {
			for (unsigned int c=0; c<nbCorresp; ++c)
				{
				i2 = i*nbCorresp+c;
				DCMatrixCol p(transModelPts,i);
	//			DCMatrixCol n(transModelNrm,i);
				DCMatrixCol q(CorrPts,i2);
				DCMatrixCol m(CorrNrm,i2);
				DCMatrixCol(At,i2) = m;
				b(i2) = - ublas::inner_prod(ublas::prod(Ri,p)-q,m);
			}
		}
		A = ublas::trans(At);
		if (useWeights) {
			assert(CorrWeights.size() == At.size2());
			for (unsigned int i=0; i<At.size2(); ++i) {
				for (unsigned int j=0; j<At.size1(); ++j) {
					At(j,i) = At(j,i) * CorrWeights(i);
				}
			}
		}
		AtA = ublas::prod(At,A);
		Atb = ublas::prod(At,b);
		DCMatrix AtAinv = invSym(AtA, regCoeff);
		cout << endl << "AtA*AtAinv: \n"<< ublas::prod(AtA,AtAinv) << endl; // Double check AtAinv
		ti = ublas::prod(AtAinv,Atb);
  } else {
		// The following minimizes point-point distances
		DVector modelCenter(3);
		DVector worldCenter(3);
		double weight = 0.0;
		double cumWeight_i = 0.0;
		double cumWeight_ic = 0.0;
		for (unsigned int i=0; i<nbModelPts; ++i) {
			for (unsigned int c=0; c<nbCorresp; ++c)
				{
				i2 = i*nbCorresp+c;
				weight = CorrWeights(i2);
				worldCenter += weight * DCMatrixCol(CorrPts,i2);
				cumWeight_ic += weight;
			}
			modelCenter += weight * ublas::prod(Ri,DCMatrixCol(transModelPts,i));
			cumWeight_i += weight;
		}
		modelCenter /= cumWeight_i;
		worldCenter /= cumWeight_ic;
		ti = worldCenter - modelCenter;
  }
	cout << endl << "trans:" << ti << flush;

	// apply current transformation to total transformation:
	R = ublas::prod(Ri,R);
	t = ublas::prod(Ri,t) + ti;

	return true;
}

