
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
#include "FMUtils.hpp"


template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
ICPLinearized(DVectorInputIterator ptBegin_, DVectorInputIterator ptEnd_,
							DVectorInputIterator nBegin_, DVectorInputIterator nEnd_,
							DInputIterator wBegin_, DInputIterator wEnd_,
							DInputIterator nConfBegin_, DInputIterator nConfEnd_,
							DInputIterator nStdDevBegin_, DInputIterator nStdDevEnd_,
							DMatrixInputIterator pCovarBegin_, DMatrixInputIterator pCovarEnd_,
							DMatrixInputIterator nCovarBegin_, DMatrixInputIterator nCovarEnd_,
							corrFuncT corrFunc, EnergyFunction &energy_,
							unsigned int nCorr, double regCoeff_, EstimationMode estMode_)
	:ICPBase()
	,estMode(estMode_)
	,searchCorrespondence(corrFunc)
	,energy(&energy_)
	,nbCorresp(nCorr)
	,regCoeff(regCoeff_)
	,useWeights(true)
	,COVAR_MAX(sqrt(DBL_MAX))
	,ERROR_MIN(1/COVAR_MAX)
	,ptBegin(ptBegin_)
	,ptEnd(ptEnd_)
	,nBegin(nBegin_)
	,nEnd(nEnd_)
  ,pCovarBegin(pCovarBegin_)
  ,pCovarEnd(pCovarEnd_)
  ,nCovarBegin(nCovarBegin_)
  ,nCovarEnd(nCovarEnd_)
	,wBegin(wBegin_)
	,wEnd(wEnd_)
	,nConfBegin(nConfBegin_)
	,nConfEnd(nConfEnd_)
	,nStdDevBegin(nStdDevBegin_)
	,nStdDevEnd(nStdDevEnd_)
  ,iteration(0)
  ,HS(6,6)
  ,Hp(6,6)
	,gS(6)
	,gp(6)
	,lastIterTransOnly(false)
{
	nbModelPts = 0;
	DVectorInputIterator pTmp=ptBegin;
	for (DVectorInputIterator nTmp=nBegin; nTmp!=nEnd; ++nTmp,++pTmp)
		++nbModelPts;
	if (pTmp != ptEnd)
		throw std::out_of_range("ICP-Constructor: number of points != number of normals"); // loop stopped on end of Normal-Iterator, so Point-Iterator should be "End" as well
	transModelPts.resize(3,nbModelPts,false);
	transModelNrm.resize(3,nbModelPts,false);
  transPCovar.resize(nbModelPts, DZeroMatrix(3));
  transNCovar.resize(nbModelPts, DZeroMatrix(3));
	CorrPts.resize(3,nbModelPts*nbCorresp,false);
	CorrNrm.resize(3,nbModelPts*nbCorresp,false);
	CorrDist = DVector(nbModelPts*nbCorresp, 0.0);
	CorrErr = DVector(nbModelPts*nbCorresp, 0.0);
	CorrWeights = DVector(nbModelPts*nbCorresp, 1.0); // assume equal full weight
	CorrNConf = DVector(nbModelPts*nbCorresp, 1.0); // assume equal full confidence
	CorrNStdDevRAD = DVector(nbModelPts*nbCorresp, 0.0); // assume equal zero uncertainty
	CorrPCovar.resize(nbModelPts*nbCorresp, DZeroMatrix(3));
  CorrNCovar.resize(nbModelPts*nbCorresp, DZeroMatrix(3));
	center = *ptBegin;
	ICPBase::reset();
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
~ICPLinearized()
{
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
void ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
reset(const DMatrix &rot, const DVector &trans)
{
	ICPBase::reset(rot, trans);
  covar = DIdMatrix(6)*COVAR_MAX;
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
bool ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
establishCorrespondences()
{
  ++iteration;
	using namespace std;
	typedef ICP::EnergyFunction::Surface Surface;
	/////////////////////////////////////////////////////////////
  // apply current transformation and search correspondences //
  /////////////////////////////////////////////////////////////
  DVectorInputIterator pi = ptBegin;
  DVectorInputIterator ni = nBegin;
  DInputIterator wi = wBegin;
  DInputIterator nci = nConfBegin;
  DInputIterator nsi = nStdDevBegin;
  DMatrixInputIterator pCi = pCovarBegin;
  DMatrixInputIterator nCi = nCovarBegin;
  unsigned int i = 0;
  unsigned int i2;
  DVector pt(3), nr(3);
  DMatrix pCovar(3,3), nCovar(3,3);
	DMatrix Rt = trans(R);
	modelCenter = DZeroVector(3);
	transModelCenter = DZeroVector(3);
	correspCenter = DZeroVector(3);
//	cout << endl << "looping over " << nbModelPts << " model points" << flush;
  unsigned int nbExcept = 0;
	while ((pi != ptEnd) && (ni != nEnd)) // for each model point
	{
//		cout << ".." << i << flush;
		pt = prod(R,*pi)+t;
		nr = prod(R,*ni);
		modelCenter += *pi;
		transModelCenter += pt;
		pCovar = (pCi != pCovarEnd) ? *pCi : DZeroMatrix(3);
		pCovar = prod(R,pCovar);
		pCovar = prod(pCovar,Rt);
		nCovar = (nCi != nCovarEnd) ? *nCi : DZeroMatrix(3);
    nCovar = prod(R,nCovar);
    nCovar = prod(nCovar,Rt);
		SearchResultsT sr = searchCorrespondence(pt, nr, *nci, nbCorresp);
//		cout << ".." << pt << flush;
		DCMatrixCol(transModelPts,i) = pt;
		DCMatrixCol(transModelNrm,i) = nr;
		transPCovar[i] = pCovar;
		transNCovar[i] = nCovar;
		for (unsigned int c=0; c<nbCorresp; ++c) { // for each correspondence of the current model point
			i2 = i*nbCorresp+c;
			try {
//				cout << "->" << DCMatrixCol(CorrPts,i2) << flush;
				DCMatrixCol(CorrPts,i2) = sr.getCorrespPoint();
				DCMatrixCol(CorrNrm,i2) = sr.getCorrespNormal();
				CorrWeights(i2) = (wi != wEnd) ? *wi : 1.0;
				if (useWeights) CorrWeights(i2) *= sr.getWeight();
				CorrDist(i2) = sr.getDist();
				CorrNConf(i2) = sr.getNormConf();
				CorrNStdDevRAD(i2) = sr.getNormStdDevRAD();
				CorrPCovar[i2] = sr.getPCovar();
				CorrNCovar[i2] = sr.getNCovar();
//				cout << ";" << weight << flush;
				Surface s1(pt, nr, (nci != nConfEnd) ? *nci : 1.0, (nsi != nStdDevEnd) ? *nsi : 0.0, pCovar, nCovar);
				Surface s2(DCMatrixCol(CorrPts,i2), DCMatrixCol(CorrNrm,i2), CorrNConf(i2), CorrNStdDevRAD(i2), CorrPCovar[i2], CorrNCovar[i2]);
        CorrErr(i2) = energy->calculate(s1,s2);
				sr.next(); // returns false if no more results available
        correspCenter += DCMatrixCol(CorrPts,i2);
			} catch (const std::exception &e) {
			  //cout << e.what();
			  ++nbExcept;
				DCMatrixCol(CorrPts,i2) = pt;
				DCMatrixCol(CorrNrm,i2) = nr;
        CorrWeights(i2) = 0.0;
        CorrDist(i2) = 0.0;
				CorrNConf(i2) = 0.0;
				CorrNStdDevRAD(i2) = M_PI/4;
        CorrPCovar[i2] = DZeroMatrix(3);
        CorrNCovar[i2] = DZeroMatrix(3);
        CorrErr(i2) = 0.0;
			}
		} // end: for each nearest neighbor to the current model point
    ++pi;
    ++ni;
		if (nci != nConfEnd) ++nci;
		if (nsi != nStdDevEnd) ++nsi;
    if (pCi != pCovarEnd) ++pCi;
    if (nCi != nCovarEnd) ++nCi;
    ++i;
	} // end: for each model point
  if (i != nbModelPts)
    throw std::logic_error("ICPLinearized::establishCorrespondences: number of points has changed!");
	cout << " nbExcept: " << nbExcept << " / " << nbModelPts*nbCorresp << flush;
  modelCenter /= (double)(nbModelPts);
	transModelCenter /= (double)(nbModelPts);
	correspCenter /= (double)(nbModelPts*nbCorresp-nbExcept);

  // calculate mean of error and distance over valid correspondences
  avgErr = 0.0;
  avgDist = 0.0;
  cumWeight = 0.0;
  vector<double> sortedCorrDst;
  for (i2=0; i2<nbModelPts*nbCorresp; ++i2) {
    if (CorrWeights[i2] > 0.0) avgErr += CorrErr[i2];
    if (CorrWeights[i2] > 0.0) avgDist += CorrDist[i2];
    if (CorrWeights[i2] > 0.0) sortedCorrDst.push_back(CorrDist[i2]);
    cumWeight += CorrWeights[i2];
  }
  if (cumWeight == 0.0)
    return false;
  sort(sortedCorrDst.begin(), sortedCorrDst.end());
  avgErr /= cumWeight;
  avgDist /= cumWeight;

  // calculate std-dev of error and distance over valid correspondences
  double stdDevDst = 0.0;
  for (i2=0; i2<nbModelPts*nbCorresp; ++i2) {
    if (CorrWeights[i2] > 0.0) stdDevDst += CorrWeights[i2]*pow(CorrDist[i2]-avgDist, 2);
  }
  stdDevDst = sqrt(stdDevDst/cumWeight);

  // define top 10% as outliers and re-weight correspondences
  unsigned int outlierIdx = (double)sortedCorrDst.size()*0.9; // 10 % outlier
  if (iteration <= 2) outlierIdx = sortedCorrDst.size()-1; // disable outlier-detection on first two iterations
//  outlierIdx = sortedCorrDst.size()-1; // disable outlier-detection
  // then, re-calculate weights based on outlier-rejection and calculate weight-based statistics
  wAvgDist=0.0;
  wAvgErr=0.0;
  cumWeight=0.0;
  outlierRatio=0.0;
  for (i2=0; i2<nbModelPts*nbCorresp; ++i2) {
    //if (CorrDst[i2]-avgDst > 2*stdDevDst)
    if (CorrDist[i2] > sortedCorrDst[outlierIdx])
      CorrWeights[i2] = 0.0;
    if (CorrWeights[i2] == 0.0)
      outlierRatio += 1.0;
    wAvgDist += CorrWeights[i2]*CorrDist[i2];
    wAvgErr += CorrWeights[i2]*CorrErr[i2];
    cumWeight += CorrWeights[i2];
	}
  outlierRatio /= (double)(nbModelPts*nbCorresp);
	wAvgDist /= cumWeight;
	wAvgErr /= cumWeight;
  cout << " outlier-ratio: " << outlierRatio << " dst@0.8: " << sortedCorrDst[outlierIdx] << ", max: " << sortedCorrDst[sortedCorrDst.size()-1] << ", stdDevDst: " << stdDevDst << " " << flush;
	return (cumWeight >= 1.0);
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
bool ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
minimizeError()
{
	using namespace std;
	using namespace matrixTools;
	//////////////////////////////////////////////////////
  //////   first step: prepare data matrix        //////
  //////////////////////////////////////////////////////
  HS=ublas::zero_matrix<double>(6,6);
  gS=ublas::zero_vector<double>(6);
  ICP::EnergyFunction::Surface s,sNN;
	double weight;
//	cout << endl << "looping over " << nbModelPts << " model points" << flush;
  DInputIterator nci = nConfBegin;
  DInputIterator nsi = nStdDevBegin;
	for (unsigned int i=0; i<nbModelPts; ++i)
	{
//		cout << ".." << i << flush;
		s.p = DCMatrixCol(transModelPts,i) - center; // Make calculations center-based (better numerical stability)
		s.n = DCMatrixCol(transModelNrm,i);
		s.nConfidence = (nci != nConfEnd) ? *nci : 1.0;
		s.nStdDevRAD = (nsi != nStdDevEnd) ? *nsi : 0.0;
		s.pCovar3D = transPCovar[i];
		s.nCovar3D = transNCovar[i];
		for (unsigned int c=0; c<nbCorresp; ++c) {
//			cout << "(" << c << ")" << flush;
			unsigned int i2 = i*nbCorresp+c;
			sNN.p = DCMatrixCol(CorrPts,i2) - center; // Make calculations center-based (better numerical stability)
			sNN.n = DCMatrixCol(CorrNrm,i2);
			sNN.nConfidence = CorrNConf(i2);
			sNN.nStdDevRAD = CorrNStdDevRAD(i2);
	    sNN.pCovar3D = CorrPCovar[i2];
	    sNN.nCovar3D = CorrNCovar[i2];
			weight = CorrWeights(i2);
			energy->approximateAt(s, sNN, Hp, gp);
      HS += weight*Hp;
      gS += weight*gp;
      for (size_t c=0; c<gS.size(); ++c) {
        assert(!std::isnan(gS(c)) && "ICPLinearized::minimizeError: gS is NAN");
        for (size_t r=0; r<gS.size(); ++r) {
          assert(!std::isnan(HS(r,c)) && "ICPLinearized::minimizeError: HS is NAN");
        }
      }
		} // end: for each nearest neighbor to the current model point
		if (nci != nConfEnd) ++nci;
		if (nsi != nStdDevEnd) ++nsi;
	} // end: for each sampled model point

	//////////////////////////////////////////////////////////
  ///////   second step: compute new transformation  ///////
  //////////////////////////////////////////////////////////
  DVector x = ublas::zero_vector<double>(6);
  DCMatrix HSinv = DIdMatrix(6,6) * 1e-6;
  switch (estMode) {
    case (Mode3R3T) : {
      HSinv = invSym(HS, regCoeff);
      x = ublas::prod(HSinv,-gS);
      //  cout << endl << "H=" << aligned_write(HS) << flush;
      //  cout << endl << "HInv=" << aligned_write(HSinv) << endl;
    } break;
    case (Mode1R3T): {
      DCMatrix HSinvPartial = invSym(DSMatrixRange(HS, ublas::range(2,6), ublas::range(2,6)), regCoeff);
      DVectorRange(x,ublas::range(2,6)) = ublas::prod(HSinvPartial,-DVectorRange(gS, ublas::range(2,6)));
      DCMatrixRange(HSinv, ublas::range(2,6), ublas::range(2,6)) = HSinvPartial; // plug partial covariance back into full matrix
    } break;
  }

	//////////////////////////////////////////////////////////
  ///////   third step: compute covariance matrix    ///////
  //////////////////////////////////////////////////////////
	double sigmasquare = 2 * getError() / (cumWeight); // -3 !!!!!
//	double sigmasquare = (ublas::inner_prod(residual,residual)+ublas::inner_prod(varN,varN))/(residual.size()-4);
	covar = HSinv*sigmasquare; // this is only valid, if the current estimation did not change anything, otherwise HSinv should be recomputed with the newly calculated transformation

	bool estimateValid = false;
	for (unsigned int i=0; i<covar.size1(); ++i) {
		if (HS(i,i) < regCoeff)
			covar(i,i) = COVAR_MAX;
		else
			estimateValid = true;
	}
//	cout << endl << "Covar:" << covar << flush;

	if (estimateValid) {
		// limit estimated angles, as system was linearized around 0
		for (int i=0; i<3; ++i) {
			if (fabs(x[i]) > 0.1) {
				cout << endl << "limiting angle " << x[i];
				x[i] = sign(x[i])*min(0.1,fabs(x[i])); // 0.1rad = 5deg
				estimateValid = false;
			}
		}
	  estimateValid &= (cumWeight >= 3); //three fully valid correspondences needed at minimum
		// in case angles had to be limited disable rotation completely (if not disabled last iteration)
		if (!estimateValid && !lastIterTransOnly) {
      cout << endl << "disabling rotation ";
		  // set rotation to zero and translation to center-distance
		  x[0] = 0;
      x[1] = 0;
      x[2] = 0;
      x[3] = correspCenter[0] - transModelCenter[0];
      x[4] = correspCenter[1] - transModelCenter[1];
      x[5] = correspCenter[2] - transModelCenter[2];
      lastIterTransOnly = true;
		} else {
      lastIterTransOnly = false;
      // TODO (1): re-calculate rotation based on Pt-Pt energy
		}
		// apply current transformation to total transformation:
		HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt(x[2],x[1],x[0],x[3],x[4],x[5],Ri,ti);
		ti += center-ublas::prod(Ri,center); // correct for center-based estimation
	  assert(!std::isnan(ti(0)) && "ICPLinearized::minimizeError: ti is NAN");
    assert(!std::isnan(ti(1)) && "ICPLinearized::minimizeError: ti is NAN");
    assert(!std::isnan(ti(2)) && "ICPLinearized::minimizeError: ti is NAN");
		R = ublas::prod(Ri,R);
		t = ublas::prod(Ri,t) + ti;
    assert(!std::isnan(t(0)) && "ICPLinearized::minimizeError: t is NAN");
    assert(!std::isnan(t(1)) && "ICPLinearized::minimizeError: t is NAN");
    assert(!std::isnan(t(2)) && "ICPLinearized::minimizeError: t is NAN");
	}

	return estimateValid;
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
bool ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
iterate()
{ 
	return (establishCorrespondences() && minimizeError());
}


template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
void ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getErrorStatistics(double &avgWError, DSMatrix &measCovar, DSMatrix &hessian) const
{
  DVectorInputIterator pi = ptBegin;
  DVectorInputIterator ni = nBegin;
  DInputIterator nci = nConfBegin;
  DInputIterator nsi = nStdDevBegin;
//  DSMatrix measCovar = ublas::zero_matrix<double>(6,6);
  measCovar  = ublas::zero_matrix<double>(6,6);
  DSMatrix measCovar_i(6,6);
  hessian = ublas::zero_matrix<double>(6,6);
  DSMatrix hessian_i(6,6);
  double weight_i;
  avgWError = 0.0;
  ICP::EnergyFunction::Surface s, sNN;
  unsigned int i=0;
  while ((pi != ptEnd) && (ni != nEnd)) // for each model point
  {
//    s.p = DCMatrixCol(transModelPts,i) - center; // Make calculations center-based (better numerical stability)
//    s.n = DCMatrixCol(transModelNrm,i);
    s.p = DCMatrixColConst(transModelPts,i);
    s.n = DCMatrixColConst(transModelNrm,i);
    s.nConfidence = (nci != nConfEnd) ? *nci : 1.0;
    s.nStdDevRAD = (nsi != nStdDevEnd) ? *nsi : 0.0;
    s.pCovar3D = transPCovar[i];
    s.nCovar3D = transNCovar[i];
    for (unsigned int c=0; c<nbCorresp; ++c) {
      unsigned int i2 = i*nbCorresp+c;
      sNN.p = DCMatrixColConst(CorrPts,i2);
      sNN.n = DCMatrixColConst(CorrNrm,i2);
      sNN.nConfidence = CorrNConf(i2);
      sNN.nStdDevRAD = CorrNStdDevRAD(i2);
      sNN.pCovar3D = CorrPCovar[i2];
      sNN.nCovar3D = CorrNCovar[i2];
      weight_i = CorrWeights(i2);
      energy->approximateCovar(s, sNN, hessian_i, measCovar_i);
      avgWError += weight_i * energy->calculate(s, sNN);
      hessian += weight_i * hessian_i;
      measCovar += weight_i * measCovar_i;
    } // end: for each nearest neighbor to the current model point
    ++ni;
    ++pi;
    if (nci != nConfEnd) ++nci;
    if (nsi != nStdDevEnd) ++nsi;
    ++i;
  } // end: for each sampled model point
  avgWError /= cumWeight;
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::DSMatrix ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getEstimatedVariance() const
{
  //return covar; // buffered covar could be used too, BUT: sometimes a "reset" also resets the covar - e.g. during sampling
	double avgError;
  DSMatrix measCovar = ublas::zero_matrix<double>(6,6);
  DSMatrix hessian = ublas::zero_matrix<double>(6,6);
  getErrorStatistics(avgError, measCovar, hessian);
//  std::cout << std::endl << "hessian:" << std::endl << aligned_write(hessian) << std::endl;

	// covar_x = 2E/n-3 * ExxI
	DSMatrix hessianI = invSym(hessian, 0.00001);
	double sigmasquare = 2 * avgError;
	sigmasquare = std::max(0.001,sigmasquare); // HACK
	DSMatrix covar_ = hessianI * sigmasquare;
//	std::cout << std::endl << "old: " << covar << std::endl << "new: " << covar_ << std::endl;
//  std::cout << std::endl << "-->covar(sigma):" << sigmasquare << std::endl << aligned_write(covar_) << std::endl;
	return covar_;
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::DSMatrix ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getAnalyticalVariance() const
{
  double avgError;
  DSMatrix measCovar = ublas::zero_matrix<double>(6,6);
  DSMatrix hessian = ublas::zero_matrix<double>(6,6);
  getErrorStatistics(avgError, measCovar, hessian);
//  std::cout << std::endl << "hessian:" << std::endl << aligned_write(hessian) << std::endl;
//  std::cout << std::endl << "measCovar:" << std::endl << aligned_write(measCovar) << std::endl;

	// covar_x = ExxI * (Ezx * covar_z * EzxT) * ExxI, Exx=hessian, (Ezx * covar_z * EzxT)=transMeasCovar
	DSMatrix hessianI = invSym(hessian, 0.00001);
//  std::cout << std::endl << "hessianInv:" << std::endl << aligned_write(hessianI) << std::endl;
	DMatrix covar_ = ublas::prod(hessianI,measCovar); // one product only is not symmetric any more
	covar_ = ublas::prod(covar_,hessianI); // this is symmetric again
//  std::cout << std::endl << "-->covar:" << std::endl << aligned_write(covar_) << std::endl;
	return ublas::symmetric_adaptor<DMatrix, ublas::upper>(covar_);
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::DSMatrix ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getSampledVariance(bool correspondenceSearch)
{
	using namespace std;

	DVector x0; getEstimate(x0);
	double b0 = getError();
	cout << endl << "sampling (corresp="<<correspondenceSearch<<") based on x=" << x0 << " with error " << b0 << endl;

	const double tVar = 1.0; // sample translation parameters up to this value
	const double rVar = M_PI/8.0; // sample rotation parameters up to this value

	/*
	// sampling based on random coordinates
	unsigned int nb = 100; // number of samples to draw
	EmpCovarEstimator covarEstimator(*this, nb, correspondenceSearch);
	cout << endl << "sampling (corresp="<<correspondenceSearch<<") based on x=" << x0 << " with error " << b0 << endl;
	DVector xd(6); // sampled transformation (difference) = (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
	for (unsigned int i=0; i<nb; ++i) {
		xd(0) = randomInRange(-rVar,rVar);
		xd(1) = randomInRange(-rVar,rVar);
		xd(2) = randomInRange(-rVar,rVar);
		xd(3) = randomInRange(-tVar,tVar);
		xd(4) = randomInRange(-tVar,tVar);
		xd(5) = randomInRange(-tVar,tVar);
		covarEstimator.generateSample(xd);
	}
	*/

  /*
  // sampling based on unscented transform
  unsigned int nb = 12; // number of samples to draw (+-stdDev per dimension)
  EmpCovarEstimator covarEstimator(*this, nb, correspondenceSearch);
  cout << endl << "sampling (corresp="<<correspondenceSearch<<") based on x=" << x0 << " with error " << b0 << endl;
  */

	// sampling based on a regular lattice:
	unsigned int nb = pow(2,6)-1; // number of samples to draw (don't sample 0,0,0,0,0,0)
	EmpCovarEstimator covarEstimator(*this, nb, correspondenceSearch);
	DVector xd(6); // sampled transformation (difference) = (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
	for (int irx = 0; irx <= 1; ++irx)
		for (int iry = 0; iry <= 1; ++iry)
			for (int irz = 0; irz <= 1; ++irz)
				for (int itx = 0; itx <= 1; ++itx)
					for (int ity = 0; ity <= 1; ++ity)
						for (int itz = 0; itz <= 1; ++itz) {
							if (irx==0 && iry==0 && irz==0 && itx==0 && ity==0 && itz==0)
								continue;
							xd(0) = irx*rVar;
							xd(1) = iry*rVar;
							xd(2) = irz*rVar;
							xd(3) = itx*tVar;
							xd(4) = ity*tVar;
							xd(5) = itz*tVar;
							covarEstimator.generateSample(xd);
						}

	return covarEstimator.getCovar(false);//true);
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
void ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
testSampledVariance(bool correspondenceSearch)
{
	using namespace std;

	DVector x0; getEstimate(x0);
	double b0 = getError();
	cout << endl << "testing covariance sampling (corresp="<<correspondenceSearch<<") based on x=" << x0 << " with error " << b0 << endl;

	const double tVar = 1.0; // sample translation parameters up to this value
	const double rVar = M_PI/8.0; // sample rotation parameters up to this value

	// sampling based on a regular lattice:
	unsigned int nb = pow(2,6)-1; // number of samples to draw (don't sample 0,0,0,0,0,0)
	EmpCovarEstimator covarEstimator(*this, nb, correspondenceSearch);
	DVector xd(6); // sampled transformation (difference) = (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
	for (int irx = 0; irx <= 1; ++irx)
		for (int iry = 0; iry <= 1; ++iry)
			for (int irz = 0; irz <= 1; ++irz)
				for (int itx = 0; itx <= 1; ++itx)
					for (int ity = 0; ity <= 1; ++ity)
						for (int itz = 0; itz <= 1; ++itz) {
							if (irx==0 && iry==0 && irz==0 && itx==0 && ity==0 && itz==0)
								continue;
							xd(0) = irx*rVar;
							xd(1) = iry*rVar;
							xd(2) = irz*rVar;
							xd(3) = itx*tVar;
							xd(4) = ity*tVar;
							xd(5) = itz*tVar;
							covarEstimator.generateSample(xd);
						}
	DMatrix covar_;
	covarEstimator.getCovar(covar_, true);

	// sampling based on random coordinates
	DMatrix covarold(1,1);
	for (unsigned int nbi = 30; nbi < 2*nb; nbi += 10) { // number of samples to draw must be > (6+5+..)=21
		EmpCovarEstimator covarEstimator(*this, nbi, correspondenceSearch);
		cout << nbi << "..." << flush;
		DVector xd(6); // sampled transformation (difference) = (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
		for (unsigned int i=0; i<nbi; ++i) {
			xd(0) = randomInRange(-2*rVar,2*rVar);
			xd(1) = randomInRange(-2*rVar,2*rVar);
			xd(2) = randomInRange(-2*rVar,2*rVar);
			xd(3) = randomInRange(-2*tVar,2*tVar);
			xd(4) = randomInRange(-2*tVar,2*tVar);
			xd(5) = randomInRange(-2*tVar,2*tVar);
			covarEstimator.generateSample(xd);
		}
		DMatrix covari;
		covarEstimator.getCovar(covari, true);
		if (covarold.size1() == covari.size1()) { // not the first iteration
			double diff = 0.0;
			for (unsigned int row = 0; row < covari.size1(); ++row)
				for (unsigned int col = 0; col < covari.size2(); ++col) {
					diff += fabs(covari(row,col)-covarold(row,col));
				}
			cout << "e=" << diff << "..." << flush;
		}
		covarold = covari;
	}

	double diff = 0.0;
	for (unsigned int row = 0; row < covar_.size1(); ++row)
		for (unsigned int col = 0; col < covar_.size2(); ++col) {
			diff += fabs(covar_(row,col)-covarold(row,col));
		}
	cout << "e=" << diff << "..." << flush;

}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
double ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getError() const
{
	return getError(R, t);
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
double ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getError(const DMatrix &Rot, const DVector &trans) const
{
  DVectorInputIterator pi = ptBegin;
  DVectorInputIterator ni = nBegin;
  DInputIterator nci = nConfBegin;
  DInputIterator nsi = nStdDevBegin;
  DMatrixInputIterator pCi = pCovarBegin;
  DMatrixInputIterator nCi = nCovarBegin;
  double error = 0.0;
	double weight;
	ICP::EnergyFunction::Surface s, sNN;
	unsigned int i=0;
	DMatrix Rt = ublas::trans(Rot);
	while ((pi != ptEnd) && (ni != nEnd)) // for each model point
	{
		s.p = prod(Rot,*pi)+trans;
		s.n = prod(Rot,*ni);
		s.pCovar3D = (pCi != pCovarEnd) ? *pCi : DZeroMatrix(3);
		s.pCovar3D = prod(Rot,s.pCovar3D);
		s.pCovar3D = prod(s.pCovar3D,Rt);
		s.nCovar3D = (nCi != nCovarEnd) ? *nCi : DZeroMatrix(3);
		s.nCovar3D = prod(Rot,s.nCovar3D);
		s.nCovar3D = prod(s.nCovar3D,Rt);
		s.nConfidence = (nci != nConfEnd) ? *nci : 1.0;
		s.nStdDevRAD = (nsi != nStdDevEnd) ? *nsi : 0.0;
		for (unsigned int c=0; c<nbCorresp; ++c) {
			unsigned int i2 = i*nbCorresp+c;
			sNN.p = DCMatrixColConst(CorrPts,i2);
			sNN.n = DCMatrixColConst(CorrNrm,i2);
			sNN.nConfidence = CorrNConf(i2);
			sNN.nStdDevRAD = CorrNStdDevRAD(i2);
      sNN.pCovar3D = CorrPCovar[i2];
      sNN.nCovar3D = CorrNCovar[i2];
			weight = CorrWeights(i2);
			error += weight * energy->calculate(s, sNN);
		} // end: for each nearest neighbor to the current model point
		++ni;
		++pi;
		if (nci != nConfEnd) ++nci;
		if (nsi != nStdDevEnd) ++nsi;
		++i;
	} // end: for each sampled model point
	return error;
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
double ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getAvgWeight() const
{
  return cumWeight / (double)(nbModelPts*nbCorresp);
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
double ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
getSignificantWeightRatio() const
{
  double cnt = 0.0;
  for (unsigned int i2=0; i2<nbModelPts*nbCorresp; ++i2) {
    if (CorrWeights[i2] > 0.5)
      cnt += 1.0;
  }
  return cnt / (double)(nbModelPts*nbCorresp);
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
EmpCovarEstimator::EmpCovarEstimator(ICPLinearized &icp_, unsigned int nbSamples, bool resetCorresp)
: icp(icp_)
	,nb(nbSamples)
	,correspondenceSearch(resetCorresp)
	,xt(6)
	,Rot(3,3)
	,trans(3)
  ,A(nb,6)
	,idx(0)
{
	icp.getEstimate(x0);
	b0 = fabs(icp.getAvgWeightError());
	expAvgError = 5*b0; //0.1;
	errMin = b0/500;
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
EmpCovarEstimator::~EmpCovarEstimator()
{
	if (correspondenceSearch) {
		HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt(x0[2],x0[1],x0[0],x0[3],x0[4],x0[5],Rot,trans);
		icp.reset(Rot, trans);
		icp.establishCorrespondences();
	}
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
void ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
EmpCovarEstimator::generateSample(const DVector &xd) {
	if (idx >= nb)
		throw std::range_error("EmpCovarEstimator::generateSample: a sample is generated out of range");

	/*
	 * The following approach directly estimates the covariance matrix
	 * by generating samples that produce an expected error value
	 * Therefore a real sample with different error is re-mapped via a quadratic
	 * error function to a virtual sample with expected error
	 * From all these samples the empirical covariance can be directly calculated
	 */
	// 1) determine transformation and get error:
	xt = x0 + xd; // (rx-roll,ry-pitch,rz-yaw,tx,ty,tz)
	HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt(xt[2],xt[1],xt[0],xt[3],xt[4],xt[5],Rot,trans);
	if (correspondenceSearch) {
		icp.reset(Rot, trans);
		icp.establishCorrespondences();
	}
	b = fabs(icp.getAvgWeightError(Rot,trans)); // average error per correspondence
	bs.push_back(b);
	b = std::max(errMin,b-b0); // make error relative, avoid negative error values
	//nx = ublas::norm_2(xd); // distance from center of quadratic error function
	// 2) estimate parameter:
	// assume pure quadratic error function: b = f(nx) = c*nx^2
	// now given: b and nx --> c = b/nx^2
	// want to know: nx' where b = expAvgError --> nx' = sqrt(expAvgError/c)
	// --> scale = nx'/nx = sqrt(expAvgError/c)/nx = sqrt(expAvgError/c/nx^2)
	//           = sqrt(expAvgError/(b/nx^2)/nx^2) = sqrt(expAvgError/b)
	scale = sqrt(expAvgError/b); // this scales the point onto the std-dev ellipsoid
	//scale = randomGauss(0.0, scale*scale);
	DMatrixRow(A,idx) = scale * xd;
	++idx;

	/*
	 * The following approach first estimates the error function by a quadratic function
	 * Then this quadratic function is inverted to get the covariance
	 * Drawbacks: LS might result in a partly negative error function --> covariance might not be positive semi-definite!
	 * Second, it is not clear which factor should be used on the inverted function to get an appropriate covariance
	 */
	// xMx = rx^2*m[1,1] + 2*rx*ry*m[2,1] + 2*rx*rz*m[3,1] + 2*rx*tx*m[4,1] + 2*rx*ty*m[5,1] + 2*rx*tz*m[6,1] + ry^2*m[2,2] + 2*ry*rz*m[3,2] + 2*ry*tx*m[4,2] + 2*ry*ty*m[5,2] + 2*ry*tz*m[6,2] + rz^2*m[3,3] + 2*rz*tx*m[4,3] + 2*rz*ty*m[5,3] + 2*rz*tz*m[6,3] + tx^2*m[4,4] + 2*tx*ty*m[5,4] + 2*tx*tz*m[6,4] + ty^2*m[5,5] + 2*ty*tz*m[6,5] + tz^2*m[6,6]
//	A(idx,0) = x(0)*x(0);
//	A(idx,1) = 2*x(0)*x(1);
//	A(idx,2) = 2*x(0)*x(2);
//	A(idx,3) = 2*x(0)*x(3);
//	A(idx,4) = 2*x(0)*x(4);
//	A(idx,5) = 2*x(0)*x(5);
//	A(idx,6) = x(1)*x(1);
//	A(idx,7) = 2*x(1)*x(2);
//	A(idx,8) = 2*x(1)*x(3);
//	A(idx,9) = 2*x(1)*x(4);
//	A(idx,10) = 2*x(1)*x(5);
//	A(idx,11) = x(2)*x(2);
//	A(idx,12) = 2*x(2)*x(3);
//	A(idx,13) = 2*x(2)*x(4);
//	A(idx,14) = 2*x(2)*x(5);
//	A(idx,15) = x(3)*x(3);
//	A(idx,16) = 2*x(3)*x(4);
//	A(idx,17) = 2*x(3)*x(5);
//	A(idx,18) = x(4)*x(4);
//	A(idx,19) = 2*x(4)*x(5);
//	A(idx,20) = x(5)*x(5);
//	b(idx) = b;
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
double ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
EmpCovarEstimator::getLastError()
{
	if (bs.empty()) return 0.0;
	return bs.back();
}

template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
matrixTools::DVector ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
EmpCovarEstimator::getLastSample()
{
	if (idx == 0) return DVector(6,0.0);
	return DMatrixRow(A,idx-1);
}
template <typename DMatrixInputIterator, typename DVectorInputIterator, typename DInputIterator, typename SearchResultsT>
ICP::DSMatrix ICP::ICPLinearized<DMatrixInputIterator, DVectorInputIterator, DInputIterator, SearchResultsT>::
EmpCovarEstimator::getCovar(bool normalizeVolume)
{
	using namespace std;
	if (idx != nb)
		cout << endl << "idx now is " << idx << " but matrix-size is " << nb;

	/*
	 * The following approach directly estimates the covariance matrix
	 * by generating samples that produce an expected error value
	 * Therefore a real sample with different error is re-mapped via a quadratic
	 * error function to a virtual sample with expected error
	 * From all these samples the empirical covariance can be directly calculated
	 */
//	std::sort(bs.begin(), bs.end());
//	cout << " errors ranging from " << bs.front() << " tp " << bs.back() << endl;
//	BOOST_FOREACH(double bi, bs) cout << bi << "  "; cout << endl;

	// calculate covariance:
	DMatrix At = ublas::trans(A);
	DSMatrix covar_;
	covar_ = ublas::prod(At,A);
	covar_ /= (double)(nb-1);

	// normalize volume:
	if (normalizeVolume) {
	  cout << "determining det of " << covar_ << flush;
		double det = determinant(covar_);
		// det(c*A) = c^d * det(A) , with d=dimensionality
		// volume of ellipsoid ~ sqrt(det(A)), in 1d equal to variance
		if (det != 0.0)
			covar_ *= 0.5/pow(det,1/6.0); // 6-dimensions
		else
			cout << "WARNING: determinant is < 0.0000001" << flush;
		//cout << "determinant is now " << determinant(covar_) << flush;
	}
	return covar_;

	/*
	 * The following approach first estimates the error function by a quadratic function
	 * Then this quadratic function is inverted to get the covariance
	 * Drawbacks: LS might result in a partly negative error function --> covariance might not be positive semi-definite!
	 * Second, it is not clear which factor should be used on the inverted function to get an appropriate covariance
	 */
	//	cout << "A=" << endl << aligned_write(A) << endl;
	//	cout << "b=" << endl << b << endl;
	//	DMatrix At = ublas::trans(A);
	//	DMatrix AtA = ublas::prod(At,A);
	//	DVector Atb = ublas::prod(At,b);
	//	cout << "AtA:" << endl << aligned_write(AtA) << endl;
	//	DMatrix AtAinv = matrixTools::invSym(AtA); // no regularization required, as enough points were sampled
	//	DVector cv = ublas::prod(AtAinv,Atb); // least squares estimate of the quadratic error function in vector-form
	//	DVector residual = b-ublas::prod(A,cv);
	//	cout << "Avg.residual = " << ublas::norm_1(residual)/(double)nb << endl;
	//	DSMatrix M(6,6); // quadratic error function in matrix from <- retrieve from vector form
	//	M(0,0) = cv(0); M(0,1) = cv(1);  M(0,2) = cv(2);  M(0,3) = cv(3);  M(0,4) = cv(4);  M(0,5) = cv(5);
	//	M(1,0) = cv(1); M(1,1) = cv(6);  M(1,2) = cv(7);  M(1,3) = cv(8);  M(1,4) = cv(9);  M(1,5) = cv(10);
	//	M(2,0) = cv(2); M(2,1) = cv(7);  M(2,2) = cv(11); M(2,3) = cv(12); M(2,4) = cv(13); M(2,5) = cv(14);
	//	M(3,0) = cv(3); M(3,1) = cv(8);  M(3,2) = cv(12); M(3,3) = cv(15); M(3,4) = cv(16); M(3,5) = cv(17);
	//	M(4,0) = cv(4); M(4,1) = cv(9);  M(4,2) = cv(13); M(4,3) = cv(16); M(4,4) = cv(18); M(4,5) = cv(19);
	//	M(5,0) = cv(5); M(5,1) = cv(10); M(5,2) = cv(14); M(5,3) = cv(17); M(5,4) = cv(19); M(5,5) = cv(20);
	////	for (unsigned int d=0; d<M.size1(); ++d) M(d,d) = max(0.0,M(d,d)); // diagonal must be positive
	//	cout << "sampled error-function:" << endl << aligned_write(M) << endl;
	//	covar_ = matrixTools::invSym(M); // turn quadratic-function matrix into covariance
	//	cout << "ID6=" << endl << aligned_write(ublas::prod(M,covar_)) << endl;
	//	cout << "sampled Covar:" << endl << aligned_write(covar_) << endl;
}

