
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

#include "HomogeneousTransformationMatrix.hpp"


template <typename DVectorInputIterator, typename SearchResultsT>
ICP::ICPChen<DVectorInputIterator, SearchResultsT>::
ICPChen(DVectorInputIterator ptBegin_, DVectorInputIterator ptEnd_, DVectorInputIterator nBegin_, DVectorInputIterator nEnd_, corrFuncT corrFunc, unsigned int nCorr, double regCoeff_)
	:ICPBase()
	,searchCorrespondence(corrFunc)
	,nbCorresp(nCorr)
	,regCoeff(regCoeff_)
	,useWeights(true)
	,COVAR_MAX(sqrt(DBL_MAX))
	,ptBegin(ptBegin_)
	,ptEnd(ptEnd_)
	,nBegin(nBegin_)
	,nEnd(nEnd_)
{
  init();
//  std::cout << "Covar-Max = " << COVAR_MAX << std::endl;
}

template <typename DVectorInputIterator, typename SearchResultsT>
ICP::ICPChen<DVectorInputIterator, SearchResultsT>::
~ICPChen()
{
}

template <typename DVectorInputIterator, typename SearchResultsT>
void ICP::ICPChen<DVectorInputIterator, SearchResultsT>::init()
{
	nbModelPts = 0;
	DVectorInputIterator pTmp=ptBegin;
	for (DVectorInputIterator nTmp=nBegin; nTmp!=nEnd; ++nTmp,++pTmp)
		++nbModelPts;
	if (pTmp != ptEnd)
		throw std::out_of_range("ICP-Constructor: number of points != number of normals"); // loop stopped on end of Normal-Iterator, so Point-Iterator should be "End" as well
//	cout << endl << "MP:" << modelPts << flush;
//	cout << endl << "MN:" << modelNrm << flush;
	center = *ptBegin;
//	center = DVector(3,0.0);
	transModelPts.resize(3,nbModelPts,false);
	transModelNrm.resize(3,nbModelPts,false);
	A.resize(nbModelPts*nbCorresp,6,false);
	At.resize(6,nbModelPts*nbCorresp,false);
	b.resize(nbModelPts*nbCorresp);
	varN.resize(nbModelPts*nbCorresp);
	AtA.resize(6,6,false);
	Atb.resize(6);
  if (useWeights)
  	weights.resize(nbModelPts*nbCorresp);
  ICPBase::reset();
}


template <typename DVectorInputIterator, typename SearchResultsT>
void ICP::ICPChen<DVectorInputIterator, SearchResultsT>::reset(const DMatrix &rot, const DVector &trans)
{
	ICPBase::reset(rot, trans);
  covar = DIdMatrix(6)*COVAR_MAX;
}

template <typename DVectorInputIterator, typename SearchResultsT>
bool ICP::ICPChen<DVectorInputIterator, SearchResultsT>::iterate()
{ 
	using namespace matrixTools;
	using namespace std;
	const double NORMDEVFAC(sin(20.0/180.0*M_PI)); // stdDev of normal angle = a -> deviation of point-normal distance = sin(a)*distance

	// matrix proxies are references, so no memory is copied
	typedef ublas::matrix_column< DCMatrix > DCMatrixCol;
	typedef ublas::matrix_range< DCMatrix > DCMatrixRange;
	typedef ublas::matrix_column< DCMatrixRange > DCMatrixRangeCol;
	typedef ublas::matrix_vector_range< DCMatrix > DCMatrixVRange;

	//////////////////////////////////////////////////////////
  //////   first step: apply current transformation   //////
  //////////////////////////////////////////////////////////
	if (nbModelPts < 3) return false;

	//	cout << endl << "applyTrans..." << flush;
  DVectorInputIterator pi = ptBegin;
  DVectorInputIterator ni = nBegin;
  unsigned int i = 0;
//  cout << endl << "applying R/t = " << R << t << endl;
	while ((pi != ptEnd) && (ni != nEnd))
	{
		DCMatrixCol(transModelPts,i) = prod(R,*pi)+t;
		DCMatrixCol(transModelNrm,i) = prod(R,*ni);
    ++pi;
    ++ni;
    ++i;
	}

	//////////////////////////////////////////////////////////
  /////////   second step: search correspondences   ////////
  /////////                and prepare data matrix  ////////
  //////////////////////////////////////////////////////////
//	cout << endl << "corresp..." << flush;
	double lastWAvgErr = wAvgErr;
  avgDist=0.0;
	wAvgDist=0.0;
  avgErr=0.0;
	wAvgErr=0.0;
	assert(At.size1() == 6);
  DCMatrixRange upperAt (At, ublas::range (0, 3), ublas::range (0, At.size2()));
  DCMatrixRange lowerAt (At, ublas::range (3, 6), ublas::range (0, At.size2()));
	double cumWeight=0.0;
	double dist,ppldist;
	double err;
	double weight;
	DVector pNN(3), nNN(3), cp(3); // temporary: next neighbor (correspondence), crossproduct
//	cout << endl << "looping over " << nbModelPts << " model points" << flush;
	for (unsigned int i=0; i<nbModelPts; ++i) {
//		cout << ".." << i << flush;
//		DCMatrixCol pt(transModelPts,i);
		DVector pt = DCMatrixCol(transModelPts,i);
		DCMatrixCol nr(transModelNrm,i);
		SearchResultsT sr = searchCorrespondence(pt, nr, 1.0, nbCorresp); // will not throw
//		cout << ".." << pt << flush;
		pt -= center;
//		cout << endl << "neighbor of " << pt << nr << " are " << flush;
		for (unsigned int c=0; c<nbCorresp; ++c) {
			unsigned int i2 = i*nbCorresp + c;
			try {
				pNN = sr.getCorrespPoint()-center;
				nNN = sr.getCorrespNormal();
//				cout << "->" << sr.getCorrespPoint() << flush;
				dist = sr.getDist();
				ppldist = ublas::inner_prod(pt-pNN,nNN);
				err = pow(ppldist,2); // squared point-to-plane error
				if (useWeights) {
					weight = sr.getWeight();
					if (err > lastWAvgErr) // this decreases influence of outliers
						weight *= lastWAvgErr/err;
				} else {
					weight = 1.0;
				}
//				cout << ";" << weight << flush;
				sr.next();
			} catch (exception &e) {
				//cerr << endl << e.what();
				pNN = pt;
				nNN = nr;
				dist = 0.0;
				ppldist = 0.0;
				err = 0.0;
				weight = 0.0;
			}
//			cout << endl << "pt" << pt;
//			cout << endl << "nr" << nr;
//			cout << endl << "pNN" << pNN;
//			cout << endl << "nNN" << nNN;
			avgDist += dist;
			avgErr += err;
			wAvgDist += weight*dist;
			wAvgErr += weight*err;
			cumWeight += weight;
//			cout << " => cp:" << pt << nNN << cp << flush;
			matrixTools::cross_product(pt, nNN, cp);
			DCMatrixRangeCol(upperAt,i2) = cp;
			DCMatrixRangeCol(lowerAt,i2) = nNN;
			b(i2) = - ppldist;
			varN(i2) = NORMDEVFAC*NORMDEVFAC*dist*dist;
			weights(i2) = weight;
		} // end: for each nearest neighbor to the current model point
	} // end: for each sampled model point
	avgDist /= (double)(nbModelPts*nbCorresp);
	avgErr /= (double)(nbModelPts*nbCorresp);
	wAvgDist /= cumWeight;
	wAvgErr /= cumWeight;
//	cout << endl << "max Weight = " << nbModelPts*nbCorresp << flush;
//	cout << endl << "cum Weight = " << cumWeight << flush;
//	cout << endl << "average Dist = " << avgDist << flush;
//	cout << endl << "weighted Dist = " << wAvgDist << flush;
	A = ublas::trans(At);
//	cout << endl << "A:" << A << flush;
//	cout << endl << "b:" << b << flush;
	if (useWeights) {
//		cout << endl << "weights = " << weights << flush;
//		A = ublas::prod(weights,A);
		// using a boost-band-matrix and prod() is really slow! Thus, multiply weight vector manually:
		assert(weights.size() == At.size2());
		for (unsigned int i=0; i<weights.size(); ++i) {
			for (unsigned int j=0; j<At.size1(); ++j) {
				At(j,i) = At(j,i) * weights(i);
			}
		}
	}
	AtA = 2*ublas::prod(At,A);
	Atb = 2*ublas::prod(At,b);
//	cout << endl << "AtA:" << AtA << flush;
//	cout << endl << "Atb:" << Atb << flush;

	//////////////////////////////////////////////////////////
  ///////   third step: compute new transformation   ///////
  //////////////////////////////////////////////////////////
//	KogniMobil::kogmo_timestamp_t then = KogniMobil::kogmo_timestamp_now();
	DCMatrix AtAinv = invSym(AtA, regCoeff);
	DVector x = ublas::prod(AtAinv,Atb);
//	KogniMobil::kogmo_timestamp_t then = KogniMobil::kogmo_timestamp_now();
//	cout << endl << "used " << then-now << "ns to compute";
//	cout << endl << "x=" << x << flush;
	DVector residual = b-ublas::prod(A,x);
	double sigmasquare = (ublas::inner_prod(residual,residual)+ublas::inner_prod(varN,varN))/(residual.size()-4);
	//sigmasquare = min(0.1,sigmasquare);
	//cout << endl << "sigmasquare="<<sigmasquare<<", AtAinv:" << AtAinv << flush;
	covar = AtAinv;
	covar = sigmasquare * AtAinv;

	bool estimateValid = false;
	for (unsigned int i=0; i<covar.size1(); ++i) {
		if (AtA(i,i) < regCoeff)
			covar(i,i) = COVAR_MAX;
		else
			estimateValid = true;
	}
//	cout << endl << "Covar:" << covar << flush;

	if (estimateValid) {
		// limit estimated angles, as system was linearized around 0
		for (int i=0; i<3; ++i) {
			if (fabs(x[i]) > 0.1)
				cout << endl << "limiting angle " << x[i];
			x[i] = sign(x[i])*min(0.1,fabs(x[i])); // 0.1rad = 5deg
		}
		// apply current transformation to total transformation:
		HomogeneousTransformationMatrix::YawPitchRollXYZ_2_Rt(x[2],x[1],x[0],x[3],x[4],x[5],Ri,ti);
		ti += center-ublas::prod(Ri,center); // correct for center-based estimation
		R = ublas::prod(Ri,R);
		t = ublas::prod(Ri,t) + ti;
//		cout << endl << "new R/t = " << R << t << endl;
//		cout << endl << "new Covar = " << covar << endl;
	}

	return estimateValid;
}


template <typename DVectorInputIterator, typename SearchResultsT>
void ICP::ICPChen<DVectorInputIterator, SearchResultsT>::getEstimatedVariance(DMatrix &covar_) const
{
	covar_ = covar;
}

template <typename DVectorInputIterator, typename SearchResultsT>
matrixTools::DMatrix ICP::ICPChen<DVectorInputIterator, SearchResultsT>::getEstimatedVariance() const
{
	return covar;
}

/*
void ICP::ICPChen<DVectorInputIterator, SearchResultsT>::test()
{
  cout << endl << "######## testing ICPChen ########" << flush;
  cout << endl << "- testing cross product..." << flush;
  DVector first(3), second(3), third;
  first(0)=1;first(1)=0;first(2)=0;
  second(0)=0;second(1)=1;second(2)=0;
	cross_product(first, second, third);
	assert(norm_2(third) == 1.0);
	assert(third(2) == 1.0);
	cout << "success" << flush;

	DCMatrix world(3,6);
	world(0,0) = 1.;	world(1,0) = 1.;	world(2,0) = 0.;
	world(0,1) = 2.;	world(1,1) = 1.5;	world(2,1) = 0.;
	world(0,2) = 3.;	world(1,2) = 2.;	world(2,2) = 0.;
	world(0,3) = 4.;	world(1,3) = 2.5;	world(2,3) = 0.;
	world(0,4) = 5.;	world(1,4) = 3.;	world(2,4) = 0.;
	world(0,5) = 6.;	world(1,5) = 3.5;	world(2,5) = 0.;
//	DCMatrix worldT = trans(world);
	DCMatrix model(3,5);
	model(0,0) = 0.5;	model(1,0) = 3.;	model(2,0) = 0.;
	model(0,1) = 1.5;	model(1,1) = 3.;	model(2,1) = 0.;
	model(0,2) = 2.5;	model(1,2) = 3.;	model(2,2) = 0.;
	model(0,3) = 3.5;	model(1,3) = 3.;	model(2,3) = 0.;
	model(0,4) = 4.5;	model(1,4) = 3.;	model(2,4) = 0.;
//	DCMatrix modelT = trans(model);
	DCMatrix modelN(3,5);
	modelN(0,0) = -1./sqrt(2);	modelN(1,0) = 1./sqrt(2);	modelN(2,0) = 0.;
	modelN(0,1) = 0.;	modelN(1,1) = 1.;	modelN(2,1) = 0.;
	modelN(0,2) = 0.;	modelN(1,2) = 1.;	modelN(2,2) = 0.;
	modelN(0,3) = 0.;	modelN(1,3) = 1.;	modelN(2,3) = 0.;
	modelN(0,4) = 1./sqrt(2);	modelN(1,4) = 1./sqrt(2);	modelN(2,4) = 0.;
//	DCMatrix modelNT = trans(modelN);
//	cout << endl << "world: " << world << flush;
//	cout << endl << "model: " << model << flush;

	cout << endl << "- testing MatrixColIterators..." << flush;
	DCMatrixColConstIterator i1 = DCMatrixColConstIterator::begin(world);
	DCMatrixColConstIterator i2 = DCMatrixColConstIterator::end(world);
	DCMatrixColConstIterator i3 = i1;
	assert(i1 == i3);
	assert(i1 != i2);
	for (unsigned int i=0;i<world.size2();++i) i1++;
	assert(i1 == i2);
	assert(i1 != i3);
	cout << "success" << flush;

	cout << endl << "- testing ICP..." << flush;
	CorrespondenceSearch corSearch(DCMatrixColConstIterator::begin(world),DCMatrixColConstIterator::end(world));
	ICP::ICPChen<DVectorInputIterator, SearchResultsT>::CorrespondanceFunction corFunc = boost::bind(&CorrespondenceSearch::findClosestNeighbor, &corSearch, _1, _2, _3);
	ICPChen icp(corFunc, model, modelN);
	icp.iterate();
	DVector xExp(6);
	xExp(0)=-0; xExp(1)=-0; xExp(2)=0.3; xExp(3)=2.024; xExp(4)=-1.467; xExp(5)=-0;
	DMatrix R; DVector t;
	icp.getEstimate(R,t);
//	cout << endl << "R=" << R << flush;
//	cout << endl << "t=" << t << flush;
	icp.getEstimate(R);
//	cout << endl << "htm=" << R << flush;
	icp.getEstimatedVariance(R);
//	cout << endl << "Covar=" << R << flush;
	DVector x;
	icp.getEstimate(x);
//	cout << endl << "x=" << x << flush;
//	cout << endl << "xExp=" << xExp << flush;
	double stateError = ublas::norm_2(x-xExp);
//	cout << endl << "state-error=" << stateError;
	assert(stateError < 0.001);

	cout << "success" << flush;
  cout << endl << "=======> successfully finished testing" << flush;
}
  */

