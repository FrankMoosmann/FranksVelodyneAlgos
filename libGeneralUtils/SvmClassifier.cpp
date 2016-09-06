/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    \file SvmClassifier.cpp

    $Author: Frank Moosmann $ ( <frank.moosmann@inrialpes.fr>)
    $Date: 2006/01/14 $
*///////////////////////////////////////////////////////////////////////////
 

#include "SvmClassifier.hpp"

#include <iostream>
#include <istream>
#include <fstream>
#include <cmath>
#include <list>
#include <string>
#include <boost/lexical_cast.hpp>

using namespace std;
 

// non-member-methods for parsing a string into enum-types with boost::program_options
void validate(boost::any& v, const vector<string>& values, SvmType*, int)
{
  namespace po = boost::program_options;
  po::validators::check_first_occurrence(v); // Make sure no previous assignment to 'v' was made.
  const string& s = po::validators::get_single_string(values); // Extract the first string from 'values'. If there is more than one string, it's an error, and exception will be thrown.
  if (s == "C_SVC") {v = boost::any(SvmType(C_SVC)); return; }
  if (s == "NU_SVC") {v = boost::any(SvmType(NU_SVC)); return; }
  if (s == "ONE_CLASS") {v = boost::any(SvmType(ONE_CLASS)); return; }
  if (s == "EPSILON_SVR") {v = boost::any(SvmType(EPSILON_SVR)); return; }
  if (s == "NU_SVR") {v = boost::any(SvmType(NU_SVR)); return; }
  throw po::validation_error(po::validation_error::invalid_option_value, "invalid value");
}
void validate(boost::any& v, const vector<string>& values, SvmKernelType*, int)
{
  namespace po = boost::program_options;
  po::validators::check_first_occurrence(v); // Make sure no previous assignment to 'v' was made.
  const string& s = po::validators::get_single_string(values); // Extract the first string from 'values'. If there is more than one string, it's an error, and exception will be thrown.
  if (s == "LINEAR") {v = boost::any(LINEAR); return; }
  if (s == "POLY") {v = boost::any(POLY); return; }
  if (s == "RBF") {v = boost::any(RBF); return; }
  if (s == "SIGMOID") {v = boost::any(SIGMOID); return; }
  if (s == "PRECOMPUTED") {v = boost::any(PRECOMPUTED); return; }
  throw po::validation_error(po::validation_error::invalid_option_value, "invalid value");
}
//void validate(boost::any& v, const vector<string>& values, SvmClassifier::SvmParams::WeightMap*, int)
//{
//  cout << "values: " << values.size() << endl;
//  namespace po = boost::program_options;
//  po::validators::check_first_occurrence(v); // Make sure no previous assignment to 'v' was made.
//  if (values.size()%2 != 0)
//    throw po::validation_error(po::validation_error::invalid_option_value, "invalid number of argument values");
//  map<ClassLabel, float> retMap;
//  for (size_t i=0; i<values.size()/2; ++i) {
//    ClassLabel c= boost::lexical_cast<ClassLabel>(values[i*2]);
//    float f = boost::lexical_cast<float>(values[i*2+1]);
//    cout << " - " << c << " : " << f << endl;
//    retMap[c] = f;
//  }
//  cout << "mapsize: " << retMap.size() << endl;
//  v = retMap;
//}

/***************************************************************************
                         Constructors / Destructors
 ***************************************************************************/
SvmClassifier::SvmParams::SvmParams()
  : parameterDescription("SVM-Parameter")
{
  svm_type = C_SVC;      // available: C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR
  kernel_type = LINEAR; // available: LINEAR [u'*v], POLY [(gamma*u'*v + coef0)^degree], RBF[exp(-gamma*|u-v|^2)], SIGMOID [tanh(gamma*u'*v + coef0)]
  degree = 5;           // degree in kernel function
  gamma = 0.01;         // gamma in kernel function (default 1/k)
  coef0 = 0;            // coef0 in kernel function (default 0)
  nu = 0.5;             // parameter nu of nu-SVC, one-class SVM, and nu-SVR
  cache_size = 500;     // in MB
  C = 100;              // cost-parameter C for C_SVC, EPSILON_SVR and NU_SVR
  eps = 1e-3;           // stopping criteria, usually use 0.00001 in nu-SVC, 0.001 in others
  p = 0.1;              // for EPSILON_SVR
  shrinking = 1;        // use the shrinking heuristics {0,1}
  probability = 0;      // do probability estimates {0,1}
  nr_weight = 0;        // for C_SVC
  weight_label = NULL;  // for C_SVC
  weight = NULL;        // for C_SVC
  weightBalanceTrainingSet = false; // if true the above 3 weight variables will be set according to training set class label distribution
  namespace po = boost::program_options;
  parameterDescription.add_options()
    ("svmType", po::value(&svm_type)->default_value(C_SVC), "available: C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR")
    ("svmKernel", po::value(&kernel_type)->default_value(LINEAR), "available: LINEAR [u'*v], POLY [(gamma*u'*v + coef0)^degree], RBF[exp(-gamma*|u-v|^2)], SIGMOID [tanh(gamma*u'*v + coef0)]")
    ("svmDegree", po::value(&degree)->default_value(degree), "degree in polynomial kernel function")
    ("svmGamma", po::value(&gamma)->default_value(gamma), "gamma in poly/rbf/sigmoid kernel function (default 1/k)")
    ("svmCoef0", po::value(&coef0)->default_value(coef0), "coef0 in poly/sigmoid kernel function (default 0)")
    ("svmNu", po::value(&nu)->default_value(nu), "parameter nu of nu-SVC, one-class SVM, and nu-SVR")
    ("svmCacheSize", po::value(&cache_size)->default_value(cache_size), "in MB")
    ("svmC", po::value(&C)->default_value(C), "cost-parameter C for C_SVC, EPSILON_SVR and NU_SVR")
    ("svmEps", po::value(&eps)->default_value(eps), "stopping criteria, usually use 0.00001 in nu-SVC, 0.001 in others")
    ("svmP", po::value(&p)->default_value(p), "for EPSILON_SVR")
    ("svmShrinking", po::value(&shrinking)->default_value(shrinking), "use the shrinking heuristics {0,1}")
    ("svmProbability", po::value(&probability)->default_value(probability), "do probability estimates {0,1}")
    ("svmClassWeights", po::value(&classWeights)->multitoken(), "specify if classes should get different weights (default 1)")
    ("svmBalanceTraining", po::value(&weightBalanceTrainingSet)->zero_tokens()->default_value(weightBalanceTrainingSet), "activate this for unbalanced training sets to automatically determine (overwrite) svmClassWeights")
    ;
}

SvmClassifier::SvmClassifier(struct svm_problem *problem, unsigned int featureDimensionality, ClassLabelHistogram *classHistogram, svm_parameter *params)
{
  prob.l = 0;
  prob.y = 0;
  prob.x = 0;
  dimensionality = featureDimensionality;
  classHistogram->getListOfClassLabels(inserter(classLabelList,classLabelList.begin())); // set up classLabelList
  svm = svm_train(problem, params);
  unsigned int classCount = classLabelList.size();
  testNode = new svm_node[dimensionality+1]; // allocate space for classification. +1 because sparse list ends with additional entry "-1"
  decValues = new double[classCount*(classCount-1)/2]; // allocate space for classification.
}

SvmClassifier::SvmClassifier(string filename)
{
  prob.l = 0;
  prob.y = 0;
  prob.x = 0;
  svm = 0;
  ifstream inFile (filename.c_str());
  if (! inFile.eof() )
  {
    string s;
    getline(inFile,s);
    //cout << s << endl;
    char classlabels [255];
    if (sscanf(s.c_str(), "SVMClassifier;%u%s", &dimensionality, classlabels) == 2) {
      string remainingLabels = classlabels;
      while (remainingLabels.length() > 0) {
        int label;
        int count = sscanf(remainingLabels.c_str(), ";%d%s", &label, classlabels);
        if (count >= 1)
          classLabelList.push_back(label);
        if (count == 2)
          remainingLabels = classlabels;
        else
          remainingLabels = "";
      }
    } else { // scanf first line
      throw runtime_error("could not load SVM classifier from file");
    }
  } // !eof
  if (! inFile.eof() ) // read svm model
  {
    string s;
    getline(inFile,s);
    char fn [256];
    if (sscanf(s.c_str(), "svmmodel=%s", fn) == 1) {
      svm = svm_load_model(fn);
    }
  }
  inFile.close();
  unsigned int classCount = classLabelList.size();
  testNode = new svm_node[dimensionality+1]; // allocate space for classification. +1 because sparse list ends with additional entry "-1"
  decValues = new double[classCount*(classCount-1)/2]; // allocate space for classification.
}


SvmClassifier::~SvmClassifier ( )
{
  if (decValues)
    delete decValues;
  if (testNode)
    delete testNode;
  if (svm)
    svm_destroy_model(svm);
  if (prob.x) {
    for (int i=0; i<prob.l; i++) {
      delete [] prob.x[i];
    }
    delete [] prob.x;
  }
  if (prob.y) {
    delete [] prob.y;
  }
}



/***************************************************************************
                         public Methods
 ***************************************************************************/

void SvmClassifier::writeToFile(string filename) const
{
  ofstream outfile;
  outfile.open(filename.c_str(), ios::out);
  outfile << "SVMClassifier;" << dimensionality;
  for(list<ClassLabel>::const_iterator p=classLabelList.begin(); p!=classLabelList.end(); ++p)
    outfile << ";" << *p;
  outfile << endl;
  if (svm != 0) {
    filename.append(".svmmodel");
    svm_save_model(filename.c_str(), svm);
    outfile << "svmmodel=" << filename << endl;
  }
  outfile.close();
}

void SvmClassifier::classify (FeatureVector &fVec, ClassLabelHistogram &hist)
{
  ClassLabelHistogram currFeatureVotes;
  unsigned int classCount = classLabelList.size();
  convertFeatureVector(&fVec, testNode);
  svm_predict_values(svm, testNode, decValues);
  //TODO: With linear kernel and few feature dimensions compared to nb of support vectors use extracted hyperplanes for classification -> faster
  unsigned int pos = 0;
  for (unsigned int i=0;i<classCount;i++) {
    for (unsigned int j=i+1;j<classCount;j++,pos++) {
      if (decValues[pos] > 0) // this is how svm_predict works
        currFeatureVotes.increaseCountOfLabel(svm->label[i],1);
      else
        currFeatureVotes.increaseCountOfLabel(svm->label[j],1);
      // this "smart" solution might give more accurate results if distance to hyperplane was used
//      currFeatureVotes.increaseCountOfLabel(svm->label[i],decValues[pos]);
//      currFeatureVotes.increaseCountOfLabel(svm->label[j],-decValues[pos]);
      //TODO: Missing: use dist to hyperplane instead of decVal: dist = decVal / |w|
      //FAQ: The distance is |decision_value| / |w|. We have |w|^2 = w^Tw = alpha^T Q alpha = 2*(dual_obj + sum alpha_i). Thus in svm.cpp please find the place where we calculate the dual objective value (i.e., the subroutine Solve()) and add a statement to print w^Tw.
      //Forum1: According to myself ;-) the "nicest" hack is to extend the struct SolutionInfo with a field to add '|w|^2' and then in the solve_c_svc, calculate the norm of this vector by using the si->obj and the value of the variable sum_alpha. Maybe you can extend the model files, so that you can save it after training.
    }
  }
  if ((ClassLabel)svm_predict(svm, testNode) != currFeatureVotes.getRealLabelWithMaxWeight()) {
    cerr << "svm classified as " << svm_predict(svm, testNode) << " but histogram as " << currFeatureVotes.getRealLabelWithMaxWeight() << ", histogram: " << currFeatureVotes << endl;
  }
  hist.integrate(currFeatureVotes);
}

vector<SvmClassifier::hyperplane> SvmClassifier::getHyperplanes() const
{
  const svm_model *model = svm;

  if (model->param.kernel_type != LINEAR)
    throw logic_error("SvmClassifier::getHyperplanes only callable if linear kernel was selected");
    
  vector<hyperplane> result;
  int i;
  int nr_class = model->nr_class;
  //int l = model->l;
  
  int *start = (int*)malloc((nr_class)*sizeof(int));
  start[0] = 0;
  for(i=1;i<nr_class;i++)
    start[i] = start[i-1]+model->nSV[i-1];

  int pos=0;                      // position index in resulting w-array
  for(i=0;i<nr_class;i++)         // for each classifier for class i
    for(int j=i+1;j<nr_class;j++) //                 against class j
    {
      // classifier (i,j): coefficients with
      // i are in sv_coef[j-1][nz_start[i]...],
      // j are in sv_coef[i][nz_start[j]...]
      int si = start[i];
      int sj = start[j];
      int ci = model->nSV[i];
      int cj = model->nSV[j];
      double *coef1 = model->sv_coef[j-1];
      double *coef2 = model->sv_coef[i];
      
      hyperplane plane;
      plane.pos = model->label[i];
      plane.neg = model->label[j];
      plane.b = - model->rho[pos];
      //cout << endl << " - processing " << plane.pos << " against " << plane.neg << ":";
      double *w=new double[dimensionality];
      for(unsigned int k=0;k<dimensionality;k++)
        w[k] = 0;
      //cout << endl << " -- processing " << ci << " positive support vectors";
      for(int k=0;k<ci;k++) { // for all SV of class i
        svm_node *sv = model->SV[si+k]; 
        for (int l=0; (sv[l].index >= 0 && sv[l].index < (int)dimensionality); ++l)
          w[sv[l].index] += coef1[si+k] * sv[l].value;
      }
      //cout << endl << " -- processing " << cj << " negative support vectors";
      for(int k=0;k<cj;k++) { // for all SV of class j
        svm_node *sv = model->SV[sj+k]; 
        for (int l=0; (sv[l].index >= 0 && sv[l].index < (int)dimensionality); ++l)
          w[sv[l].index] += coef2[sj+k] * sv[l].value;
      }
      //convert w into vector
      for(unsigned int k=0;k<dimensionality;k++) {
        plane.w.push_back(w[k]);
      }
      result.push_back(plane);
      pos++;
      delete[] w;
    }
  free(start);
  return result;
}

/***************************************************************************
                         private Methods
 ***************************************************************************/

/*!
  \brief Method for converting FeatureVectors into svm_problem

  ATTENTION: this method allocates memory for parameters "prob" and "x_space".
             memory must be deallocated later by using free(prob.x), free(prob.y) and free(x_space)

  \param featureVec the list with FeatureVectors
*/
void SvmClassifier::convertFeatureVectors (list<FeatureVector*>* featureVec) {
  /* struct svm_problem
    int l;
    double *y;
    struct svm_node **x;
  */
  int fCount = featureVec->size();
  prob.l = fCount;
  prob.y = new double [fCount];
  prob.x = new svm_node* [fCount];
  int i = 0;
  for(list<FeatureVector*>::const_iterator featI=featureVec->begin(); featI!=featureVec->end(); ++featI) {
    FeatureVector* fVec = *featI;
    int compactSize = 0;
    for (unsigned int j=0; j<dimensionality; j++) {
      if (fVec->getElement(j) != 0) {
        compactSize++;
      }
    }
    prob.y[i] = (fVec->getTag<TagTrueClassLabel>()).label;
    prob.x[i] = new svm_node[compactSize+1];
    convertFeatureVector(fVec, prob.x[i]);
    i++;
  }
}

/*!
  \brief Method for converting one FeatureVector into a svm_node
  
  ASSUMPTION: enough memory for svm_nodes is already allocated (FEATURE_NUMBER+1)

  \param fVec pointer to the FeatureVector
  \param x pointer to the memory where the data will be stored (must already be allocated!!!)
 */
void SvmClassifier::convertFeatureVector (FeatureVector* fVec, struct svm_node *x) {
  /* struct svm_node
  int index;
  double value;
  */
  //cout << endl << "processing feature vector of size" << fVec->getSize() << flush;
  int fullSize = min(dimensionality,fVec->getSize());
  for (int i=0; i<fullSize; i++) {
    if (fVec->getElement(i) != 0) {
      x->index = i;
      x->value = fVec->getElement(i);
      x++;
    } else {
      //cout << " " << i << "X";
    }
  }
  x->index = -1;
}

