/*
 * FeatureSet.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: moosmann
 */

#include "FeatureSet.hpp"

#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <boost/typeof/std/utility.hpp> //BOOST_AUTO

#include <TagTrueClassLabel.hpp>

using namespace std;

FeatureSet::FeatureSet() {
}

FeatureSet::~FeatureSet() {
}

void FeatureSet::clear() {
  featVecList.clear();
}

void FeatureSet::push_back(FeatureVector f) {
  featVecList.push_back(f);
}

FeatureSet FeatureSet::clone() const
{
  FeatureSet other;
  for (FeatureSet::const_iterator i = begin(); i != end(); ++i) {
    FeatureVector cv = i->clone();
    // check that all data is the same:
    if (cv.size() != i->size())
      throw std::runtime_error("ssssssssss");
    for (unsigned int idx=0; idx<cv.size(); ++idx) {
      if (cv[idx] != (*i)[idx]) {
        cout << "(" << cv[idx] << "!=" << (*i)[idx] << ")" << flush;
        throw std::runtime_error("??????????");
      }
    }
    if (cv.getTag<TagTrueClassLabel>().label != i->getTag<TagTrueClassLabel>().label) {
      cout << "(" << cv.getTag<TagTrueClassLabel>().label << "!=" << i->getTag<TagTrueClassLabel>().label << ")" << flush;
      throw std::runtime_error("!!!!!!!!!!!");
    }
    other.push_back(cv);
  }
  return other;
}

void FeatureSet::splice(FeatureSet &other) {
  featVecList.splice(featVecList.end(), other.featVecList);
}

void FeatureSet::splitUniform(double splitFactor, FeatureSet &second) {
  double counter = 0.0;
  BOOST_AUTO(fvi, featVecList.begin());
  while (fvi != featVecList.end()) {
    counter += splitFactor;
    if (counter >= 1.0) {
      counter -= 1.0;
      second.featVecList.push_back(*fvi);
      fvi = featVecList.erase(fvi);
      // list::splice might be more efficient but invalidates the iterators -> more complicated code
    } else {
      ++fvi;
    }
  }
}

void FeatureSet::saveToFile(string fileName) {
  ofstream datafile(fileName.c_str());
  BOOST_FOREACH(const FeatureVector &fv, featVecList) {
    for (unsigned int i=0; i<fv.size(); ++i) {
      datafile << fv[i] << " ";
    }
    if (fv.hasTag<TagTrueClassLabel>()) {
      datafile << fv.getTag<TagTrueClassLabel>().label << " ";
    }
    datafile << endl;
  }
}

//void FeatureVector::writeIntoStream(list<FeatureVector*> &featureList, ostream &data, ostream &info, int maxcount)
//{
//  //\TODO: write Info-Objects into stream too!!!
//  int count = 0;
//  for(list<FeatureVector*>::const_iterator q=featureList.begin(); (q!=featureList.end()); ++q)
//  {
//    FeatureVector *featI = *q;
//    //featurelist << featI->windowPosX << " " << featI->windowPosY << " " << featI->windowSize << " 0 ";
//    for(unsigned int i = 0; i<featI->getSize(); i++)
//      data << featI->getElement(i) << " ";
//    data << endl;
//    cerr << "Writing Info-Objects of Features not supported yet" << endl;
////    info << featI->classLabel << " ";
////    info << featI->windowToObjectRatio << " ";
////    info << featI->windowPosX << " ";
////    info << featI->windowPosY << " ";
////    info << featI->windowSize << " ";
////    info << featI->imageWidth << " ";
////    info << featI->imageHeight << " ";
//      info << endl;
//    //cout << "." << flush;
//    count++;
//    if ((maxcount > 0) && (count >= maxcount))
//      break;
//  }
//}
//
//
//int FeatureVector::readFromStream(istream &data, istream &info, list<FeatureVector*> &featureList, int nbToRead, std::map<unsigned int,bool> *featureDimensionsIgnored)
//{
//  //\TODO: read Info-Objects from stream too!!!
//  string sdata;
//  string sinfo;
//  int nbRead = 0;
//  //cout << endl << "reading...";
//  while ( ((nbRead < nbToRead) || (nbToRead == 0)) && (!data.eof()) && (!info.eof())) {
//    getline(data,sdata);
//    getline(info,sinfo);
//    FeatureVector *fv = new FeatureVector();
//
//    istringstream ssdata(sdata);
//    FeatureData value;
//    int i=0;
//    while (!ssdata.eof()) {
//      ssdata >> value;
//      if (ssdata.fail())
//        break;
//      else {
//        if ((!featureDimensionsIgnored) || (!(*featureDimensionsIgnored)[i])) {
//          fv->addElement(value);
//        } else {
//          //cout << "skipping dimension "<<i<<endl;
//        }
//      }
//      i++;
//    }
//    TagImagePos tg;
//    cerr << "Reading Info-Objects of Features not supported yet" << endl;
////    if (sscanf(sinfo.c_str(), "%d %lf %u %u %u %u %u ",
////                            &(fv->classLabel),
////                            &(tg.windowToObjectRatio),
////                            &(tg.windowPosX),
////                            &(tg.windowPosY),
////                            &(tg.windowSize),
////                            &(tg.imageWidth),
////                            &(tg.imageHeight))
////        >= 1) {  //at least classLabel must have been read successfully
////      fv->setTag<TagImagePos>(tg);
//      featureList.push_back(fv);
////      nbRead++;
////    } else {
////      delete fv;
////    }
//  }
//  //cout << noRead;
//  return nbRead;
//}


//int FeatureVector::readFeatureVectorsFromFile (string fileName, string directory, list<FeatureVector*> &featureList, int nb, std::map<unsigned int,bool> *featureDimensionsIgnored) {
//  //\TODO: test this new converstion via boost::filesystem!!!
//  using namespace boost::filesystem;
//  path filename(fileName);
//  path basename = path(directory) / filename.file_string();
//  string dataname = basename.string() + ".desc";
//  string infoname = basename.string() + ".info";
//  //string dataname = directory + extractFileName(fileName) + ".desc";
//  //string infoname = directory + extractFileName(fileName) + ".info";
//  //cout << "reading file " << dataname << " / .info" << endl;
//  ifstream dataFile(dataname.c_str());
//  ifstream infoFile(infoname.c_str());
//  int result = 0;
//  if ((dataFile) && (infoFile)) {
//    result = FeatureVector::readFromStream(dataFile, infoFile, featureList, nb, featureDimensionsIgnored);
//  }
//  dataFile.close();
//  infoFile.close();
//  return result;
//}
//
//void FeatureVector::saveFeatureVectorsToFile(string fileName, list<FeatureVector*> &featureList, string directory) {
//  //\TODO: test this new conversion via boost::filesystem!!!
//  using namespace boost::filesystem;
//  path filename(fileName);
//  path basefilename = path(directory) / filename.file_string();
//  string dataname = basefilename.string() + ".desc";
//  cout << endl << "saving features to " << dataname;
//  string infoname = basefilename.string() + ".info";
//  cout << endl << "saving info to " << infoname;
//  //string dataname = directory + extractFileName(fileName) + ".desc";
//  //string infoname = directory + extractFileName(fileName) + ".info";
//  ofstream datafile(dataname.c_str(), ios::out);
//  ofstream infofile(infoname.c_str(), ios::out);
//  FeatureVector::writeIntoStream(featureList, datafile, infofile);
//  datafile.close();
//  infofile.close();
//}
