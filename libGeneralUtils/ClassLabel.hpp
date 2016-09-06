#ifndef CLASS_LABEL_H
#define CLASS_LABEL_H

/*!
  \file   ClassLabel.h
  \brief  Provides a definitions for a class label (for classification purposes)
  \author  Frank Moosmann (<frank.moosmann@kit.edu>)
  \date    07.11.2005
*/

#include <climits>
#include <string>

typedef int ClassLabel;
const int NULL_CLASS_LABEL = INT_MAX;

//! for visualization purposes, associate one color to each class. returns lower case string
std::string getClassColorString(const ClassLabel c);

//! for visualization purposes, associate one color to each class. returns RGB-values in range [0..255]
void getClassColorRGB(const ClassLabel c, unsigned int &R, unsigned int &G, unsigned int &B);

#endif /*CLASS_LABEL_H*/
