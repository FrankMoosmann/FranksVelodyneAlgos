/*
 * ConsoleProgressBar.cpp
 *
 *  Created on: Feb 16, 2010
 *      Author: moosmann
 */

#include "ConsoleProgressBar.hpp"
#include <cmath>
#include <iostream>
#include <boost/format.hpp>

using namespace std;

const char ConsoleProgressBar::rotChars[4] = {'\\', '|', '/', '-'};

ConsoleProgressBar::ConsoleProgressBar(unsigned int finishAtCount, unsigned int barLength_)
  : barLength(min(barLength_,253u))
   ,maxVal(finishAtCount)
   ,currVal(0)
   ,currRotIdx(0)
{
  buffer[barLength+2] = 0; //end-character
  updateDisplay(false);
}

ConsoleProgressBar::~ConsoleProgressBar() {

}

void ConsoleProgressBar::increment() {
  if (currVal < maxVal)
    ++currVal;
  updateDisplay(true);
}

void ConsoleProgressBar::finish() {
  currVal = maxVal;
  updateDisplay(true);
}

void ConsoleProgressBar::refresh() {
  updateDisplay(true);
}

void ConsoleProgressBar::operator() (int nb) {
  currVal = min(currVal+nb, maxVal);
  updateDisplay(true);
}

void ConsoleProgressBar::updateDisplay(bool clear) {
  unsigned int pos = barLength*currVal/maxVal; // relative position within the brackets [0..max]
  buffer[0] = '[';
  for (unsigned int i=0; i<pos; ++i)
    buffer[i+1] = '=';
  currRotIdx = (currRotIdx+1)%4;
  buffer[pos+1] = rotChars[currRotIdx];
  for (unsigned int i=pos+1; i<barLength; ++i)
    buffer[i+1] = ' ';
  buffer[barLength+1] = ']';
  float percent = 100.0f*(float)currVal/(float)maxVal;
//  cout << boost::format("\r%s %5.1f%%") % buffer % percent << flush;
  if (clear)
    for (unsigned int i=0; i<(barLength+2+1+5+1);++i) // + 2 brackets + space + 5 digits + "%"
      cout << "\b";
  cout << boost::format("%s %5.1f%%") % buffer % percent << flush;
}
