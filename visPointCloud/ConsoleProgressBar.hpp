/*
 * ConsoleProgressBar.h
 *
 *  Created on: Feb 16, 2010
 *      Author: moosmann
 */

#ifndef CONSOLEPROGRESSBAR_H_
#define CONSOLEPROGRESSBAR_H_

class ConsoleProgressBar {
public:
  ConsoleProgressBar(unsigned int finishAtCount, unsigned int barLength = 40);
  virtual ~ConsoleProgressBar();

  void increment();
  void finish();
  void refresh();
  void operator() (int nb=1);

private:
  static const char rotChars[4];
  const unsigned int barLength; // number of characters within []
  const unsigned int maxVal; // target value
  unsigned int currVal; // value of progress
  unsigned int currRotIdx; // index for rotating dash
  char buffer[256];
  char percent[7];

  void updateDisplay(bool clear); // set clear=false the first time it is called
};

#endif /* CONSOLEPROGRESSBAR_H_ */
