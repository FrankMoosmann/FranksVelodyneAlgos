// Author: Andreas Geiger, 2011

#ifndef TIMER_HPP
#define TIMER_HPP

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/time.h>

class Timer {
  
public:
  
  Timer() {reset();}
  
  ~Timer() {}
  
  void start (std::string title = "") {
    //std::cout << desc.size()+1 << ") " << title << std::endl;
    stop();
    desc.push_back(title);
    push_back_curr_time(starttimes);
  }
  
  void stop () {
    if (stoptimes.size()<starttimes.size())
      push_back_curr_time(stoptimes);
  }
  
  void plot () {
    stop();
    float total_time = 0;
    size_t maxDescLength = 10;
    for (int32_t i=0; i<(int32_t)desc.size(); i++) {
      maxDescLength = std::max(maxDescLength, desc[i].length());
    }
    for (int32_t i=0; i<(int32_t)desc.size(); i++) {
      float curr_time = getTimeDifferenceMilliseconds(starttimes[i],stoptimes[i]);
      total_time += curr_time;
      std::cout.width(maxDescLength+2);
      std::cout << desc[i] << " ";
      std::cout << std::fixed << std::setprecision(1) << std::setw(6);
      std::cout << curr_time;
      std::cout << " ms" << std::endl;
    }
    //std::cout << "========================================" << std::endl;
    std::cout << "Total time ";
    std::cout << std::fixed << std::setprecision(1) << std::setw(6);
    std::cout << total_time;
    std::cout << " ms" << std::endl;
  }
  
  void reset () {
    desc.clear();
    starttimes.clear();
    stoptimes.clear();
  }
  
private:
  std::vector<std::string>  desc;
  std::vector<timeval>      starttimes;
  std::vector<timeval>      stoptimes;
  
  static void push_back_curr_time (std::vector<timeval> &vec) {
    timeval curr_time;
    gettimeofday(&curr_time,0);
    vec.push_back(curr_time);
  }
  
  static float getTimeDifferenceMilliseconds(timeval a,timeval b) {
    return ((float)(b.tv_sec -a.tv_sec ))*1e+3 +
           ((float)(b.tv_usec-a.tv_usec))*1e-3;
  }
};

#endif
