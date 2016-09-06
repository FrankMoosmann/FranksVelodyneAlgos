// -*- c++ -*-
//------------------------------------------------------------------------------
//                          TimeVal.h
//------------------------------------------------------------------------------
//  Copyright (c) 1999 by Vladislav Grinchenko
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Library General Public
//  License as published by the Free Software Foundation; either
//  version 2 of the License, or (at your option) any later version.
//------------------------------------------------------------------------------
//  Created: 09/28/1999
//------------------------------------------------------------------------------
#ifndef TIME_VAL_H
#define TIME_VAL_H

#include <sys/time.h>    // gettimeofday(3)
#include <unistd.h>    // gettimeofday(3)

#include <string> 
using std::string;

namespace ASSA {

/** @file TimeVal.h
 *
 * Class TimeVal is a wrapper around UNIX timeval structure.
 * 
 */
class TimeVal : public timeval
{
public:
  enum { 
    gmt,          /**< GMT */
    loc            /**< Local Time Zone */
  };

  /** Default constructor. Sets time to 0 sec. 0 usecs. To get 
      current time, use TimeVal (gettimeofday());
   */
  TimeVal ();

  /** Constructor from seconds/microseconds pair.
   */
  TimeVal (long sec_, long msec_);

  /** Constructor from double.
   */
  TimeVal (double d_);

  /** Constructor from <TT> struct timeval</TT>.
   */
  TimeVal (const timeval& tv_);

  /** Copy constructor.
   */
  TimeVal (const TimeVal& tv_); 

  /** Implicit conversion to double.
   */
  operator double () const;

  /// Set seconds.
  void sec (long sec_) { tv_sec = sec_; }

  /// Get secons.
  long sec (void) const { return tv_sec; }

  /// Set microseconds
  void msec (long msec_) { tv_usec = msec_; }

  /// Get microseconds
  long msec (void) const { return tv_usec; }

  /** Convert tv_usec's microseconds (=1/1,000,000 sec) 
    to milliseconds (=1/1,000 sec).
  */
  long millisec () const;

  /// Set timezone.
  void tz (int tz_) { m_tz = tz_; }

  /// Get timezone.
  int tz (void) const { return m_tz; }

  TimeVal& operator= (const TimeVal& tv_);

  /// Addition.
  TimeVal& operator+= (const TimeVal& rhs_);

  /// Substraction.
  TimeVal& operator-= (const TimeVal& rhs_);

  /// Addition.
  friend TimeVal operator+(const TimeVal& lhs_, const TimeVal& rhs_);

  /// Substraction.
  friend TimeVal operator-(const TimeVal& lhs_, const TimeVal& rhs_);

  /// Comparison.
  bool operator< (const TimeVal& rhs_) const;

  /// Equality.
  bool operator==(const TimeVal& rhs_) const;

  /// Comparison.
  friend bool operator> (const TimeVal& lhs_, const TimeVal& rhs_);

  /// Comparison.
  friend bool operator!=(const TimeVal& lhs_, const TimeVal& rhs_);

  /// Comparison.
  friend bool operator<=(const TimeVal& lhs_, const TimeVal& rhs_);

  /// Comparison.
  friend bool operator>=(const TimeVal& lhs_, const TimeVal& rhs_);

  /** Format timeval structure into readable format.
      Default format is CCYY/DDD HH:MM:SS.MMM which is de fasco
      for the software. To get something different, pass fmt_
      format string as specified by strftime(3). Popular format
      is "%c" which will return something like: 
      "Fri Oct 1 10:54:27 1999". Note that timezone aspect of
      formatting time is controlled by tz() member function.
      @param fmt_ Format string as in strftime(3)
      @return Formatted string.
  */
  string fmtString (const char* fmt_ = NULL) const;

  /** Format timeval structure in readable format HH:MM:SS 
   */
  string fmt_hh_mm_ss () const;

  /** Format timeval structure in readable format HH:MM:SS.MLS
   */
  string fmt_hh_mm_ss_mls () const;

  /** Format timeval structure in readable format MM:SS 
   */
  string fmt_mm_ss () const;

  /** Format timeval structure in readable format MM:SS.MLS
   */
  string fmt_mm_ss_mls () const;

  /** Format timeval structure in readable format SS.MLS
   */
  string fmt_ss_mls () const;

  /** Static that returns zero timeval: {0,0}
   */
  static TimeVal zeroTime () { return m_zero; }

  /** Shields off underlying OS differences in getting
      current time.
      @return time of the day as timeval
  */
  static TimeVal gettimeofday ();

protected:
  /// Internal initialization common to most constructors.
  void init (long, long, int);

private:
  /// Normalization after arithmetic operation.
  void normalize ();

private:
  /// Time zone
  int m_tz;

  /// Zero time value
  static TimeVal m_zero;
};
//------------------------------------------------------------------------------
// Inlines
//------------------------------------------------------------------------------

inline void
TimeVal::
init (long s_, long ms_, int tz_)
{
  tv_sec = s_;
  tv_usec = ms_;
  m_tz = tz_;
  normalize ();
}

inline 
TimeVal::
TimeVal ()
{
  init (0, 0, gmt);
}

inline 
TimeVal::
TimeVal (long sec_, long msec_) 
{
  init (sec_, msec_, gmt);
}

inline 
TimeVal::
TimeVal (double d_)
  : m_tz (gmt)
{
  long l = long(d_);
  tv_sec = l;
  tv_usec = (long) ((d_ - double(l))*1000000.0);
  normalize();
}

inline 
TimeVal::
TimeVal (const timeval& tv_)
{
  init (tv_.tv_sec, tv_.tv_usec, gmt);
}

inline 
TimeVal::
TimeVal (const TimeVal& tv_)
{
  init (tv_.tv_sec, tv_.tv_usec, tv_.m_tz);
}

inline 
TimeVal::operator double () const
{ 
  return tv_sec + tv_usec / 1000000.0;
}

inline long
TimeVal::
millisec () const
{
  return (msec () % 1000000) / 1000;
}

inline string
TimeVal::
fmt_hh_mm_ss () const
{
  return fmtString ("%T");
}

inline string
TimeVal::
fmt_mm_ss () const
{
  return fmtString ("%M:%S");
}

//------------------------------------------------------------------------------
// Friend functions
//------------------------------------------------------------------------------

inline TimeVal&
TimeVal::
operator=(const TimeVal& tv_)
{
  init (tv_.tv_sec, tv_.tv_usec, tv_.m_tz);
  return *this;
}

inline TimeVal
operator+(const TimeVal& lhs_, const TimeVal& rhs_)
{
  TimeVal temp(lhs_);
  temp += rhs_;
  temp.normalize ();
  return temp;
}

inline TimeVal
operator-(const TimeVal& lhs_, const TimeVal& rhs_)
{
  TimeVal temp(lhs_);
  temp -= rhs_;
  temp.normalize ();
  return temp;
}

inline bool 
TimeVal::
operator<(const TimeVal& rhs_) const
{
  return (tv_sec < rhs_.tv_sec
    || (tv_sec == rhs_.tv_sec && tv_usec < rhs_.tv_usec) ) ;
}

inline bool 
TimeVal::
operator==(const TimeVal& rhs_) const
{
  return !(*this < rhs_ || rhs_ < *this);
}

inline bool
operator> (const TimeVal& lhs_, const TimeVal& rhs_)
{
  return rhs_ < lhs_;
}

inline bool 
operator!=(const TimeVal& lhs_, const TimeVal& rhs_)
{
  return !( lhs_ == rhs_ );
}

inline bool
operator<=(const TimeVal& lhs_, const TimeVal& rhs_)
{
  return !(rhs_ < lhs_);
}

inline bool
operator>=(const TimeVal& lhs_, const TimeVal& rhs_)
{
  return !(lhs_ < rhs_);
}

} // end namespace ASSA

#endif /* TIME_VAL_H */  
