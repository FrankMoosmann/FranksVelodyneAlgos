// -*- c++ -*-
//------------------------------------------------------------------------------
//                              TimeVal.cpp
//------------------------------------------------------------------------------
//  Copyright (c) 2000 by Vladislav Grinchenko
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Library General Public
//  License as published by the Free Software Foundation; either
//  version 2 of the License, or (at your option) any later version.
//------------------------------------------------------------------------------
//  Created: 09/28/99
//------------------------------------------------------------------------------
#include <time.h>    // localtime(3), gmtime(3)
#include <stdio.h>    // sprintf(3)
#include <cstring>

#include "TimeVal.hpp"

#if defined (WIN32)
#   include <windows.h>
#endif

using namespace ASSA;

//------------------------------------------------------------------------------
// Static definitions
//------------------------------------------------------------------------------

TimeVal TimeVal::m_zero;       // zero time
static const long ONE_SECOND = 1000000;

#ifndef __GNUC__
#define EPOCHFILETIME (116444736000000000i64)
#else
#define EPOCHFILETIME (116444736000000000LL)
#endif

//------------------------------------------------------------------------------
// Wrappers
//------------------------------------------------------------------------------
TimeVal 
TimeVal::
gettimeofday () 
{ 
  timeval tv;
  
#ifdef WIN32
    FILETIME        ft;
    LARGE_INTEGER   li;
    __int64         t;
    static int      tzflag;

    GetSystemTimeAsFileTime(&ft);
    li.LowPart  = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;
    t  = li.QuadPart;       /* In 100-nanosecond intervals */
    t -= EPOCHFILETIME;     /* Offset to the Epoch time */
    t /= 10;                /* In microseconds */
    tv.tv_sec  = (long)(t / 1000000);
    tv.tv_usec = (long)(t % 1000000); 
#else
  ::gettimeofday (&tv, 0);
#endif
  return tv;
}

//------------------------------------------------------------------------------
// Membe functions
//------------------------------------------------------------------------------

TimeVal&
TimeVal::
operator+=(const TimeVal& rhs_)
{
  tv_sec += rhs_.tv_sec;
  tv_usec += rhs_.tv_usec;

  if (tv_usec >= ONE_SECOND) {
    tv_usec -= ONE_SECOND;
    tv_sec++;
  } 
  else if (tv_sec >= 1 && tv_usec < 0) {
    tv_usec += ONE_SECOND;
    tv_sec--;
  }
  normalize ();
  return *this;
}

TimeVal&
TimeVal::
operator-=(const TimeVal& rhs_)
{
  tv_sec -= rhs_.tv_sec;
  tv_usec -= rhs_.tv_usec;

  if (tv_usec < 0) {
    tv_usec += ONE_SECOND;
    tv_sec--;
  } 
  else if (tv_usec >= ONE_SECOND) {
    tv_usec -= ONE_SECOND;
    tv_sec++;
  }
  normalize ();
  return *this;
}

void
TimeVal::
normalize ()
{
  if (tv_usec >= ONE_SECOND) {
    do {
      tv_sec++;
      tv_usec -= ONE_SECOND;
    }
    while (tv_usec >= ONE_SECOND);
  }
  else if (tv_usec <= -ONE_SECOND) {
    do {
      tv_sec--;
      tv_usec += ONE_SECOND;
    }
    while (tv_usec <= -ONE_SECOND);
  }

  if (tv_sec >= 1 && tv_usec < 0) {
    tv_sec--;
    tv_usec += ONE_SECOND;
  }
  else if (tv_sec < 0 && tv_usec > 0) {
    tv_sec++;
    tv_usec -= ONE_SECOND;
  }
}


//------------------------------------------------------------------------------
// All possible variation of HH:MM:SS.MMM I could think of:
//------------------------------------------------------------------------------

string 
TimeVal::
fmtString (const char* fmt_) const
{
  struct tm ct;
  char buf[80];
  memset (buf, 0, 80);

  if (m_tz == gmt)
    ct = *( localtime ((const time_t*) &tv_sec) );
  else
    ct = *( gmtime ((const time_t*) &tv_sec) );

  if (fmt_ == NULL) {
    strftime (buf, 80, "%Y/%j %H:%M:%S", &ct);
    sprintf (buf + strlen(buf),
       ".%03ld", (tv_usec %1000000)/1000);
  }
  else {
    strftime(buf, 80, fmt_, &ct);
  }
  return string (buf);
}

string
TimeVal::
fmt_hh_mm_ss_mls () const
{
    struct tm ct;
    char buf [80];
  memset (buf, 0, 80);

    if (m_tz == gmt)
        ct = *( localtime ((const time_t*) &tv_sec) );
    else
        ct = *( gmtime ((const time_t*) &tv_sec) );

  strftime (buf, 80, "%H:%M:%S", &ct);
  sprintf (buf + strlen(buf), ".%03ld", millisec ());

    return string (buf);
}

string
TimeVal::
fmt_mm_ss_mls () const
{
    struct tm ct;
    char buf [80];
  memset (buf, 0, 80);

    if (m_tz == gmt)
        ct = *( localtime ((const time_t*) &tv_sec) );
    else
        ct = *( gmtime ((const time_t*) &tv_sec) );

  strftime (buf, 80, "%M:%S", &ct);
  sprintf (buf + strlen(buf), ".%03ld", millisec ());

    return string (buf);
}

string
TimeVal::
fmt_ss_mls () const
{
    struct tm ct;
    char buf [80];
  memset (buf, 0, 80);

    if (m_tz == gmt)
        ct = *( localtime ((const time_t*) &tv_sec) );
    else
        ct = *( gmtime ((const time_t*) &tv_sec) );

  strftime (buf, 80, "%S", &ct);
  sprintf (buf + strlen(buf), ".%03ld", millisec ());

    return string (buf);
}


