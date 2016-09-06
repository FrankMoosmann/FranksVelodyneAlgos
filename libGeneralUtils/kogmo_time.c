/* KogMo-RTDB: Real-time Database for Cognitive Automobiles
 * Copyright (c) 2003-2007 Matthias Goebl <matthias.goebl*goebl.net,mg*tum.de>
 *     Lehrstuhl fuer Realzeit-Computersysteme (RCS)
 *     Technische Universitaet Muenchen (TUM)
 * Licensed under the GNU Lesser General Public License v3, see the file COPYING.
 */
/*! \file kogmo_time.c
 * \brief Implementation of Time Handling Functions
 *
 * Copyright (c) 2005,2006 Matthias Goebl <mg*tum.de>
 *     Lehrstuhl fuer Realzeit-Computersysteme (RCS)
 *     Technische Universitaet Muenchen (TUM)
 */

// In this piece of source code, I tried to follow the GNU Coding Standards,
// see http://www.gnu.org/prep/standards/standards.html


#include <stdio.h>  /* printf */
#include <string.h> /* strlen strncpy memset */
#include <stdlib.h> /* strtoul */
#include <time.h>
#include "kogmo_time.h"

#ifdef MACOSX
#include <sys/time.h>
#include <mach/clock.h>

kogmo_timestamp_t
kogmo_timestamp_now (void)
{
  kogmo_timestamp_t ts; 
  struct timeval tval;
  int err;
  err = gettimeofday ( &tval, NULL);
  if (err != 0)
    return 0;
  ts = (kogmo_timestamp_t) tval.tv_sec*KOGMO_TIMESTAMP_TICKSPERSECOND
       + tval.tv_usec*1000/KOGMO_TIMESTAMP_NANOSECONDSPERTICK;
  return ts;
} 

#else

kogmo_timestamp_t
kogmo_timestamp_now (void)
{
  kogmo_timestamp_t ts;
  struct timespec tspec;
  int err;
  // TODO: This gives unix time, but according to definition in kogmo_time.h
  // this should be TAI, without leap seconds.
  // TODO: use CLOCK_MONOTONIC with offset by CLOCK_REALTIME.
  // But first we would need a concept for a inter-vehicle time base.
  // Most likely, the vehicles will sync their time to GPS and there
  // will be no time skips.
  err = clock_gettime (CLOCK_REALTIME, &tspec);
  if (err != 0)
    return 0;
  ts = (kogmo_timestamp_t) tspec.tv_sec*KOGMO_TIMESTAMP_TICKSPERSECOND
       + tspec.tv_nsec/KOGMO_TIMESTAMP_NANOSECONDSPERTICK;
  return ts;
} 
#endif


inline int64_t
kogmo_timestamp_diff_ns (kogmo_timestamp_t ts_begin,
                         kogmo_timestamp_t ts_end)
{
  return (ts_end-ts_begin) * KOGMO_TIMESTAMP_NANOSECONDSPERTICK;
}


inline double
kogmo_timestamp_diff_secs (kogmo_timestamp_t ts_begin,
                           kogmo_timestamp_t ts_end)
{
  return (double)(ts_end-ts_begin) / KOGMO_TIMESTAMP_TICKSPERSECOND;
}


inline kogmo_timestamp_t
kogmo_timestamp_add_ns (kogmo_timestamp_t ts, int64_t ns)
{
  return ts + ns / KOGMO_TIMESTAMP_NANOSECONDSPERTICK;
}


inline kogmo_timestamp_t
kogmo_timestamp_add_secs (kogmo_timestamp_t ts, double secs)
{
  return ts + secs * KOGMO_TIMESTAMP_TICKSPERSECOND;
}


#ifdef MACOSX
int
kogmo_timestamp_to_string (kogmo_timestamp_t ts, kogmo_timestamp_string_t str)
{
  struct tm timetm;
  time_t secs;
  if ( snprintf( str, sizeof (kogmo_timestamp_string_t), "%lli", ts ) < 0 )
    return(-1);
  return 0;
}

#else
int
kogmo_timestamp_to_string (kogmo_timestamp_t ts, kogmo_timestamp_string_t str)
{
  struct tm timetm;
  time_t secs;
  str[0] = '\0'; // null-termination, to be sure
  secs = ts / KOGMO_TIMESTAMP_TICKSPERSECOND;
  if ( localtime_r( &secs, &timetm) == NULL )
    return(-1);
  if ( strftime( str, sizeof (kogmo_timestamp_string_t) -1,
                 "%Y-%m-%d %H:%M:%S", &timetm ) == 0 )
    return(-1);
  if ( snprintf( str + strlen (str),
                 sizeof (kogmo_timestamp_string_t) - strlen (str),
                 ".%09lli",
                 (long long int) (ts % KOGMO_TIMESTAMP_TICKSPERSECOND) ) < 0 )
    return(-1);
  str[sizeof(kogmo_timestamp_string_t)-1] = '\0'; // to be sure
  return 0;
}
#endif


kogmo_timestamp_t
kogmo_timestamp_from_string (char *str)
{
  kogmo_timestamp_t ts;
  int err;
  struct tm tmx;
  time_t secs;
  char ns_string[10]="";
  char *ns_string_end;
  unsigned long int ns=0;
  int i;
  long long int lli;

  if ( str == NULL )
    return 0;
  memset ( &tmx, 0, sizeof(struct tm));
  ns = 0;
  err = sscanf ( str, "%d-%d-%d%*[ _tT]%d:%d:%d.%9s",
                 &tmx.tm_year, &tmx.tm_mon , &tmx.tm_mday,
                 &tmx.tm_hour, &tmx.tm_min, &tmx.tm_sec,
                 ns_string );
  // we need at least a date (time will be 00:00:00.0 then)
  if ( err < 3 )
    {
      // for even more comfort, this also accepts times relative to now,
      // +/- followed by an offset in secondstimestamps (e.g. +12.3 or -7.25)
      if ( str[0] == '+' || str[0] == '-' )
        {
          double off;
          err = sscanf ( str, "%lf", &off);
          if ( err != 1 )
            return 0;
          ts = kogmo_timestamp_add_secs ( kogmo_timestamp_now (), off);
          return ts;
        }
      // for your comfort, this also accepts timestamps (>9999)
      err = sscanf ( str, "%lli", &lli);
      ts = (kogmo_timestamp_t) lli;
      if ( err != 1 || ts <= 9999)
        return 0;
      return ts;
    }
  // the ranges of those value are a horrible confusion, see mktime(3)
  tmx.tm_year -= 1900;
  tmx.tm_mon -= 1;
  tmx.tm_isdst = -1;
  secs = mktime (&tmx);
  if ( secs < 0 )
    return 0;
  ns = strtoul(ns_string, &ns_string_end, 10);
  // calculate the correct decimal fraction (9 digits)
  // this is: ns *= pow ( 10, 9 - strlen (ns_string) );
  // but prevent dependency from math-library
  for(i=ns_string_end-ns_string;i<9;i++)
    ns *= 10;
  ts = (kogmo_timestamp_t)ns +
       (kogmo_timestamp_t)secs * KOGMO_TIMESTAMP_TICKSPERSECOND;
  return ts;
}
