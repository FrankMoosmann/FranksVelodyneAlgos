/*! \file kogmo_time.h
 * \brief Definition of Time Handling Types and Functions (C-Interface)
 *
 * Copyright (c) 2005,2006 Matthias Goebl <mg*tum.de>
 *     Lehrstuhl fuer Realzeit-Computersysteme (RCS)
 *     Technische Universitaet Muenchen (TUM)
 */

/*! \defgroup kogmo_time Time Handling
 * \brief Types and Functions for Time Handling in Cognitive Automobiles.
 */
/*@{*/

// In this piece of source code, I tried to follow the GNU Coding Standards,
// see http://www.gnu.org/prep/standards/standards.html
// This file could also be used as a reference for doxygen use.


#ifndef KOGMO_TIME_H
#define KOGMO_TIME_H

#include <inttypes.h> /* int64_t */

#ifdef __cplusplus
 namespace KogniMobil {
 extern "C" {
 #define _const const   //!< constant argument for inclusion in C++
#else /* no C++ */
 #define _const         //!< constant argument, but irrelevant in C
#endif



/*! \brief Absolute Timestamp (in Ticks since the "Epoch").
 *
 * This is the reference type that to be used for all timestamps. \n
 * It can represent points in time from the year 1970 until 2262
 * with a nanosecond resolution.
 *
 * Background:\n\n
 * 
 * The Epoch starts 00:00:00 UTC, January 1, 1970, see also
 * (POSIX.1 and time(2), http://de.wikipedia.org/wiki/Epoch). \n
 * This is a signed 64-Bit value that will overflow in the year 2262. \n
 * ( 2^63/1000000000/60/60/24/365 = 292.471.. )
 *
 * Why a 64-Bit value value? Aren't 32 Bits enough?\n
 * A 32-Bit (unsigned) value would overflow after
 *  - less than 5 seconds with a nanosecond resolution
 *    ( 2^32/1000/1000/10e00=4.294... )
 *  - less than 2 hours with a microsecond resolution
 *    ( 2^32/1000/60/60/1000=1.193... )
 *  - less than 50 days with millisecond resolution
 *    ( 2^32/1000/60/60/24 = 49.710.. )
 *
 * So only a millisecond granularity would be usefull, but
 *  - we can use it only for a relative time with a defined
 *    begin (e.g. system boot)
 *  - We have to remember this start time (and possibly
 *    save it in out logs, maybe communicate it to other vehicles)
 *  - We wouldn't win a "Around the World in Eighty Days" race ;-)
 *
 * So please
 * - use this timestamp type
 * - get the timestamp via kogmo_timestamp_now()
 * - work with the timestamps using kogmo_timestamp_difference_in_seconds()
 *   and kogmo_timestamp_add_seconds()
 * - use KOGMO_TIMESTAMP_TICKSPERSECOND just for your information
 * - DO NOT implement own calculations with nanoseconds
 *
 * There are convenience functions with float-times in seconds,
 * these can also be used in calculations.\n
 * For exchanging times you should use kogmo_time_t.
 *
 * For any questions, please contact me!
 * Matthias Goebl <mg*tum.de>
 */
typedef int64_t kogmo_timestamp_t;


/*! \brief Ticks per Second.
 * Value, by which the result from kogmo_timestamp_now() will be incremented
 * every second.
 *
 * As we use Nanoseconds, its equal to 1000000000. \n
 * But *PLEASE* dont hardcode 1000000000 anywere, use this constant if necessary!\n
 * Better use the following functions to handle timestamps.
 */
#define KOGMO_TIMESTAMP_TICKSPERSECOND 1000000000
#define KOGMO_TIMESTAMP_NANOSECONDSPERTICK (1000000000/KOGMO_TIMESTAMP_TICKSPERSECOND)


/*! \brief Get absolute Timestamp for current Time.
 * \returns zero on errors
 */
kogmo_timestamp_t
kogmo_timestamp_now (void);


/*! \brief Calcute difference between two Timestamps (in Nanoseconds).
 *
 * \param ts_begin First Timestamp, begin of the interval
 * \param ts_end Second Timestamp, end of the interval
 * \returns Length of the interval in nanoseconds (64-Bit integer value)
 */
int64_t
kogmo_timestamp_diff_ns (kogmo_timestamp_t ts_begin,
                         kogmo_timestamp_t ts_end);


/*! \brief Calcute difference between two Timestamps (in Seconds, floating point).
 *
 * This may simplify your calculations.\n
 * For exchanging timestamps you should use kogmo_timestamp_t instead.
 *
 * \param ts_begin First Timestamp, begin of the interval
 * \param ts_end Second Timestamp, end of the interval
 * \returns Length of the interval in Seconds (floating point value)
 */
double
kogmo_timestamp_diff_secs (kogmo_timestamp_t ts_begin,
                           kogmo_timestamp_t ts_end);


/*! \brief Add time offset to Timestamps (in Nanoseconds).
 *
 * \param ts Timestamp
 * \param ns Nanoseconds to add, can be negative (64-Bit integer value)
 * \returns New Timestamp
 */
kogmo_timestamp_t
kogmo_timestamp_add_ns (kogmo_timestamp_t ts, int64_t ns);


/*! \brief Add time offset to Timestamps (in Seconds, floating point).
 *
 * \param ts Timestamp
 * \param secs Seconds to add, can be negative (floating point value)
 * \returns New Timestamp
 */
kogmo_timestamp_t
kogmo_timestamp_add_secs (kogmo_timestamp_t ts, double secs);


/*! \brief Length of the Output of kogmo_timestamp_to_string()
 */
#define KOGMO_TIMESTAMP_STRINGLENGTH 30

/*! \brief Type of a pre-allocated String to receive the Output of
 * kogmo_timestamp_to_string()
 */
typedef char kogmo_timestamp_string_t[KOGMO_TIMESTAMP_STRINGLENGTH];


/*! \brief Convert an Absolute Timestamp to a String,
 * formated according to the ISO-Standard.
 *
 * Format: 2006-04-26 23:54:27.832241000 (YYYY-MM-DD hh:mm:ss.sssssssss)
 * \param ts Timestamp
 * \param str Pre-Allocated String for the Output
 * \returns non-zero on errors
 *
 * The format follows ISO 8601 "Data elements and interchange formats,
 * Information interchange, Representation of dates and times",
 * with the addition of subseconds as ".sssssssss". \n
 * See http://en.wikipedia.org/wiki/ISO_8601 and
 * http://hydracen.com/dx/iso8601.htm.
 */
int
kogmo_timestamp_to_string (kogmo_timestamp_t ts, kogmo_timestamp_string_t str);


/*! \brief Convert an ISO-Time-String to an Absolute Timestamp.
 *
 * \param str Time as String, see kogmo_timestamp_to_string() for format
 * \returns Timestamp, zero on errors
 *
 * For a incomplete ISO8601-Time the missing elements will be assumed to 0.\n
 * For your comfort, str can also be a numeric timestamp in a string (>9999).
 */
kogmo_timestamp_t
kogmo_timestamp_from_string (_const char *str);



#ifdef __cplusplus
}; /* extern "C" */
}; /* namespace KogniMobil */
#endif

#endif /* KOGMO_TIME_H */
/*@}*/
