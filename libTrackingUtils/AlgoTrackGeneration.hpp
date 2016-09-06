#ifndef ALGOTRACKGENERATION_H_
#define ALGOTRACKGENERATION_H_

#include "LidarFrame.hpp"
#include "AlgoRegistration.hpp"

/*!
 * \brief method to predict tracks from last frame
 * 
 * makes copy of the tracks of lastFrame into current frame and predicts the state
 * deletes tracks, if covariance is too high or if track is out of range
 * 
 * \param lf Last frame tracks are copied from
 * \param cf Current frame tracks are copied into and updated
 */
void predictTracks(FrameSPtr lf, FrameSPtr cf);

/*!
 * \brief method to decide about create/merge/split tracks
 * 
 * \param frame  The frame used for reading and storing
 * \param reg    Registration object for error-calculation
 * \param coutOnly if true, all values are just written to cout for training data collection or debugging
 */
void decideMerge(FrameSPtr frame, TrackRegistration::SPtr reg, bool coutOnly);

/*!
 * \brief method to create/merge/split tracks, delete is handled in prediction
 *
 * \param frame  The frame used for reading and storing
 */
void mergeSplitCreateTracks(FrameSPtr frame, TrackRegistration::SPtr reg);


/*!
 * \brief method to dump the details of all long-term-tracks into a text file
 *
 * \param frame  The frame containing the tracks
 */
void outputLTTTrackDetails(FrameSPtr frame, bool skipWorld = false);

#endif /*ALGOTRACKGENERATION_H_*/
