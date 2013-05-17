//
//  retibe.h
//
/*
 Copyright (c) 2013, Ruofeng Chen
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 - Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 - The names of its contributors may not be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef retibe_retibe_h
#define retibe_retibe_h

#include "m_pd.h" // 

#define FS 44100.f
#define MIN_PERIOD 15
#define MAX_PERIOD 120

static t_class *retibe_tilde_class;
typedef struct _retibe_tilde {
	t_object  x_obj;
	t_sample x_f;
	
	// useful variables
	t_int hopSize;
	t_int frameSize;
	t_int spectrumSize;
	t_int windowSize; // used to perform dynamic programming
	t_int featureSize; // number of subbands
	t_int differenceSize; // used to compute advanced difference
	t_int acfSize; // used in appending autocorrelation function computation
	t_int lowerPeriod;
	t_int upperPeriod;
	t_int bestPeriod; // CAUTION: index+lowerPeriod=actualPeriod!
	t_float tightness; // the shape of transition profile
	t_float transitionWeight; // how much we emphasize transition
	t_sample* frameSignal;
	t_sample* hopSignal; // the buffer to be appended to frameSignal
	t_sample* spectrum;
    t_sample* subbandDiff;
	t_sample* df; // 1-D score for potential onsets
	t_sample* filterbanks; // a precomputed filter to accelerate the computation of subbandFeatures
	t_sample* subbandFeatures; // a reduced spectrum on octave bands
	t_sample* acf; // a partial autocorrelation score in lag [lowerPeriod, upperPeriod] CAUTION: index+lowerPeriod=actualPeriod!
	t_sample* transitionProfile; // for dynamic programming
	t_sample* cumscore; // dynamic programming cumulative score
	t_int* backlink; // dynamic programming backtracking path
	t_int* beatProjection; // beat score for the future
	t_float periodDeviationSlope;
	t_float minBeatIntervalFactor;
	t_int offset; // find beat how many frames in advance
	t_int isMusicSize; // the window size to classify music/non-music
	t_float* zcr; // used to classify music/non-music
	t_float* power; // used to classify music/non-music
	t_float zcrTH;
	t_float powerTH;
    t_float* tempoPeriodPrior;
    t_int useWindowCenteredAt120BPM;
	
	// outlets
	t_outlet* bang;
	t_outlet* resetBang;
	t_outlet* debug;
	
	// internal states
	t_int sampleCount; // we need this because the frame size in pd is 64 samples but we want more
	t_int recomputeCumscore; // if tempo changes, recompute cumscore
	t_int afterFoundBeatCount; // we don't want beats to be found too often
	t_int resetFlag1;
	t_int resetFlag2;
    t_int classifyMusicFlag;
    t_int debugCount;
    t_int needToWriteToFileFlag;
} t_retibe_tilde;



#endif
