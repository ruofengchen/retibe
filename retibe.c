//
//  retibe.c
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

#include <stdio.h>
#include "retibe.h"
#include <math.h>
#include <assert.h>

// compute a value of a triangle window function of specified size
t_float triangle(t_int t, t_int size)
{
	return 1.f - fabsf(t - (size - 1.f) / 2.f) * 2.f / (size - 1.f);
}

t_float hann(t_int t, t_int size)
{
    return 0.5f * (1.f - cosf(M_2_PI * t / (size - 1.f)));
}

t_float centeredTriangle(t_int t, t_int fadeToZero)
{
	return 1.f - fabsf(t) / fadeToZero;
}

// evaluate on t/p
t_float bell_shape(t_int t, t_int p, t_float tightness)
{
	return expf( -0.5f * powf( tightness * logf( (t_float)t / p ), 2.f ) );
}

t_float bell_shape2(t_int t, t_int p)
{
    return (t_float)t / (p * p) * expf(- (t_float) (t * t) / (2 * p * p));
}

t_float normpdf(t_int t, t_int mu, t_int sigma)
{
    return 0.3989 * expf(-0.5f * (t - mu) * (t - mu) / (sigma * sigma));
}

// autocorrelation function in a smaller range [lower, upper] and return the increased current autocorrelation size
t_int update_acf(t_sample* s, t_int windowSize, t_int acfSize, t_int lower, t_int upper, t_sample* acf)
{
	t_int i;
	for (i=lower; i<=upper; i++) {
		if (i < acfSize) {
			t_float sum = acf[i-MIN_PERIOD] * (acfSize-i-1);
			sum += s[windowSize-1] * s[windowSize-1-i];
			sum /= (acfSize-i);
			acf[i-MIN_PERIOD] = sum;
		}
	}
	return acfSize + 1;
}

t_int compute_tempo_given_old_period(t_sample* acf, t_int lowerPeriod, t_int upperPeriod, t_int oldPeriod, t_float slope)
{
    
    // find maximum value within the region
    t_float maxVal = 0.f;
    t_int maxIdx = 0;
    t_int i;
    for (i=lowerPeriod-MIN_PERIOD; i<=upperPeriod-MIN_PERIOD; i++) {
        t_float weight = 1.f - slope * fabsf(i - oldPeriod); // we want tempo to stay as constant as possible, so we control how much they can deviate
        if (acf[i] * weight > maxVal) {
            maxVal = acf[i] * weight;
            maxIdx = i;
        }
    }
    
    if (maxVal == 0.f) {
        return 43-MIN_PERIOD; // BPM 120 - most probable tempo, right?
    } else {
        return maxIdx;
    }
	
}

t_int compute_tempo_given_old_period_and_prior(t_sample* acf, t_int lowerPeriod, t_int upperPeriod, t_int oldPeriod, t_float slope, t_float* prior)
{
	// find maximum value within the region
	t_float maxVal = 0.f;
	t_int maxIdx = 0;
    t_int i;
	for (i=lowerPeriod-MIN_PERIOD; i<=upperPeriod-MIN_PERIOD; i++) {
		t_float weight = 1.f - slope * fabsf(i - oldPeriod); // we want tempo to stay as constant as possible, so we control how much they can deviate
//        t_float weight = periodTransition[oldPeriod*(MAX_PERIOD-MIN_PERIOD+1)+i];
		if (acf[i] * weight * prior[i] > maxVal) {
			maxVal = acf[i] * weight * prior[i];
			maxIdx = i;
		}
	}
    
	if (maxVal == 0.f) {
		return 43-MIN_PERIOD; // BPM 120 - most probable tempo, right?
	} else {
		return maxIdx;
	}
}

// given localscore and transition profile, compute a cumulative score and backtracked path
void dynamic_programming(t_retibe_tilde* x)
{
	t_int i, j;
	// if tempo changes, recompute cumscore and backlink
	if (x->recomputeCumscore) {
		for (i=0; i<x->windowSize; i++) {
			// find the max transition score
			t_float maxVal = -1;
			t_float maxIdx = 0;
			for (j=x->lowerPeriod-MIN_PERIOD; j<=x->upperPeriod-MIN_PERIOD; j++) {
				if (i-j-MIN_PERIOD >= 0 && x->cumscore[i-j-MIN_PERIOD] * x->transitionProfile[j] > maxVal) {
					maxVal = x->cumscore[i-j-MIN_PERIOD] * x->transitionProfile[j];
					maxIdx = j;
				}
			}
			
			if (maxVal == -1) {
				// initialization
				x->cumscore[i] = x->df[i];
				x->backlink[i] = i;
			}
			else {
				// dynamic programming routine
				x->cumscore[i] = maxVal * x->transitionWeight + x->df[i] * (1-x->transitionWeight);
				x->backlink[i] = i - maxIdx - MIN_PERIOD;
			}
		}
	}
	// if tempo stays the same, simply append cumscore and backlink
	else {
		// shift previously computed cumscore and backlink
		for (i=0; i<x->windowSize-1; i++) {
			x->cumscore[i] = x->cumscore[i+1];
			x->backlink[i] = x->backlink[i+1] - 1;
		}
		
		// it used to be a for loop
		t_int i = x->windowSize-1;
		
		// find the max transition score
		t_float maxVal = -1;
		t_float maxIdx = 0;
		for (j=x->lowerPeriod-MIN_PERIOD; j<=x->upperPeriod-MIN_PERIOD; j++) {
			if (i-j-MIN_PERIOD >= 0 && x->cumscore[i-j-MIN_PERIOD] * x->transitionProfile[j] > maxVal) {
				maxVal = x->cumscore[i-j-MIN_PERIOD] * x->transitionProfile[j];
				maxIdx = j;
			}
		}
		
		if (maxVal == -1) {
			// initialization
			x->cumscore[i] = x->df[i];
			x->backlink[i] = i;
		}
		else {
			// dynamic programming routine
			x->cumscore[i] = maxVal * x->transitionWeight + x->df[i] * (1-x->transitionWeight);
			x->backlink[i] = i - maxIdx - MIN_PERIOD;
		}
	}

	// sum up two consecutive beat cumscores, find the max
	t_int lastBeat = x->windowSize-1;
	t_float maxLastBeatScore = 0.f;
	for (i=1; i<x->bestPeriod+MIN_PERIOD; i++) {
		t_int t = x->windowSize - i;
        t_float peakScoreSum = x->cumscore[t] + x->cumscore[x->backlink[t]];
		if (peakScoreSum > maxLastBeatScore) {
			maxLastBeatScore = peakScoreSum;
			lastBeat = t;
		}
	}
	
	// project to the future
	t_int projectedBeat = lastBeat + x->bestPeriod + MIN_PERIOD - x->windowSize;
	for (i=0; i<MAX_PERIOD-1; i++) {
		x->beatProjection[i] = x->beatProjection[i+1];
	}
	x->beatProjection[MAX_PERIOD-1] = 0;
	if (0 <= projectedBeat && projectedBeat < x->upperPeriod) {
		x->beatProjection[projectedBeat]++;
	}
	
	x->afterFoundBeatCount++;
	if (x->afterFoundBeatCount > x->minBeatIntervalFactor * (x->bestPeriod+MIN_PERIOD)) { // this is to avoid too rapid beats
		// test if the maximum value of beat projection is now (index=offset)
		t_int findBeat = 1;
		for (i=0; i<MAX_PERIOD; i++) {
			if (i != x->offset && x->beatProjection[i] >= x->beatProjection[x->offset]) {
				findBeat = 0;
				break;
			}
		}
		
		if (findBeat) {
			// fire!
			x->afterFoundBeatCount = 0;
			outlet_bang(x->bang);
		}
	}
}

// music/non-music classification
t_int is_music(t_retibe_tilde* x)
{
    t_int i;
    // care about silent head - so it's good at first launch
    t_int zeroCount = 0;
    for (i=0; i<x->frameSize; i++) {
        if (x->frameSignal[i] == 0) {
            zeroCount ++;
        }
    }
    
	// zero-crossing rate
	t_int zeroCrossCount = 0;
	for (i=zeroCount; i<x->frameSize-1; i++) {
		if (x->frameSignal[i] * x->frameSignal[i+1] < 0) {
			zeroCrossCount++;
		}
	}
	t_float zcr = (t_float) zeroCrossCount / (x->frameSize-zeroCount);
	
	// mean-square power
	t_float power = 0.f;
	for (i=zeroCount; i<x->frameSize; i++) {
		power += x->frameSignal[i] * x->frameSignal[i];
	}
	power /= (x->frameSize-zeroCount);
	
	// put them into statistic arrays
	for (i=0; i<x->isMusicSize-1; i++) {
		x->zcr[i] = x->zcr[i+1];
		x->power[i] = x->power[i+1];
	}
	x->zcr[x->isMusicSize-1] = zcr;
	x->power[x->isMusicSize-1] = power;
    
    // care about silent head
    t_int silenceCount = 0;
    for (i=0; i<x->isMusicSize; i++) {
        if (x->zcr[i] == 0) {
            silenceCount++;
        }
    }
	
	// weight more on recent values
	t_float zcrMean = 0.f;
	t_float powerMean = 0.f;
	for (i=silenceCount; i<x->isMusicSize; i++) {
		zcrMean += (i-silenceCount+1) * x->zcr[i];
		powerMean += (i-silenceCount+1) * x->power[i];
	}
	
    t_int modifiedSize = x->isMusicSize - silenceCount;
	zcrMean /= (modifiedSize * (1+modifiedSize) * 0.5);
	powerMean /= (modifiedSize * (1+modifiedSize) * 0.5);
	
	if (zcrMean > x->zcrTH && powerMean > x->powerTH) return 1;
	else return 0;
}

void reset(t_retibe_tilde* x)
{
    t_int i;
	printf("reset!\n");
    
	outlet_bang(x->resetBang);
	
	for (i=0; i<x->frameSize; i++) {
		x->frameSignal[i] = 0.f;
	}
	for (i=0; i<x->hopSize; i++) {
		x->hopSignal[i] = 0.f;
	}
	for (i=0; i<x->frameSize; i++) {
		x->spectrum[i] = 0.f;
	}
	for (i=0; i<x->windowSize; i++) {
		x->df[i] = 0.f;
	}
	for (i=0; i<x->differenceSize * x->featureSize; i++) {
		x->subbandFeatures[i] = 0.f;
	}
	for (i=0; i<MAX_PERIOD - MIN_PERIOD + 1; i++) {
		x->acf[i] = 0.f;
	}
	for (i=0; i<MAX_PERIOD - MIN_PERIOD + 1; i++) {
		x->transitionProfile[i] = 0.f;
	}
	for (i=0; i<x->windowSize; i++) {
		x->cumscore[i] = 0.f;
	}
	for (i=0; i<x->windowSize; i++) {
		x->backlink[i] = i;
	}
	for (i=0; i<MAX_PERIOD; i++) {
		x->beatProjection[i] = 0;
	}
	
	x->sampleCount = 0;
	x->afterFoundBeatCount = 0;
	x->recomputeCumscore = 1;
	x->acfSize = 0;
    x->bestPeriod = 43-MIN_PERIOD;
}

// analyze spectrum and do some magic
void beat_tracking(t_retibe_tilde* x)
{
    t_int i, j, k;
    
	// fft (don't touch spectrum[0]!)
	for (i=0; i<x->frameSize; i++)
	{
		x->spectrum[i] = x->frameSignal[i];
	}
	mayer_realfft(x->frameSize, x->spectrum);
	for (i=1; i<x->spectrumSize; i++)
	{
		x->spectrum[i] = (x->spectrum[i] * x->spectrum[i]) + (x->spectrum[x->frameSize-i] * x->spectrum[x->frameSize-i]); // magnitude
		x->spectrum[i] = log10f(x->spectrum[i] + 1); // convert to log spectrum
//        x->spectrum[i] = x->spectrum[i];
	}
	
	// compute subband feature by filterbank multiplication
	for (i=0; i<x->differenceSize-1; i++) {
		for (j=0; j<x->featureSize; j++) {
			x->subbandFeatures[i*x->featureSize+j] = x->subbandFeatures[(i+1)*x->featureSize+j];
		}
	}
	for (j=0; j<x->featureSize; j++) { // can probably save some computation here
		t_float sum = 0.f;
		for (k=0; k<x->spectrumSize; k++) {
			sum += x->spectrum[k] * x->filterbanks[j*x->spectrumSize+k];
		}
		x->subbandFeatures[(x->differenceSize-1)*x->featureSize+j] = sum;
	}
	
	// compute detection function using difference between current energy and previous mean energy, and finally HWR it
	for (i=0; i<x->windowSize-1; i++) {
		x->df[i] = x->df[i+1];
	}
    
    // for each subband, compute advanced distance and HWR it
    for (i=0; i<x->featureSize; i++) {
        // find the first non-zero subbandFeature along time
        t_int numNonZero = 0;
        t_float mean = 0.f;
        for (j=x->differenceSize-2; j>=0; j--) {
            if (x->subbandFeatures[j*x->featureSize+i] != 0) {
                mean += x->subbandFeatures[j*x->featureSize+i];
                numNonZero++;
            }
            else {
                break;
            }
        }
        if (numNonZero == 0) {
            x->subbandDiff[i] = 0.f;
        }
        else {
            mean /= numNonZero;
            x->subbandDiff[i] = x->subbandFeatures[(x->differenceSize-1)*x->featureSize+i] - mean;
            x->subbandDiff[i] = x->subbandDiff[i] > 0.f? x->subbandDiff[i] : 0.f;
        }
    }
    
    t_float sum = 0.f;
    for (i=0; i<x->featureSize; i++) {
        sum += x->subbandDiff[i];
    }
	x->df[x->windowSize-1] = sum;
    
	// compute autocorrelation function
	x->acfSize = update_acf(x->df, x->windowSize, x->acfSize, MIN_PERIOD, MAX_PERIOD, x->acf);

	// compute period (1/tempo) from detection function. period + MIN_PERIOD = real period
	t_int oldBestPeriod = x->bestPeriod;
    if (x->useWindowCenteredAt120BPM) {
        x->bestPeriod = compute_tempo_given_old_period_and_prior(x->acf, x->lowerPeriod, x->upperPeriod, oldBestPeriod, x->periodDeviationSlope, x->tempoPeriodPrior);
    }
    else {
        x->bestPeriod = compute_tempo_given_old_period(x->acf, x->lowerPeriod, x->upperPeriod, oldBestPeriod, x->periodDeviationSlope);
    }
	
	x->recomputeCumscore = oldBestPeriod!=x->bestPeriod?1:0;

	// compute transition profile
    if (x->recomputeCumscore) {
        for (i=x->lowerPeriod-MIN_PERIOD; i<=x->upperPeriod-MIN_PERIOD; i++) {
            x->transitionProfile[i] = bell_shape(i+MIN_PERIOD, x->bestPeriod+MIN_PERIOD, x->tightness);
        }
    }
	
	// find beats in df with the contraint of transition profile
	dynamic_programming(x);
	
	outlet_float(x->debug, x->df[x->windowSize-1]);
}

// prepare for frame level analysis
t_int* retibe_tilde_perform(t_int* w)
{
    t_int i;

	// casting from signal array to appropriate types
	t_retibe_tilde* x = (t_retibe_tilde*) w[1];
	t_sample* in = (t_sample *) w[2];
	t_int bufferSize = (t_int) w[3]; // this is supposed to be 64
	
	// test if the input buffer is all zeros
	t_int nonzeroCount = 0;
	for (i=0; i<bufferSize; i++) {
        if (in[i] != 0) {
            nonzeroCount++;
        }
    }
	if (!nonzeroCount) {
        
		x->resetFlag1 = x->resetFlag1==0?1:2;
		if (x->resetFlag1==1) {
			reset(x);
		}
		return w+4;
	}
	x->resetFlag1 = 0;
	
	// not all zeros, good, continue
	if (x->sampleCount < x->hopSize) {
		
		// append to the tail
		for (i=0; i<bufferSize; i++,x->sampleCount++) {
			x->hopSignal[x->sampleCount] = in[i];
		}
	}
	
    // the following should be performed right after the buffer is fully filled
	if (x->sampleCount==x->hopSize) {
        
		// reset the counter
		x->sampleCount = 0;		
		
		// shift frameSignal and fill hopSignal into it
		for (i=0; i<x->frameSize-x->hopSize; i++) {
			x->frameSignal[i] = x->frameSignal[i+x->hopSize];
		}
		for (i=0; i<x->hopSize; i++) {
			x->frameSignal[i+x->frameSize-x->hopSize] = x->hopSignal[i];
		}
		
		// if not noise, track beat
		if (!x->classifyMusicFlag || is_music(x)) {
			beat_tracking(x);
		}
		else {
            // printf("going to noise\n");
			x->resetFlag2 = x->resetFlag2==0?1:2;
			if (x->resetFlag2==1) {
				reset(x);
			}
			return w+4;
		}
		x->resetFlag2 = 0;

	}
	return w+4;
}

void retibe_tilde_dsp(t_retibe_tilde *x, t_signal **sp)
{
	dsp_add(retibe_tilde_perform, 3, x,
			sp[0]->s_vec, sp[0]->s_n);
}

void* retibe_tilde_new()
{
    t_int i, j;
	t_retibe_tilde *x = (t_retibe_tilde *)pd_new(retibe_tilde_class);
	
	// set up inlets and outlets
	x->bang = outlet_new(&x->x_obj, &s_bang);
	x->resetBang = outlet_new(&x->x_obj, &s_bang);
	x->debug = outlet_new(&x->x_obj, &s_float);
	
	// default hopSize and frameSize
	x->hopSize = 512;
	x->frameSize = 4096;
	x->spectrumSize = x->frameSize / 2;
	x->windowSize = 512; // 512*512/44100~6 second
	x->sampleCount = 0;
	x->resetFlag1 = 0;
	x->resetFlag2 = 0;
	
	x->zcrTH = 0.01f;
	x->powerTH = 0.0001f;
	x->isMusicSize = 64; // 64*512/44100~0.7 second
	x->frameSignal = (t_sample *)t_getbytes(x->frameSize * sizeof(t_sample));
	x->hopSignal = (t_sample *)t_getbytes(x->hopSize * sizeof(t_sample));
	x->spectrum = (t_sample *)t_getbytes(x->frameSize * sizeof(t_sample));
    x->subbandDiff = (t_sample *)t_getbytes(x->featureSize * sizeof(t_sample));
	x->df = (t_sample *)t_getbytes(x->windowSize * sizeof(t_sample));
	x->zcr = (t_sample *)t_getbytes(x->isMusicSize * sizeof(t_sample));
	x->power = (t_sample *)t_getbytes(x->isMusicSize * sizeof(t_sample));
	for (i=0; i<x->isMusicSize; i++) {
		x->zcr[i] = 0.f;
	}
	for (i=0; i<x->isMusicSize; i++) {
		x->power[i] = 0.f;
	}
	
	// filterbank multiplication computation
	x->featureSize = 4;
	x->filterbanks = t_getbytes(x->spectrumSize * x->featureSize * sizeof(t_sample)); // it's a matrix stored in an array	
	t_float freq_start = 20.f;
	t_float freq_end = 14000.f;
	t_float exponent = powf(freq_end / freq_start, 1.f / (x->featureSize + 1));
	for (i=0; i<x->featureSize; i++) {
		t_int band_start = ceilf(freq_start * powf(exponent, i) * 2 / FS * x->spectrumSize); // spectrum[0] is useless DC, we won't count it into df
		t_int band_end = ceilf(freq_start * powf(exponent, i + 2) * 2 / FS * x->spectrumSize);
		
		// compute filterbanks
		for (j=0; j<x->spectrumSize; j++) {
			if (band_start <= j && j <= band_end) {
				x->filterbanks[i * x->spectrumSize + j] = hann(j-band_start+1, band_end-band_start+3) / (band_end-band_start+1);
			}
		}
	}
	
	// subband feature preparation
	x->differenceSize = 10;
	assert(x->differenceSize>1);
	x->subbandFeatures = (t_sample *)t_getbytes(x->differenceSize * x->featureSize * sizeof(t_sample)); // we need the previous few features to compute difference
	
	// dynamic programming preparation
	x->lowerPeriod = 15; // fastest: 30->172 bpm 20->258 bpm
	x->upperPeriod = 119; // slowest: 60->86 bpm 80->64 bpm
    assert(x->lowerPeriod >= MIN_PERIOD && x->upperPeriod <= MAX_PERIOD);
	x->periodDeviationSlope = 0.01f; // 0.08 is too much for some songs
	x->tightness = 5.f; // the width of transition profile
	x->transitionWeight = 0.9f; // how much we trust the tempo versus how much we trust local df
	x->acf = (t_sample *)t_getbytes((MAX_PERIOD - MIN_PERIOD + 1) * sizeof(t_sample));
	x->transitionProfile = (t_sample *)t_getbytes((MAX_PERIOD - MIN_PERIOD + 1) * sizeof(t_sample));
	x->cumscore = (t_sample *)t_getbytes(x->windowSize * sizeof(t_sample));
	x->backlink = (t_int *)t_getbytes(x->windowSize * sizeof(t_int));
	x->beatProjection = (t_int *)t_getbytes(MAX_PERIOD * sizeof(t_int));
	x->minBeatIntervalFactor = 0.3f; // it times best period to make sense
	x->offset = 3;
    x->classifyMusicFlag = 1;
    
    x->useWindowCenteredAt120BPM = 1;
    x->tempoPeriodPrior = (t_sample *)t_getbytes((MAX_PERIOD - MIN_PERIOD + 1) * sizeof(t_sample));
    for (i = 0; i < MAX_PERIOD - MIN_PERIOD + 1; i++) {
        x->tempoPeriodPrior[i] = bell_shape2(i+MIN_PERIOD, 43);
    }
    t_float sum = 0.f;
    for (i = 0; i < MAX_PERIOD - MIN_PERIOD + 1; i++) {
        sum += x->tempoPeriodPrior[i];
    }
    for (i = 0; i < MAX_PERIOD - MIN_PERIOD + 1; i++) {
        x->tempoPeriodPrior[i] /= sum;
    }
    
	// initialize everything to default
	reset(x);
	
	return (void *)x;
}

void retibe_tilde_free(t_retibe_tilde* x)
{
	t_freebytes(x->frameSignal, x->frameSize * sizeof(t_sample));
	t_freebytes(x->hopSignal, x->hopSize * sizeof(t_sample));
	t_freebytes(x->spectrum, x->frameSize * sizeof(t_sample));
    t_freebytes(x->subbandDiff, x->featureSize * sizeof(t_sample));
	t_freebytes(x->df, x->windowSize * sizeof(t_sample));
	t_freebytes(x->zcr, x->isMusicSize * sizeof(t_sample));
	t_freebytes(x->power, x->isMusicSize * sizeof(t_sample));
	t_freebytes(x->filterbanks, x->spectrumSize * x->featureSize * sizeof(t_sample));
	t_freebytes(x->subbandFeatures, x->differenceSize * x->featureSize * sizeof(t_sample));
	t_freebytes(x->acf, (MAX_PERIOD - MIN_PERIOD + 1) * sizeof(t_sample));
	t_freebytes(x->transitionProfile, (MAX_PERIOD - MIN_PERIOD + 1) * sizeof(t_sample));
	t_freebytes(x->cumscore, x->windowSize * sizeof(t_sample));
	t_freebytes(x->backlink, x->windowSize * sizeof(t_int));
	t_freebytes(x->beatProjection, MAX_PERIOD * sizeof(t_int));
	return;
}

// message functions

void retibe_tilde_set_zcrTH(t_retibe_tilde* x, t_float number) {
	if (number <=0) return;
	x->zcrTH = number;
}

void retibe_tilde_set_powerTH(t_retibe_tilde* x, t_float number) {
	if (number <=0) return;
	x->powerTH = number;
}

void retibe_tilde_set_offset(t_retibe_tilde* x, t_float number) {
	if (number <0) return;
	x->offset = (t_int)number;
}

void retibe_tilde_set_minBeatIntervalFactor(t_retibe_tilde* x, t_float number) {
	if (number <=0) return;
	x->minBeatIntervalFactor = number;
}

void retibe_tilde_set_periodDeviationSlope(t_retibe_tilde* x, t_float number) {
	if (number <0) return;
	x->periodDeviationSlope = number;
}

void retibe_tilde_set_classifyMusicFlag(t_retibe_tilde* x, t_float number) {
	x->classifyMusicFlag = number!=0?1:0;
}

void retibe_tilde_reset(t_retibe_tilde* x) {
	reset(x);
}

void retibe_tilde_set_lowerPeriod(t_retibe_tilde* x, t_float number) {
    x->lowerPeriod = number > MIN_PERIOD ? number : MIN_PERIOD;
}

void retibe_tilde_set_upperPeriod(t_retibe_tilde* x, t_float number) {
    x->upperPeriod = number < MAX_PERIOD ? number : MAX_PERIOD;
}

void retibe_tilde_set_useWindowCenteredAt120BPM(t_retibe_tilde* x, t_float number) {
	x->useWindowCenteredAt120BPM = number!=0?1:0;
}

void retibe_tilde_setup(void) {
//	post("retibe is loaded.");
	retibe_tilde_class = class_new(gensym("retibe~"), (t_newmethod)retibe_tilde_new, (t_method)retibe_tilde_free, sizeof(t_retibe_tilde), CLASS_DEFAULT, 			A_DEFFLOAT, 0);
	class_addmethod(retibe_tilde_class, (t_method)retibe_tilde_dsp, gensym("dsp"), 0);
	CLASS_MAINSIGNALIN(retibe_tilde_class, t_retibe_tilde, x_f);
	
	class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_zcrTH, gensym("zcr"), A_FLOAT, 0 );	
	class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_powerTH, gensym("power"), A_FLOAT, 0 );	
	class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_offset, gensym("offset"), A_FLOAT, 0 );	
	class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_minBeatIntervalFactor, gensym("intervalfactor"), A_FLOAT, 0 );	
	class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_periodDeviationSlope, gensym("periodstability"), A_FLOAT, 0 );
    class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_reset, gensym("reset"), 0 );    
    class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_classifyMusicFlag, gensym("classifymusic"), A_FLOAT, 0 );    
    class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_lowerPeriod, gensym("lowerperiod"), A_FLOAT, 0 );
    class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_upperPeriod, gensym("upperperiod"), A_FLOAT, 0 );
    class_addmethod( retibe_tilde_class, (t_method)retibe_tilde_set_useWindowCenteredAt120BPM, gensym("useprior"), A_FLOAT, 0 );
}