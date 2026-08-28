#include "math.h"
#include "common.h"

/*
	   This file contains interpAzOffsets, interpAzSigma, interpRangeOffsetInMeters, interpRangeSigma, which are used to interpolate the offsets fields.
*/

static void computeCoordsForInterp(inputImageStructure *inputImage, Offsets *offsets, double range, double azimuth,
								   double *rangeOff, double *azimuthOff, double *imageLength, double *normAzimuth);

/*
   Convert multi-look range/azimuth pixel coordinates to SLC pixel coordinates.
   range, azimuth : multi-look image coordinates
   *slcRange, *slcAzimuth : SLC pixel coordinates (center of first ML pixel convention)
*/
void computeSLCFromMLCoords(inputImageStructure *inputImage, double range, double azimuth,
                            double *slcRange, double *slcAzimuth)
{
    double rgFirstCenter = (inputImage->nRangeLooks - 1) * 0.5;
    double azFirstCenter = (inputImage->nAzimuthLooks - 1) * 0.5;
    *slcRange   = rgFirstCenter + range   * inputImage->nRangeLooks;
    *slcAzimuth = azFirstCenter + azimuth * inputImage->nAzimuthLooks;
}

#define MAXSIG 0.2
/*
   Interpolate azimuth offsets map for velocity generation and apply baseline/geometry corrections
*/
float interpAzOffset(double range, double azimuth, Offsets *offsets, inputImageStructure *inputImage, double Range, double theta,
					 float azSLPixSize)
{
	float result;
	float zeroOffset;
	double alongTrack;
	double imageLength, normAzimuth;
	double azimuthML, rangeML;
	double azimuthOff, rangeOff;
	/*
	  Compute azimuth in image coord stuff
	*/

	azimuthML = azimuth; /* compute coords will convert input azimuth to range offsets, so save */
	rangeML = range;	 /* compute coords will convert input azimuth to range offsets, so save */

	computeCoordsForInterp(inputImage, offsets, range, azimuth, &rangeOff, &azimuthOff, &imageLength, &normAzimuth);
	alongTrack = normAzimuth * offsets->doffdx;
	/*
		Interpolate
	*/

	result = bilinearInterp((float **)offsets->da, rangeOff, azimuthOff, offsets->nr, offsets->na, -0.9999 * LARGEINT, (float)-LARGEINT);
	if (result < -0.9999 * LARGEINT)
		return -LARGEINT;
	if (offsets->deltaB != DELTABNONE)
	{
		zeroOffset = svAzOffset(inputImage, offsets, rangeML, azimuthML); /* Offset from SV */
	}
	else
		zeroOffset = 0.0;
	
	zeroOffset += (float)offsets->c1 + Range * sin(theta) * offsets->dbcds - Range * cos(theta) * offsets->dbhds;
	
	/* Apply scaling corrections */
	/* if (inputImage->lookDir == LEFT)
		result *= -1.0; */
	result *= azSLPixSize;
	result -= zeroOffset;
	result -= alongTrack;
	//fprintf(stderr," Zero offset %f\n", result);	
/*
if(rangeOff > 0 && rangeOff < offsets->nr && azimuthOff > 0 && azimuthOff < offsets->na)
{
	fprintf(stderr, "interpAzOffset: range %f azimuth %f rangeOff %f azimuthOff %f normAzimuth %f alongTrack %f %i %i\n",
		 range, azimuth, rangeOff, azimuthOff, normAzimuth, alongTrack, offsets->nr, offsets->na); *
	fprintf(stderr, "interpAzOffset: da %f %f %f %f %f %f %f %f\n", result, zeroOffset, alongTrack, offsets->c1, offsets->dbcds, offsets->dbhds, Range, theta	);
	*/
	

	return result;
}

/*
   Interpolate azimuth sigma  map for velocity generation
*/
float interpAzSigma(double range, double azimuth, Offsets *offsets, inputImageStructure *inputImage, double Range, double theta, float azSLPixSize)
{
	float result;
	double alongTrack;
	double imageLength, normAzimuth;
	double sigmaStreaks;
	double rangeOff, azimuthOff;
	/*
	  Compute azimuth in image coord stuff
	*/
	computeCoordsForInterp(inputImage, offsets, range, azimuth, &rangeOff, &azimuthOff, &imageLength, &normAzimuth);
	alongTrack = normAzimuth * offsets->doffdx;
	/*
	 Interpolate
	*/
	result = bilinearInterp((float **)offsets->sa, rangeOff, azimuthOff, offsets->nr, offsets->na, -0.99999 * LARGEINT, (float)-LARGEINT);
	/*if (result < -0.99999*LARGEINT) return -LARGEINT;	*/
	/* Added 10/20/2017 to cap errors at MAXSIG of a pixel - sigma streaks could make it larger */
	if (result > MAXSIG || result < 0)
		result = MAXSIG;
	/* Crudely emprical noise floor for sigmaStreaks  - neglible but can increment if needed */
	sigmaStreaks = max(offsets->sigmaStreaks, 0.001);
	result = sqrt(result * result + sigmaStreaks * sigmaStreaks);
	result *= azSLPixSize;
	return result;
}

/*
   Interpolate range offset  map for velocity generation and apply baseline/geometry corrections
*/
float interpRangeOffsetInMeters(double range, double azimuth, Offsets *offsets, inputImageStructure *inputImage,
						double Range, double thetaD, float rSLPixSize, double theta, double *demError)
{
	float result;
	float zeroOffset;
	double bn, bp, bSq;
	double bnS, bpS;
	double xsq;
	double imageLength, normAzimuth, azimuthML;
	double rangeOff, azimuthOff;

	/*	  Compute azimuth in image coord stuff	*/
	azimuthML = azimuth; /* compute coords will convert input azimuth to range offsets, so save */
	computeCoordsForInterp(inputImage, offsets, range, azimuth, &rangeOff, &azimuthOff, &imageLength, &normAzimuth);
	/*
	   Baseline or deltaBaseline (SV case)
	*/
	xsq = normAzimuth * normAzimuth;
	bn = offsets->bn + offsets->dBn * normAzimuth + offsets->dBnQ * xsq;
	bp = offsets->bp + offsets->dBp * normAzimuth + offsets->dBpQ * xsq;

	/*
	 State vector baseline ?
	*/
	if (offsets->deltaB != DELTABNONE)
	{
		svInterpBnBp(inputImage, offsets, azimuthML, &bnS, &bpS);
		bn = bnS + bn;
		bp = bpS + bp;
	}
	bSq = bn * bn + bp * bp;
	/*
	  Note this can be derived from Eq 7, JGlac 1996, page 566. Solution for quadratic equation
	*/
	zeroOffset = sqrt(pow(Range, 2.0) - 2.0 * Range * (bn * sin(thetaD) + bp * cos(thetaD)) + bSq) - Range + offsets->rConst;
	/* Changed 3/1/16 from 30 to 15 for the nominal DEM error */
	*demError = fabs(bn) * 15.0 / (Range * sin(theta));
	/*
	  Interpolate
	*/
	result = bilinearInterp((float **)offsets->dr, rangeOff, azimuthOff, offsets->nr, offsets->na, -0.9999 * LARGEINT, (float)-LARGEINT);
	if (result < -0.9999 * LARGEINT)
		return -LARGEINT;
	/*
	   apply corrections
	*/
	result *= rSLPixSize; /* This puts offset in meters */
	result -= zeroOffset; /* Substract the offset in meters */
	return result;
}

/*
   Interpolate range sigma  map for velocity generation
*/
float interpRangeSigma(double range, double azimuth, Offsets *offsets, inputImageStructure *inputImage, double Range, double thetaD, float rSLPixSize)
{
	float result;
	double imageLength, normAzimuth;
	double rangeOff, azimuthOff;
	/*
	  Compute azimuth in image coord stuff
	*/
	computeCoordsForInterp(inputImage, offsets, range, azimuth, &rangeOff, &azimuthOff, &imageLength, &normAzimuth);
	/*
	 Interpolate
	*/
	result = bilinearInterp((float **)offsets->sr, rangeOff, azimuthOff, offsets->nr, offsets->na, -0.9999 * LARGEINT, (float)-LARGEINT);
	/*
		If sigma range is set, it represents a minimum error
	*/
	if (result > MAXSIG || result < 0.0)
		result = MAXSIG;
	result = max(offsets->sigmaRange, result);
	result *= rSLPixSize; /* put in units of meters */
	return result;
}

static void computeCoordsForInterp(inputImageStructure *inputImage, Offsets *offsets, double range, double azimuth,
								   double *rangeOff, double *azimuthOff, double *imageLength, double *normAzimuth)
{
	double slcRange, slcAzimuth;
	computeSLCFromMLCoords(inputImage, range, azimuth, &slcRange, &slcAzimuth);
	*imageLength = (double)inputImage->azimuthSize;
	*normAzimuth = (azimuth - 0.5 * (*imageLength)) / (*imageLength);
	*rangeOff   = (slcRange   - offsets->rO) / offsets->deltaR;
	*azimuthOff = (slcAzimuth - offsets->aO) / offsets->deltaA;
	return;
}

/*
  Variance (m^2 of slant range) to add to the range-offset error budget for
  error the .sr band cannot represent.

  interpRangeSigma() returns the .sr band -- a LOCAL neighbourhood scatter that
  Cullst computes after removing a local plane -- so it describes matching noise
  and is structurally blind to long-wavelength error (ionosphere, orbit ramps).
  Measured against Sentinel-1 over stable ground, the offsets formal error came
  out ~6x too small and, worse, SMALLER than the phase formal error, inverting
  the inverse-variance weighting in computeVxy() so crossing offsets dragged the
  combined solution.  See mosaicSource/CLAUDE.md.

  Two ways to supply the missing term, and the difference matters:

  -rSigmaConst X   adds a fixed X metres to EVERY frame.  Raises the offsets
                   budget relative to phase -- which is the miscalibration that
                   corrupts the combined product -- while leaving the relative
                   weighting BETWEEN offset frames untouched.

  -rSigmaResidual  adds each frame's own rparams tie-point fit residual
                   (offsets.sigmaRresidual).  DEFAULT ON.  The exact analogue of
                   the phase path's min(6*PI, vhParam->sigma), and it adapts to
                   the sensor automatically, which a fixed constant cannot.

  How the two terms divide the work is SENSOR-DEPENDENT, and that is the point:

    NISAR   ionosphere inflates the residual to ~0.22 m, ~5x the .sr term, so it
            dominates and the budget becomes essentially per-frame.  That is
            correct: such data really is frame-limited, and noisy frames get
            down-weighted wholesale.
    TSX     X-band, virtually no ionosphere, so .sr captures essentially all the
            noise and should mirror the residual.  Worst case the two
            double-count, i.e. sigma overestimated by ~sqrt(2).
    S1      a mix of the two regimes.

  So the weighting is per-frame where the frame-level error dominates and
  per-pixel where it does not, sliding between them on its own.  A sqrt(2)
  overestimate in the X-band limit is an acceptable price for that.

  Note the side effect, which is expected rather than a defect: .sr varies
  spatially within a frame (measured p90/p10 ~ 7x in sigma, ~48x in weight), and
  adding a term much larger than it flattens that spatial discrimination -- 99%
  of it for NISAR, 76% for an S1-like 0.04 m residual.  Where the frame-level
  error genuinely dominates, flat weighting IS the right answer.

  -rSigmaConst X   adds a fixed X metres instead.  Diagnostic only: it cannot
                   adapt across sensors and flattens the spatial weighting just
                   as thoroughly.  Takes precedence when set.

  The sigma<0 "no solution" sentinel is filtered by the callers before this is
  reached.
*/
double rangeAccuracyVar(Offsets *offsets)
{
	if (rSigmaConst > 0.0)
	{
		return rSigmaConst * rSigmaConst;
	}
	if (rSigmaResidual == TRUE && offsets->sigmaRresidual > 0.0)
	{
		return offsets->sigmaRresidual * offsets->sigmaRresidual;
	}
	return 0.0;
}
