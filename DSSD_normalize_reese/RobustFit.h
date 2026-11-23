#ifndef DSSD_ROBUST_FIT_H
#define DSSD_ROBUST_FIT_H

#include "NormConfig.h"
#include <array>
#include <limits>
#include <map>
#include <vector>

constexpr double kInvalid = std::numeric_limits<double>::quiet_NaN();
struct PairSample
{
	double x, back;
};
enum PixelStatus
{
	TooFewEvents, NarrowEnergyRange, PixelFitFailed, PoorPixelFit, Estimated
};
struct Pixel
{
	std::vector<PairSample> sample;
	ULong64_t events = 0;
	PixelStatus status = TooFewEvents;
	double slope = kInvalid, offset = kInvalid;
	double centerBack = kInvalid, centerX = kInvalid;
	double varianceSlope = kInvalid, varianceCenter = kInvalid, covariance = kInvalid;
	double width = kInvalid, inlierFraction = 0.0, span = 0.0;
};
using Pixels = std::map<std::pair<int, int>, Pixel>;

struct StripFit
{
	std::array<double, TOTAL_CH> gain, offset;
	int reference = -1, strips = 0, pixels = 0;
	double objective = kInvalid;
	StripFit() { gain.fill(kInvalid); offset.fill(kInvalid); }
};

void EstimatePixels(Pixels& pixels);
StripFit FitStrips(const Pixels& pixels);

#endif
