#include "RobustFit.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>

namespace
{
double Quantile(std::vector<double> values, double fraction)
{
	const std::size_t index = static_cast<std::size_t>((values.size() - 1) * fraction);
	std::nth_element(values.begin(), values.begin() + index, values.end());
	return values[index];
}

void EstimatePixel(Pixel& pixel, double noiseRatio, unsigned seed)
{
	using namespace NormConfig;
	const auto& points = pixel.sample;
	if (points.size() < MIN_PIXEL_EVENTS)
		return;
	std::vector<double> xs, ys;
	for (const auto& p : points)
	{
		xs.push_back(p.x);
		ys.push_back(p.back);
	}
	pixel.span = Quantile(xs, 0.9) - Quantile(xs, 0.1);
	if (pixel.span < MIN_X_SPAN)
	{
		pixel.status = NarrowEnergyRange;
		return;
	}
	pixel.centerBack = Quantile(ys, 0.5);
	const double scaleX = std::max(NOISE_FLOOR, pixel.span);
	double bestSlope = 1.0, bestCenter = Quantile(xs, 0.5);
	double bestScore = std::numeric_limits<double>::infinity();
	std::mt19937 random(seed);
	std::uniform_int_distribution<std::size_t> choose(0, points.size() - 1);
	std::vector<double> residuals(points.size());

	// 最小中位绝对正交残差给出稳健初值，避免尾部事例控制初始直线。
	for (int trial = 0; trial < SEED_TRIALS; ++trial)
	{
		const auto& a = points[choose(random)];
		const auto& b = points[choose(random)];
		const double dy = b.back - a.back;
		if (std::abs(b.x - a.x) < pixel.span * 0.25 || dy == 0)
			continue;
		const double slope = (b.x - a.x) / dy;
		if (!(slope > 0 && std::isfinite(slope)))
			continue;
		const double center = a.x + slope * (pixel.centerBack - a.back);
		const double denominator = std::hypot(1.0, slope * noiseRatio);
		for (std::size_t i = 0; i < points.size(); ++i)
			residuals[i] = std::abs(points[i].x - center - slope * (points[i].back - pixel.centerBack)) / denominator;
		const double score = Quantile(residuals, 0.5);
		if (score < bestScore)
		{
			bestScore = score;
			bestSlope = slope;
			bestCenter = center;
		}
	}
	if (!std::isfinite(bestScore))
	{
		pixel.status = PixelFitFailed;
		return;
	}
	const double scale = std::max(NOISE_FLOOR, 1.4826 * bestScore);
	const double centerOrigin = bestCenter;
	const auto objective = [&](const double* par)
	{
		const double slope = std::exp(par[0]);
		const double center = centerOrigin + scaleX * par[1];
		const double denominator = scale * std::hypot(1.0, slope * noiseRatio);
		double loss = 0.0;
		for (const auto& p : points)
		{
			const double r = (p.x - center - slope * (p.back - pixel.centerBack)) / denominator;
			// soft-L1，r^2 较小时稳定计算，目标在原点具有二次曲率。
			loss += 2.0 * r * (r / (std::hypot(1.0, r) + 1.0));
		}
		return loss;
	};
	ROOT::Math::Functor function(objective, 2);
	ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
	minimizer.SetFunction(function);
	minimizer.SetPrintLevel(-1);
	minimizer.SetTolerance(1e-5);
	minimizer.SetMaxFunctionCalls(3000);
	minimizer.SetVariable(0, "log_slope", std::log(bestSlope), 0.001);
	minimizer.SetVariable(1, "center_shift", 0.0, 0.001);
	if (!minimizer.Minimize() || !minimizer.Hesse() || minimizer.CovMatrixStatus() != 3)
	{
		pixel.status = PixelFitFailed;
		return;
	}
	pixel.slope = std::exp(minimizer.X()[0]);
	pixel.centerX = centerOrigin + scaleX * minimizer.X()[1];
	pixel.offset = pixel.centerX - pixel.slope * pixel.centerBack;
	pixel.varianceSlope = std::pow(pixel.slope, 2) * minimizer.CovMatrix(0, 0);
	pixel.varianceCenter = scaleX * scaleX * minimizer.CovMatrix(1, 1);
	pixel.covariance = pixel.slope * scaleX * minimizer.CovMatrix(0, 1);
	pixel.width = scale * std::hypot(1.0, pixel.slope * noiseRatio);
	int inliers = 0;
	xs.clear();
	for (const auto& p : points)
	{
		if (std::abs(p.x - pixel.centerX - pixel.slope * (p.back - pixel.centerBack)) < 3 * pixel.width)
		{
			++inliers;
			xs.push_back(p.x);
		}
	}
	pixel.inlierFraction = static_cast<double>(inliers) / points.size();
	const double determinant = pixel.varianceSlope * pixel.varianceCenter - pixel.covariance * pixel.covariance;
	pixel.status = PoorPixelFit;
	if (inliers < MIN_PIXEL_EVENTS || pixel.inlierFraction < MIN_INLIER_FRACTION ||
	    !(pixel.slope > 0 && determinant > 0 && std::isfinite(determinant)) ||
	    std::sqrt(pixel.varianceSlope) / pixel.slope > MAX_REL_SLOPE_ERROR)
		return;
	if (Quantile(xs, 0.9) - Quantile(xs, 0.1) < MIN_X_SPAN)
	{
		pixel.status = NarrowEnergyRange;
		return;
	}
	pixel.status = Estimated;
}
}

void EstimatePixels(Pixels& pixels)
{
	int done = 0, accepted = 0;
	for (auto& [channels, pixel] : pixels)
	{
		const bool yh = channels.second >= DSSD_YH_CH_LOW;
		EstimatePixel(pixel, yh ? NormConfig::SIGMA_YH_OVER_X : NormConfig::SIGMA_Y_OVER_X,
		              NormConfig::SEED + channels.first * TOTAL_CH + channels.second);
		accepted += pixel.status == Estimated;
		if (++done % 500 == 0)
			std::cout << "像素拟合 " << done << '/' << pixels.size() << "，有效 " << accepted << std::endl;
	}
}

StripFit FitStrips(const Pixels& pixels)
{
	std::array<std::vector<int>, TOTAL_CH> neighbors;
	std::array<ULong64_t, TOTAL_CH> support{};
	for (const auto& [channels, pixel] : pixels)
	{
		if (pixel.status != Estimated)
			continue;
		const auto [x, back] = channels;
		neighbors[x].push_back(back);
		neighbors[back].push_back(x);
		support[x] += pixel.events;
		support[back] += pixel.events;
	}
	std::array<bool, TOTAL_CH> visited{};
	std::vector<int> active;
	ULong64_t bestSupport = 0;
	const int requested = NormConfig::REFERENCE_X;
	if (requested < -1 || requested > DSSD_X_CH_HIGH - DSSD_X_CH_LOW)
		throw std::runtime_error("REFERENCE_X 超出 X 局部条号范围");
	for (int start = DSSD_X_CH_LOW; start <= DSSD_X_CH_HIGH; ++start)
	{
		if (visited[start] || neighbors[start].empty())
			continue;
		std::vector<int> component{start};
		visited[start] = true;
		ULong64_t componentSupport = 0;
		for (std::size_t i = 0; i < component.size(); ++i)
		{
			const int channel = component[i];
			componentSupport += support[channel];
			for (const int next : neighbors[channel])
			{
				if (!visited[next])
				{
					visited[next] = true;
					component.push_back(next);
				}
			}
		}
		if (requested >= 0)
		{
			if (std::find(component.begin(), component.end(), requested + DSSD_X_CH_LOW) != component.end())
				active = std::move(component);
		}
		else if (component.size() > active.size() ||
		         (component.size() == active.size() && componentSupport > bestSupport))
		{
			bestSupport = componentSupport;
			active = std::move(component);
		}
	}
	if (active.empty())
		throw std::runtime_error("参考条所在分量没有有效像素，请检查统计量、能量跨度和筛选参数");
	StripFit fit;
	fit.strips = active.size();
	fit.reference = requested >= 0 ? requested + DSSD_X_CH_LOW : active.front();
	if (requested < 0)
	{
		for (const int ch : active)
		{
			if (ch <= DSSD_X_CH_HIGH && support[ch] > support[fit.reference])
				fit.reference = ch;
		}
	}
	fit.gain[fit.reference] = 1.0;
	fit.offset[fit.reference] = 0.0;
	std::vector<int> queue{fit.reference};
	for (std::size_t i = 0; i < queue.size(); ++i)
	{
		const int current = queue[i];
		for (const int next : neighbors[current])
		{
			if (std::isfinite(fit.gain[next]))
				continue;
			const bool fromX = current <= DSSD_X_CH_HIGH;
			const auto& pixel = pixels.at(fromX ? std::make_pair(current, next) : std::make_pair(next, current));
			if (fromX)
			{
				fit.gain[next] = fit.gain[current] * pixel.slope;
				fit.offset[next] = fit.offset[current] + fit.gain[current] * pixel.offset;
			}
			else
			{
				fit.gain[next] = fit.gain[current] / pixel.slope;
				fit.offset[next] = fit.offset[current] - fit.gain[next] * pixel.offset;
			}
			queue.push_back(next);
		}
	}
	std::array<int, TOTAL_CH> index;
	index.fill(-1);
	int npar = 0;
	for (const int channel : active)
	{
		if (channel != fit.reference)
		{
			index[channel] = npar;
			npar += 2;
		}
	}
	struct Constraint
	{
		int front, back;
		const Pixel* pixel;
	};
	std::vector<Constraint> constraints;
	for (const auto& [channels, pixel] : pixels)
	{
		if (pixel.status == Estimated && std::isfinite(fit.gain[channels.second]))
			constraints.push_back({index[channels.first], index[channels.second], &pixel});
	}
	const auto evaluate = [&](const double* par, double* gradient)
	{
		if (gradient)
			std::fill_n(gradient, npar, 0.0);
		double loss = 0.0;
		for (const auto& c : constraints)
		{
			const auto& p = *c.pixel;
			const double kx = c.front >= 0 ? std::exp(par[c.front]) : 1.0;
			const double bx = c.front >= 0 ? par[c.front + 1] * NormConfig::OFFSET_SCALE : 0.0;
			const double kb = std::exp(par[c.back]);
			const double bb = par[c.back + 1] * NormConfig::OFFSET_SCALE;
			const double s = kb / kx;
			const double center = (bb - bx) / kx + s * p.centerBack;
			const double ds = s - p.slope, dc = center - p.centerX;
			// 保留斜率与中心幅度的协方差，采用二维 Mahalanobis 残差。
			const double a = p.covariance / p.varianceSlope;
			const double conditional = p.varianceCenter - a * p.covariance;
			const double v = dc - a * ds;
			const double q = ds * ds / p.varianceSlope + v * v / conditional;
			const double root = std::sqrt(1.0 + q / 9.0);
			loss += 2.0 * q / (root + 1.0);
			if (gradient)
			{
				const double gs = 2.0 * (ds / p.varianceSlope - a * v / conditional) / root;
				const double gc = 2.0 * v / conditional / root;
				gradient[c.back] += gs * s + gc * s * p.centerBack;
				gradient[c.back + 1] += gc * NormConfig::OFFSET_SCALE / kx;
				if (c.front >= 0)
				{
					gradient[c.front] -= gs * s + gc * center;
					gradient[c.front + 1] -= gc * NormConfig::OFFSET_SCALE / kx;
				}
			}
		}
		return loss;
	};
	const auto objective = [&](const double* p) { return evaluate(p, nullptr); };
	const auto gradient = [&](const double* p, double* g) { evaluate(p, g); };
	ROOT::Math::GradFunctor function(objective, npar, gradient);
	ROOT::Minuit2::Minuit2Minimizer minimizer(ROOT::Minuit2::kMigrad);
	minimizer.SetFunction(function);
	minimizer.SetPrintLevel(-1);
	minimizer.SetStrategy(1);
	minimizer.SetTolerance(1e-4);
	minimizer.SetMaxFunctionCalls(50000);
	for (const int channel : active)
	{
		const int i = index[channel];
		if (i >= 0)
		{
			minimizer.SetVariable(i, "log_gain_" + std::to_string(channel), std::log(fit.gain[channel]), 0.001);
			minimizer.SetVariable(i + 1, "offset_" + std::to_string(channel),
			                      fit.offset[channel] / NormConfig::OFFSET_SCALE, 0.001);
		}
	}
	if (!minimizer.Minimize() || !std::isfinite(minimizer.MinValue()))
		throw std::runtime_error("全局拟合未收敛，Minuit2 状态 " + std::to_string(minimizer.Status()));
	for (const int channel : active)
	{
		const int i = index[channel];
		if (i >= 0)
		{
			fit.gain[channel] = std::exp(minimizer.X()[i]);
			fit.offset[channel] = minimizer.X()[i + 1] * NormConfig::OFFSET_SCALE;
			if (!(fit.gain[channel] > 0 && std::isfinite(fit.gain[channel]) && std::isfinite(fit.offset[channel])))
				throw std::runtime_error("全局拟合产生非法系数");
		}
	}
	fit.pixels = constraints.size();
	fit.objective = minimizer.MinValue();
	return fit;
}
