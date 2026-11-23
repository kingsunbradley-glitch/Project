#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH2D.h>
#include <TParameter.h>
#include <TString.h>
#include <TObjString.h>
#include <memory>
#include <random>
#include <sstream>
#include "RobustFit.h"
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

namespace
{
namespace fs = std::filesystem;

struct Readout
{
	const char* branch;
	int first;
	int last;
	double minAmplitude;
	double saturation;
};

// 道址、饱和值、符合时间窗来自共享 setup.h；幅度窗属于本程序的选择配置。
constexpr Readout kX{"DSSDX", DSSD_X_CH_LOW, DSSD_X_CH_HIGH, NormConfig::MIN_X, E_SATURATION_DSSDX};
constexpr Readout kY{"DSSDY", DSSD_Y_CH_LOW, DSSD_Y_CH_HIGH, NormConfig::MIN_Y, E_SATURATION_DSSDY};
constexpr Readout kYH{"DSSDYH", DSSD_YH_CH_LOW, DSSD_YH_CH_HIGH, NormConfig::MIN_YH, E_SATURATION_DSSDYH};
struct Hit
{
	int channel;
	double amplitude;
	ULong64_t time;
};

class MapReadout
{
public:
	MapReadout(TTreeReader& reader, Readout layout)
		: layout_(layout),
		  multiplicity_(reader, Form("%s_mul", layout.branch)),
		  channels_(reader, Form("%s_Ch", layout.branch)),
		  amplitudes_(reader, Form("%s_E", layout.branch)),
		  times_(reader, Form("%s_Ts", layout.branch))
	{
		if (NormConfig::REJECT_PILEUP && reader.GetTree()->GetBranch(Form("%s_PU", layout.branch)))
			pileup_ = std::make_unique<TTreeReaderArray<UShort_t>>(reader, Form("%s_PU", layout.branch));
	}

	bool SingleHit(Hit& hit)
	{
		if (*multiplicity_ != 1)
			return false;
		if (channels_.GetSize() != 1 || amplitudes_.GetSize() != 1 || times_.GetSize() != 1 || (pileup_ && pileup_->GetSize() != 1))
			throw std::runtime_error(std::string("分支数组长度与多重性不符 ") + layout_.branch);

		const double ch = channels_[0];
		const double amplitude = amplitudes_[0];
		if (!(ch >= 0 && ch <= layout_.last - layout_.first && ch == std::floor(ch)
		      && amplitude > layout_.minAmplitude && amplitude < NormConfig::MAX_AMPLITUDE
		      && amplitude != layout_.saturation) || (pileup_ && (*pileup_)[0] != 0))
			return false;
		hit = {layout_.first + static_cast<int>(ch), amplitude, times_[0]};
		return true;
	}

private:
	Readout layout_;
	TTreeReaderValue<UShort_t> multiplicity_;
	TTreeReaderArray<Double_t> channels_;
	TTreeReaderArray<Double_t> amplitudes_;
	TTreeReaderArray<ULong64_t> times_;
	std::unique_ptr<TTreeReaderArray<UShort_t>> pileup_;
};

bool Coincident(const Hit& x, const Hit& back)
{
	return std::max(x.time, back.time) - std::min(x.time, back.time) < COIN_SI_WINDOW;
}

template <typename Consume>
Long64_t ReadPairs(const std::vector<fs::path>& inputs, Long64_t maxEntries, Consume consume)
{
	Long64_t entries = 0;
	for (const auto& path : inputs)
	{
		TFile input(path.c_str(), "READ");
		auto* tree = input.Get<TTree>("tr_map");
		if (input.IsZombie() || !tree)
			throw std::runtime_error("找不到 tr_map " + path.string());
		int yhBranches = 0;
		for (const auto* suffix : {"mul", "Ch", "E", "Ts"})
		{
			for (const auto* side : {"DSSDX", "DSSDY"})
				if (!tree->GetBranch(Form("%s_%s", side, suffix)))
					throw std::runtime_error("缺少分支 " + std::string(side) + "_" + suffix + " in " + path.string());
			yhBranches += tree->GetBranch(Form("DSSDYH_%s", suffix)) != nullptr;
		}
		if (yhBranches != 0 && yhBranches != 4)
			throw std::runtime_error("DSSDYH 分支不完整 " + path.string());
		for (const auto layout : {kX, kY, kYH})
		{
			if (NormConfig::REJECT_PILEUP && (layout.first != kYH.first || yhBranches) &&
			    !tree->GetBranch(Form("%s_PU", layout.branch)))
				std::cout << path.filename() << " 中 " << layout.branch
				          << "_PU 缺省，此组使用单条、幅度和时间选择" << std::endl;
		}
		if (tree->GetEntries() == 0)
			continue;
		TTreeReader reader(tree);
		MapReadout front(reader, kX), back(reader, kY);
		std::unique_ptr<MapReadout> parallel;
		if (yhBranches)
			parallel = std::make_unique<MapReadout>(reader, kYH);
		while (reader.Next())
		{
			if (maxEntries && entries >= maxEntries)
				return entries;
			++entries;
			Hit x{}, y{}, yh{};
			if (!front.SingleHit(x))
				continue;
			const bool hasY = back.SingleHit(y) && Coincident(x, y);
			const bool hasYH = parallel && parallel->SingleHit(yh) && Coincident(x, yh);
			// Y 和 YH 是同一物理条的并行读出；两者都符合时要求物理条号一致。
			if (hasY && hasYH && y.channel - kY.first != yh.channel - kYH.first)
				continue;
			if (hasY)
				consume(x, y);
			if (hasYH)
				consume(x, yh);
		}
		if (reader.GetEntryStatus() != TTreeReader::kEntryBeyondEnd)
			throw std::runtime_error("读取 tr_map 失败 " + path.string());
	}
	return entries;
}

struct PairPlots
{
	TH2D raw, normalized, differenceVsEnergy;
	TH1D spectrumX, spectrumBack;
	TH1D rawResidual, normalizedResidual;

	explicit PairPlots(const std::string& pair)
		: raw((pair + "_before").c_str(), ";X amplitude;Back amplitude",
		      400, 0, NormConfig::MAX_AMPLITUDE, 400, 0, NormConfig::MAX_AMPLITUDE),
		  normalized((pair + "_after").c_str(), ";X normalized;Back normalized",
		             400, 0, NormConfig::MAX_AMPLITUDE, 400, 0, NormConfig::MAX_AMPLITUDE),
		  differenceVsEnergy((pair + "_difference_vs_energy").c_str(), ";X normalized;X - Back normalized",
		                     400, 0, NormConfig::MAX_AMPLITUDE, 1000, -1000, 1000),
		  spectrumX((pair + "_spectrum_x").c_str(), ";X normalized;Pairs", 16050, 0, NormConfig::MAX_AMPLITUDE),
		  spectrumBack((pair + "_spectrum_back").c_str(), ";Back normalized;Pairs", 16050, 0, NormConfig::MAX_AMPLITUDE),
		  rawResidual((pair + "_residual_before").c_str(), ";2(X-Back)/(X+Back);Counts", 2000, -0.1, 0.1),
		  normalizedResidual((pair + "_residual_after").c_str(), ";2(X-Back)/(X+Back);Counts", 2000, -0.1, 0.1)
	{
	}

	void Fill(double x, double back, double sx, double sb, double bx, double bb)
	{
		raw.Fill(x, back);
		rawResidual.Fill(2 * (x - back) / (x + back));
		x = sx * x + bx;
		back = sb * back + bb;
		normalized.Fill(x, back);
		if (x + back != 0) normalizedResidual.Fill(2 * (x - back) / (x + back));
		differenceVsEnergy.Fill(x, x - back);
		spectrumX.Fill(x);
		spectrumBack.Fill(back);
	}

	void Write()
	{
		raw.Write();
		normalized.Write();
		rawResidual.Write();
		normalizedResidual.Write();
		differenceVsEnergy.Write(); spectrumX.Write(); spectrumBack.Write();
	}
};

void WriteOutput(const std::vector<fs::path>& inputs, const Pixels& pixels, const StripFit& fit,
                 const fs::path& rootPath, const fs::path& datPath, Long64_t maxEntries, Long64_t entries, const std::string& fitError)
{
	TH1::AddDirectory(false);
	PairPlots xy("XY"), xyh("XYH");
	if (fit.strips) ReadPairs(inputs, maxEntries, [&](const Hit& x, const Hit& back)
	{
		const double sx = fit.gain[x.channel], sb = fit.gain[back.channel];
		if (std::isfinite(sx) && std::isfinite(sb))
			(back.channel <= kY.last ? xy : xyh).Fill(x.amplitude, back.amplitude, sx, sb, fit.offset[x.channel], fit.offset[back.channel]);
	});

	TFile output(rootPath.c_str(), "NEW");
	if (output.IsZombie())
		throw std::runtime_error("无法创建诊断 ROOT 文件 " + rootPath.string());
	std::ofstream table;
	if (fit.strips) table.open(datPath);
	table << std::setprecision(12);

	int xChannel, backChannel, status;
	ULong64_t events, sampleSize;
	double slope, offset, varianceSlope, varianceCenter, covariance, centerBack, centerX, fraction, width, span;
	bool used;
	TTree pixelTree("pixel_fit", "Robust orthogonal pixel fits; local curvature covariance");
	pixelTree.Branch("x", &xChannel);
	pixelTree.Branch("back", &backChannel);
	pixelTree.Branch("events", &events);
	pixelTree.Branch("slope", &slope);
	pixelTree.Branch("offset", &offset);
	pixelTree.Branch("sample_size", &sampleSize);
	pixelTree.Branch("center_back", &centerBack); pixelTree.Branch("center_x", &centerX);
	pixelTree.Branch("var_slope", &varianceSlope); pixelTree.Branch("var_center", &varianceCenter);
	pixelTree.Branch("cov_slope_center", &covariance);
	pixelTree.Branch("inlier_fraction", &fraction); pixelTree.Branch("width", &width); pixelTree.Branch("span_x", &span);
	pixelTree.Branch("status", &status);
	pixelTree.Branch("used", &used);
	for (const auto& [channels, pixel] : pixels)
	{
		xChannel = channels.first;
		backChannel = channels.second;
		events = pixel.events;
		slope = pixel.slope;
		offset = pixel.offset; sampleSize = pixel.sample.size();
		centerBack = pixel.centerBack; centerX = pixel.centerX;
		varianceSlope = pixel.varianceSlope; varianceCenter = pixel.varianceCenter; covariance = pixel.covariance;
		fraction = pixel.inlierFraction; width = pixel.width; span = pixel.span;
		status = pixel.status;
		used = pixel.status == Estimated && std::isfinite(fit.gain[backChannel]);
		pixelTree.Fill();
	}

	int channel;
	double gain;
	bool reference;
	TTree stripTree("strip_fit", "Relative gains on a common scale");
	stripTree.Branch("channel", &channel);
	stripTree.Branch("gain", &gain);
	stripTree.Branch("offset", &offset);
	stripTree.Branch("status", &status);
	stripTree.Branch("used", &used);
	stripTree.Branch("reference", &reference);
	for (const auto layout : {kX, kY, kYH})
	{
		for (channel = layout.first; channel <= layout.last; ++channel)
		{
			gain = fit.gain[channel];
			offset = fit.offset[channel];
			used = std::isfinite(gain);
			status = used ? 0 : 1;
			for (const auto& [channels, pixel] : pixels)
			{
				if (!used && (channels.first == channel || channels.second == channel))
					status = std::max(status, pixel.status == Estimated ? 3 : 2);
			}
			reference = (channel == fit.reference);
			stripTree.Fill();
			if (used)
				table << channel << '\t' << offset << '\t' << gain << '\n';
		}
	}
	if (fit.strips) table.close();
	if (fit.strips && !table)
		throw std::runtime_error("写入系数表失败 " + datPath.string());

	xy.Write();
	xyh.Write();
	pixelTree.Write();
	stripTree.Write();
	TParameter<int>("reference_x", fit.reference).Write();
	TParameter<int>("constraint_ndf", fit.strips ? 2 * fit.pixels - 2 * (fit.strips - 1) : 0).Write();
	TParameter<double>("robust_objective", fit.objective).Write();
	TParameter<Long64_t>("entries_scanned", entries).Write();
	TObjString(fitError.c_str()).Write("fit_error");
	TObjString("0=too_few_events;1=narrow_energy_range;2=fit_failed;3=poor_fit;4=estimated").Write("pixel_status_codes");
	TObjString("0=normalized;1=no_pairs;2=no_valid_pixel;3=outside_reference_component_or_global_failure").Write("strip_status_codes");
	std::ostringstream config;
	config << "setup=" << NORM_SETUP_PATH << "\n"
	       << "amplitude_min=" << NormConfig::MIN_X << ',' << NormConfig::MIN_Y << ',' << NormConfig::MIN_YH
	       << "\namplitude_max=" << NormConfig::MAX_AMPLITUDE << "\ncoin_ns=" << COIN_SI_WINDOW
	       << "\nreject_pileup=" << NormConfig::REJECT_PILEUP << "\nmin_pixel_events=" << NormConfig::MIN_PIXEL_EVENTS
	       << "\nmin_x_span=" << NormConfig::MIN_X_SPAN << "\nmax_sample=" << NormConfig::MAX_PIXEL_SAMPLE
	       << "\nseed=" << NormConfig::SEED << "\nseed_trials=" << NormConfig::SEED_TRIALS
	       << "\nnoise_floor=" << NormConfig::NOISE_FLOOR
	       << "\nmin_inlier_fraction=" << NormConfig::MIN_INLIER_FRACTION
	       << "\nmax_relative_slope_error=" << NormConfig::MAX_REL_SLOPE_ERROR
	       << "\nrequested_reference_x=" << NormConfig::REFERENCE_X << "\nsigma_back_over_x="
	       << NormConfig::SIGMA_Y_OVER_X << ',' << NormConfig::SIGMA_YH_OVER_X << "\nmax_entries=" << maxEntries << "\n";
	for (const auto& input : inputs) config << "input=" << fs::absolute(input).string() << '\n';
	TObjString(config.str().c_str()).Write("configuration");
	output.Flush();
	if (output.TestBit(TFile::kWriteError))
		throw std::runtime_error("写入诊断 ROOT 文件失败 " + rootPath.string());
}

void PublishOutput(const std::vector<fs::path>& inputs, const Pixels& pixels, const StripFit& fit,
                   const fs::path& outputDir, const std::string& basename, Long64_t maxEntries, Long64_t entries, const std::string& fitError)
{
	fs::create_directories(outputDir);
	const fs::path rootPath = outputDir / (basename + ".root");
	const fs::path datPath = outputDir / (basename + ".dat");
	if (fs::exists(rootPath) || fs::exists(datPath))
		throw std::runtime_error("输出文件已存在 " + basename);

	std::string pattern = (outputDir / ".normalize-XXXXXX").string();
	std::vector<char> tempName(pattern.begin(), pattern.end());
	tempName.push_back(0);
	if (!mkdtemp(tempName.data()))
		throw std::runtime_error("无法创建临时输出目录");
	const fs::path temp = tempName.data();
	bool rootPublished = false;
	try
	{
		WriteOutput(inputs, pixels, fit, temp / "diagnostics.root", temp / "parameters.dat", maxEntries, entries, fitError);
		fs::create_hard_link(temp / "diagnostics.root", rootPath);
		rootPublished = true;
		if (fit.strips) fs::create_hard_link(temp / "parameters.dat", datPath);
	}
	catch (...)
	{
		std::error_code ignored;
		if (rootPublished)
			fs::remove(rootPath, ignored);
		fs::remove_all(temp, ignored);
		throw;
	}
	fs::remove_all(temp);
	std::cout << "诊断 " << fs::absolute(rootPath) << std::endl;
	if (fit.strips) std::cout << "系数 " << fs::absolute(datPath) << std::endl;
}

int ParseRun(const std::string& text)
{
	if (text.empty() || text.size() > 5 || text.find_first_not_of("0123456789") != std::string::npos)
		throw std::runtime_error("run 号应为 0 到 99999 的整数");
	return std::stoi(text);
}
} // namespace

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		std::cout << "用法 ./Normalize_DSSD start stop [--input-dir DIR] [--prefix ID]\n"
		             "       [--output-dir DIR] [--max-entries N]\n";
		return argc == 2 && std::string(argv[1]) == "--help" ? 0 : 1;
	}
	try
	{
		const int first = ParseRun(argv[1]), last = ParseRun(argv[2]);
		if (last < first) throw std::runtime_error("终止 run 小于起始 run");
		fs::path inputDir = OUTPUT_FILE_PATH, outputDir = ".";
		std::string prefix = FILENAME_PREFIX;
		Long64_t maxEntries = 0;
		for (int i = 3; i < argc; i += 2)
		{
			if (i + 1 >= argc) throw std::runtime_error("选项缺少参数");
			const std::string option = argv[i];
			if (option == "--input-dir") inputDir = argv[i + 1];
			else if (option == "--output-dir") outputDir = argv[i + 1];
			else if (option == "--prefix") prefix = argv[i + 1];
			else if (option == "--max-entries")
			{
				const std::string value = argv[i + 1];
				if (value.empty() || value.find_first_not_of("0123456789") != std::string::npos)
					throw std::runtime_error("max-entries 应为非负整数");
				maxEntries = std::stoll(value);
			}
			else throw std::runtime_error("未知选项 " + option);
		}
		if (prefix.find('/') != std::string::npos || prefix.find(static_cast<char>(92)) != std::string::npos)
			throw std::runtime_error("prefix 应为文件名前缀");
		if (!(NormConfig::SIGMA_Y_OVER_X > 0 && NormConfig::SIGMA_YH_OVER_X > 0) ||
		    (NormConfig::MAX_PIXEL_SAMPLE && NormConfig::MAX_PIXEL_SAMPLE < NormConfig::MIN_PIXEL_EVENTS))
			throw std::runtime_error("噪声比或抽样上限非法");
		std::vector<fs::path> inputs;
		for (int run = first; run <= last; ++run)
		{
			const fs::path input = inputDir / Form("%s%05d_map.root", prefix.c_str(), run);
			if (fs::exists(input)) inputs.push_back(input);
			else std::cerr << "缺失 run " << input << '\n';
		}
		if (inputs.empty()) throw std::runtime_error("找不到输入 ROOT 文件");
		Pixels pixels;
		std::mt19937_64 random(NormConfig::SEED);
		ULong64_t pairs = 0;
		const auto entries = ReadPairs(inputs, maxEntries, [&](const Hit& x, const Hit& back)
		{
			auto& pixel = pixels[{x.channel, back.channel}];
			++pixel.events; ++pairs;
			const PairSample sample{x.amplitude, back.amplitude};
			if (!NormConfig::MAX_PIXEL_SAMPLE || pixel.sample.size() < NormConfig::MAX_PIXEL_SAMPLE)
				pixel.sample.push_back(sample);
			else
			{
				std::uniform_int_distribution<ULong64_t> choose(0, pixel.events - 1);
				const auto i = choose(random);
				if (i < pixel.sample.size()) pixel.sample[i] = sample;
			}
		});
		std::cout << "扫描 " << entries << "，有效配对 " << pairs << "，像素 " << pixels.size() << std::endl;
		EstimatePixels(pixels);
		StripFit fit;
		std::string fitError;
		try { fit = FitStrips(pixels); }
		catch (const std::exception& error) { fitError = error.what(); }
		const std::string basename = Form("DSSD_Normalize_Robust_%s_%05d_%05d", prefix.c_str(), first, last);
		PublishOutput(inputs, pixels, fit, outputDir, basename, maxEntries, entries, fitError);
		if (!fitError.empty())
		{
			std::cerr << fitError << '\n';
			return 3;
		}
		const int total = (kX.last - kX.first + 1) + (kY.last - kY.first + 1) + (kYH.last - kYH.first + 1);
		std::cout << "参考 X 全局道 " << fit.reference << "，拟合像素 " << fit.pixels
		          << "，已归一化道 " << fit.strips << '/' << total << '\n';
		if (fit.strips < total) std::cout << "其余通道见 strip_fit 的 used=false；系数表只含有效道。\n";
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}
