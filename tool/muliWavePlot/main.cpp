#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include "CAEN_DGTZ.h"
#include "Data_Map.h"

namespace
{

constexpr double ENERGY_THRESHOLD = 1000.0;
constexpr std::size_t DETECTOR_COUNT = 4;

struct WaveEntry
{
	int runnum;
	Long64_t map_entry;
};

struct DetectorSpec
{
	const char *panel_name;
	const char *wave_name;
	int channel_offset;
};

struct DecodedWave
{
	std::string name;
	int idx;
	int channel;
	double energy;
	std::vector<double> samples;
};

struct PanelResult
{
	DetectorSpec spec;
	int multiplicity;
	int energy_candidate_count;
	std::vector<DecodedWave> waves;
	std::vector<std::string> errors;
	std::string unavailable_reason;
};

struct WaveCandidate
{
	std::size_t panel_index;
	Particle_Data *data;
	int idx;
};

const std::array<DetectorSpec, DETECTOR_COUNT> DETECTOR_SPECS = {{
	{"XWave", "XWave", DSSD_X_CH_LOW},
	{"YWave", "YWave", DSSD_Y_CH_LOW},
	{"YHWave", "YHWave", DSSD_YH_CH_LOW},
	{"SSDWave", "SSDWave", SSD_CH_LOW}
}};

const std::array<int, 10> WAVE_COLORS = {{
	kBlue + 1,
	kRed + 1,
	kGreen + 2,
	kMagenta + 1,
	kOrange + 7,
	kCyan + 2,
	kViolet + 1,
	kAzure + 7,
	kPink + 7,
	kBlack
}};

std::string Trim(const std::string &input)
{
	const std::string whitespace = " \t\r\n";
	const std::size_t first = input.find_first_not_of(whitespace);
	if(first == std::string::npos)
		return "";
	const std::size_t last = input.find_last_not_of(whitespace);
	return input.substr(first, last - first + 1);
}

bool ParseSignedInteger(const std::string &token, long long &value)
{
	errno = 0;
	char *end = nullptr;
	const long long parsed = std::strtoll(token.c_str(), &end, 10);
	if(errno == ERANGE || end == token.c_str() || *end != '\0')
		return false;
	value = parsed;
	return true;
}

bool ReadWaveEntries(const std::string &path, std::vector<WaveEntry> &entries, bool &had_error)
{
	std::ifstream input(path.c_str());
	if(!input.is_open())
	{
		std::cerr << "Cannot open wave-entry file: " << path << std::endl;
		return false;
	}

	std::string line;
	std::size_t line_number = 0;
	while(std::getline(input, line))
	{
		++line_number;
		line = Trim(line);
		if(line.empty() || line[0] == '#')
			continue;

		std::istringstream stream(line);
		std::string run_token;
		std::string entry_token;
		stream >> run_token >> entry_token;

		if(run_token == "runnum" && entry_token == "map_entry")
			continue;

		if(run_token.empty() || entry_token.empty())
		{
			std::cerr << path << ':' << line_number
				<< ": expected two columns: runnum map_entry" << std::endl;
			had_error = true;
			continue;
		}

		std::string extra_token;
		if(stream >> extra_token)
		{
			if(extra_token.empty() || extra_token[0] != '#')
			{
				std::cerr << path << ':' << line_number
					<< ": unexpected data after map_entry" << std::endl;
				had_error = true;
				continue;
			}
		}

		long long run_value = 0;
		long long entry_value = 0;
		if(!ParseSignedInteger(run_token, run_value)
			|| !ParseSignedInteger(entry_token, entry_value)
			|| run_value < 0
			|| run_value > INT_MAX
			|| entry_value < 0)
		{
			std::cerr << path << ':' << line_number
				<< ": invalid non-negative runnum or map_entry" << std::endl;
			had_error = true;
			continue;
		}

		WaveEntry entry;
		entry.runnum = static_cast<int>(run_value);
		entry.map_entry = static_cast<Long64_t>(entry_value);
		entries.push_back(entry);
	}

	if(entries.empty())
	{
		std::cerr << "No valid runnum/map_entry records found in " << path << std::endl;
		return false;
	}
	return true;
}

std::string BuildRunPath(const char *directory, int runnum, const char *suffix)
{
	std::ostringstream path;
	path << directory;
	const std::string directory_string = directory;
	if(directory_string.empty() || directory_string[directory_string.size() - 1] != '/')
		path << '/';
	path << FILENAME_PREFIX << std::setw(5) << std::setfill('0') << runnum << suffix;
	return path.str();
}

void SetAllPanelsUnavailable(std::array<PanelResult, DETECTOR_COUNT> &panels,
	const std::string &reason)
{
	for(std::size_t i = 0; i < panels.size(); ++i)
		panels[i].unavailable_reason = reason;
}

Particle_Data *GetDetectorData(Data_Map &map_reader, std::size_t panel_index)
{
	switch(panel_index)
	{
		case 0:
			return &map_reader.Data_DSSDX;
		case 1:
			return &map_reader.Data_DSSDY;
		case 2:
			return &map_reader.Data_DSSDYH;
		case 3:
			return &map_reader.Data_SSD;
		default:
			return nullptr;
	}
}

// The GUI decoder performs this same lookup, but its trace-only loop does not
// advance ch_idx. Keep the shared low-level decoder and perform the aggregate
// navigation here so DGTZ_ChAggregateNo values greater than zero work too.
int DecodeTrace(CAEN_DGTZ &decoder,
	uint64_t entry,
	int dgtz_type,
	int board_aggregate_number,
	int current_channel,
	int channel_aggregate_number,
	DPP_PHA_Event_t &pha_data,
	uint16_t *trace)
{
	if(dgtz_type < 0 || dgtz_type >= 2 || current_channel < 0
		|| current_channel >= MAXNCH_V1724 || board_aggregate_number < 0
		|| channel_aggregate_number < 0 || trace == nullptr)
		return -1;

	decoder.SetGetWave(true);
	if(decoder.GetEntry(static_cast<Long64_t>(entry)) <= 0)
		return -1;

	if(decoder.DG_data[dgtz_type].nDG == 0)
		return -1;

	uint32_t bank_size = static_cast<uint32_t>(decoder.DG_data[dgtz_type].nDG);
	if(bank_size < 10)
		return -1;

	std::memset(&pha_data, 0, sizeof(DPP_PHA_Event_t));
	uint32_t *board = &(decoder.DG_data[dgtz_type].DG[2]);
	bank_size -= 2;

	uint32_t board_aggregate_count = 0;
	if(decoder.GetBoardAggregateNum(board, bank_size, &board_aggregate_count) != 0)
		return -1;

	for(uint32_t board_index = 0; board_index < board_aggregate_count; ++board_index)
	{
		Board_Aggregate_Info_t board_info;
		if(decoder.GetBoardAggregateInfo(board, &board_info) != 0
			|| board_info.BoardAggregateSize == 0)
			return -1;

		if(static_cast<int>(board_index) == board_aggregate_number)
		{
			uint32_t *channel = board + 4;
			for(int channel_index = 0; channel_index < MAXNCH_V1724; ++channel_index)
			{
				if(!((board_info.ChannelMask >> channel_index) & 0x1))
					continue;

				Channel_Aggregate_Info_t channel_info;
				if(decoder.GetChannelAggregateInfo(channel, dgtz_type, &channel_info) != 0
					|| channel_info.ChannelAggregateSize == 0)
					return -1;

				if(channel_index == current_channel)
				{
					uint32_t *event_data = channel + (channel_info.FI ? 2 : 1);
					const uint32_t data_size = (channel_info.ET ? 1U : 0U)
						+ (channel_info.ES ? channel_info.SampleNumber / 2U : 0U)
						+ (channel_info.EE ? 1U : 0U)
						+ (channel_info.E2 == 1 ? 1U : 0U);
					uint32_t data_position = channel_info.FI ? 2U : 1U;
					int aggregate_index = 0;

					if(data_size == 0)
						return -1;

					while(data_position < channel_info.ChannelAggregateSize)
					{
						if(aggregate_index == channel_aggregate_number)
						{
							int decode_result = -1;
							if(dgtz_type == 0)
								decode_result = decoder.GetDPPPHAEvent_V1724(
									event_data, data_size, &pha_data, channel_info);
							else if(dgtz_type == 1)
								decode_result = decoder.GetDPPPHAEvent_V1730(
									event_data, data_size, &pha_data, channel_info);

							if(decode_result != 0)
								return -1;
							pha_data.SampleNo = static_cast<uint16_t>(
								std::min<uint32_t>(channel_info.SampleNumber, NSAMPLE_MAX));
							decoder.Get_Trace(trace);
							return 0;
						}

						event_data += data_size;
						data_position += data_size;
						++aggregate_index;
					}
				}

				channel += channel_info.ChannelAggregateSize;
			}
		}

		board += board_info.BoardAggregateSize;
	}

	return -1;
}

void DrawMessage(const PanelResult &panel, const std::vector<std::string> &messages)
{
	TLatex text;
	text.SetNDC(kTRUE);
	text.SetTextAlign(22);
	text.SetTextFont(42);
	text.SetTextSize(0.055);
	text.DrawLatex(0.5, 0.82, panel.spec.panel_name);

	text.SetTextSize(0.04);
	double y = 0.58;
	for(std::size_t i = 0; i < messages.size() && i < 6; ++i)
	{
		text.DrawLatex(0.5, y, messages[i].c_str());
		y -= 0.08;
	}
}

void DrawPanel(const PanelResult &panel)
{
	gPad->SetGrid();

	if(!panel.unavailable_reason.empty())
	{
		DrawMessage(panel, {panel.unavailable_reason});
		return;
	}
	if(panel.multiplicity == 0)
	{
		DrawMessage(panel, {"No data (mul = 0)"});
		return;
	}
	if(panel.energy_candidate_count == 0)
	{
		std::ostringstream message;
		message << "No waveform with E > " << ENERGY_THRESHOLD;
		DrawMessage(panel, {message.str()});
		return;
	}
	if(panel.waves.empty())
	{
		std::vector<std::string> messages;
		messages.push_back("No waveform could be decoded");
		for(std::size_t i = 0; i < panel.errors.size() && i < 4; ++i)
			messages.push_back(panel.errors[i]);
		DrawMessage(panel, messages);
		return;
	}

	double minimum = std::numeric_limits<double>::max();
	double maximum = std::numeric_limits<double>::lowest();
	std::size_t maximum_samples = 0;
	for(std::size_t wave_index = 0; wave_index < panel.waves.size(); ++wave_index)
	{
		const DecodedWave &wave = panel.waves[wave_index];
		maximum_samples = std::max(maximum_samples, wave.samples.size());
		for(std::size_t sample = 0; sample < wave.samples.size(); ++sample)
		{
			minimum = std::min(minimum, wave.samples[sample]);
			maximum = std::max(maximum, wave.samples[sample]);
		}
	}

	if(!std::isfinite(minimum) || !std::isfinite(maximum))
	{
		DrawMessage(panel, {"Decoded waveform contains no samples"});
		return;
	}

	const double span = std::max(1.0, maximum - minimum);
	const double margin = span * 0.08;
	const double x_max = maximum_samples > 1 ? static_cast<double>(maximum_samples - 1) : 1.0;
	const TString frame_title = TString::Format("%s;Sample;ADC", panel.spec.panel_name);
	gPad->DrawFrame(0.0, minimum - margin, x_max, maximum + margin, frame_title.Data());

	const int column_count = std::min<int>(3,
		std::max<int>(1, static_cast<int>((panel.waves.size() + 7) / 8)));
	const int row_count = static_cast<int>((panel.waves.size() + column_count - 1) / column_count);
	const double legend_y_low = std::max(0.12, 0.90 - row_count * 0.045);
	TLegend legend(0.48, legend_y_low, 0.89, 0.90);
	legend.SetNColumns(column_count);
	legend.SetBorderSize(0);
	legend.SetFillStyle(0);
	legend.SetTextSize(0.027);

	for(std::size_t wave_index = 0; wave_index < panel.waves.size(); ++wave_index)
	{
		const DecodedWave &wave = panel.waves[wave_index];
		TGraph graph(static_cast<int>(wave.samples.size()));
		graph.SetName(wave.name.c_str());
		graph.SetTitle(wave.name.c_str());
		graph.SetLineColor(WAVE_COLORS[wave_index % WAVE_COLORS.size()]);
		graph.SetLineWidth(2);
		for(std::size_t sample = 0; sample < wave.samples.size(); ++sample)
			graph.SetPoint(static_cast<int>(sample), static_cast<double>(sample), wave.samples[sample]);

		TGraph *drawn_graph = static_cast<TGraph *>(graph.DrawClone("L SAME"));
		drawn_graph->SetName(wave.name.c_str());
		const TString legend_text = TString::Format("%s  Ch=%d  E=%.1f",
			wave.name.c_str(), wave.channel, wave.energy);
		legend.AddEntry(drawn_graph, legend_text.Data(), "l");
	}
	legend.DrawClone();

	if(!panel.errors.empty())
	{
		TLatex warning;
		warning.SetNDC(kTRUE);
		warning.SetTextColor(kRed + 1);
		warning.SetTextSize(0.032);
		const TString message = TString::Format("%zu waveform(s) failed to decode",
			panel.errors.size());
		warning.DrawLatex(0.12, 0.12, message.Data());
	}
}

void WriteCanvas(TFile &output,
	std::size_t chain_number,
	const WaveEntry &entry,
	const std::array<PanelResult, DETECTOR_COUNT> &panels)
{
	output.cd();
	const TString canvas_name = TString::Format("ChainNum_%zu", chain_number);
	const TString canvas_title = TString::Format(
		"%s: run=%d, map_entry=%lld",
		canvas_name.Data(), entry.runnum, static_cast<long long>(entry.map_entry));
	TCanvas canvas(canvas_name.Data(), canvas_title.Data(), 1400, 1000);
	canvas.Divide(2, 2);

	for(std::size_t panel_index = 0; panel_index < panels.size(); ++panel_index)
	{
		canvas.cd(static_cast<int>(panel_index + 1));
		DrawPanel(panels[panel_index]);
	}

	canvas.Modified();
	canvas.Update();
	canvas.Write(canvas_name.Data(), TObject::kOverwrite);
}

bool ProcessWaveEntry(const WaveEntry &entry,
	std::size_t chain_number,
	Data_Map &map_reader,
	CAEN_DGTZ &decoder,
	TFile &output)
{
	std::array<PanelResult, DETECTOR_COUNT> panels;
	for(std::size_t i = 0; i < panels.size(); ++i)
	{
		panels[i].spec = DETECTOR_SPECS[i];
		panels[i].multiplicity = -1;
		panels[i].energy_candidate_count = 0;
	}

	bool success = true;
	const std::string map_path = BuildRunPath(OUTPUT_FILE_PATH, entry.runnum, "_map.root");
	std::unique_ptr<TFile> map_file(TFile::Open(map_path.c_str(), "READ"));
	if(!map_file || map_file->IsZombie())
	{
		std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
			<< ": cannot open map file " << map_path << std::endl;
		SetAllPanelsUnavailable(panels, "Map file unavailable");
		WriteCanvas(output, chain_number, entry, panels);
		return false;
	}

	TTree *map_tree = nullptr;
	map_file->GetObject("tr_map", map_tree);
	if(map_tree == nullptr)
	{
		std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
			<< ": tree tr_map is missing in " << map_path << std::endl;
		SetAllPanelsUnavailable(panels, "Tree tr_map unavailable");
		WriteCanvas(output, chain_number, entry, panels);
		return false;
	}

	const Long64_t map_entry_count = map_tree->GetEntries();
	if(entry.map_entry < 0 || entry.map_entry >= map_entry_count)
	{
		std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
			<< ": map entry is outside [0, " << map_entry_count << ')' << std::endl;
		SetAllPanelsUnavailable(panels, "Map entry out of range");
		WriteCanvas(output, chain_number, entry, panels);
		return false;
	}

	if(map_reader.Init(map_tree) < 0 || map_reader.GetEntry(entry.map_entry) <= 0)
	{
		std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
			<< ": failed to read map entry" << std::endl;
		SetAllPanelsUnavailable(panels, "Map entry could not be read");
		WriteCanvas(output, chain_number, entry, panels);
		return false;
	}

	std::vector<WaveCandidate> candidates;
	for(std::size_t panel_index = 0; panel_index < panels.size(); ++panel_index)
	{
		Particle_Data *data = GetDetectorData(map_reader, panel_index);
		PanelResult &panel = panels[panel_index];
		panel.multiplicity = data == nullptr ? 0 : static_cast<int>(data->mul);
		if(data == nullptr)
			continue;

		int multiplicity = panel.multiplicity;
		if(multiplicity > MAX_SIGNAL_EVENT)
		{
			std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
				<< " detector=" << panel.spec.panel_name << ": mul=" << multiplicity
				<< " exceeds MAX_SIGNAL_EVENT=" << MAX_SIGNAL_EVENT << std::endl;
			multiplicity = MAX_SIGNAL_EVENT;
			success = false;
		}

		for(int idx = 0; idx < multiplicity; ++idx)
		{
			if(data->E[idx] <= ENERGY_THRESHOLD)
				continue;
			++panel.energy_candidate_count;
			candidates.push_back({panel_index, data, idx});
		}
	}

	if(!candidates.empty())
	{
		const std::string raw_path = BuildRunPath(INPUT_FILE_PATH, entry.runnum, ".root");
		std::unique_ptr<TFile> raw_file(TFile::Open(raw_path.c_str(), "READ"));
		TTree *trigger_tree = nullptr;
		if(raw_file && !raw_file->IsZombie())
			raw_file->GetObject("Trigger", trigger_tree);

		if(!raw_file || raw_file->IsZombie() || trigger_tree == nullptr)
		{
			std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
				<< ": cannot open Trigger tree in " << raw_path << std::endl;
			for(std::size_t candidate_index = 0; candidate_index < candidates.size(); ++candidate_index)
			{
				const WaveCandidate &candidate = candidates[candidate_index];
				PanelResult &panel = panels[candidate.panel_index];
				const std::string wave_name = std::string(panel.spec.wave_name)
					+ "_" + std::to_string(candidate.idx);
				panel.errors.push_back(wave_name + ": raw file unavailable");
				std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
					<< " detector=" << panel.spec.panel_name << " idx=" << candidate.idx
					<< ": raw file or Trigger tree unavailable" << std::endl;
			}
			success = false;
		}
		else
		{
			decoder.Init(trigger_tree);
			for(std::size_t candidate_index = 0; candidate_index < candidates.size(); ++candidate_index)
			{
				const WaveCandidate &candidate = candidates[candidate_index];
				PanelResult &panel = panels[candidate.panel_index];
				Particle_Data &data = *candidate.data;
				const int idx = candidate.idx;
				const std::string wave_name = std::string(panel.spec.wave_name)
					+ "_" + std::to_string(idx);
				const int channel = static_cast<int>(data.Ch[idx]) + panel.spec.channel_offset;
				const int dgtz_type = static_cast<int>(data.DGTZ_Type[idx]);
				const int decoder_channel = channel % MAXNCH_V1724;

				std::array<uint16_t, NSAMPLE_MAX> trace = {{0}};
				DPP_PHA_Event_t pha_data = {};
				const int decode_result = DecodeTrace(
					decoder,
					data.DGTZ_EntryNo[idx],
					dgtz_type,
					static_cast<int>(data.DGTZ_BdAggregateNo[idx]),
					decoder_channel,
					static_cast<int>(data.DGTZ_ChAggregateNo[idx]),
					pha_data,
					trace.data());

				const std::size_t sample_count = std::min<std::size_t>(pha_data.SampleNo, NSAMPLE_MAX);
				if(decode_result != 0 || sample_count == 0)
				{
					panel.errors.push_back(wave_name + ": decode failed");
					std::cerr << "runnum=" << entry.runnum << " map_entry=" << entry.map_entry
						<< " detector=" << panel.spec.panel_name << " idx=" << idx
						<< ": waveform decode failed" << std::endl;
					success = false;
					continue;
				}

				DecodedWave wave;
				wave.name = wave_name;
				wave.idx = idx;
				wave.channel = channel;
				wave.energy = data.E[idx];
				wave.samples.resize(sample_count);
				for(std::size_t sample = 0; sample < sample_count; ++sample)
					wave.samples[sample] = static_cast<double>(trace[sample]);
				panel.waves.push_back(std::move(wave));
			}
		}
	}

	WriteCanvas(output, chain_number, entry, panels);
	return success;
}

} // namespace

int main(int argc, char *argv[])
{
	if(argc > 3)
	{
		std::cerr << "Usage: " << argv[0] << " [waveEntry.dat] [wave.root]" << std::endl;
		return 2;
	}

	const std::string input_path = argc >= 2 ? argv[1] : "waveEntry.dat";
	const std::string output_path = argc >= 3 ? argv[2] : "wave.root";

	bool had_error = false;
	std::vector<WaveEntry> entries;
	if(!ReadWaveEntries(input_path, entries, had_error))
		return 1;

	gROOT->SetBatch(kTRUE);
	TH1::AddDirectory(kFALSE);
	std::unique_ptr<TFile> output(TFile::Open(output_path.c_str(), "RECREATE"));
	if(!output || output->IsZombie())
	{
		std::cerr << "Cannot create output ROOT file: " << output_path << std::endl;
		return 1;
	}

	std::unique_ptr<Data_Map> map_reader(new Data_Map());
	std::unique_ptr<CAEN_DGTZ> decoder(new CAEN_DGTZ());
	for(std::size_t entry_index = 0; entry_index < entries.size(); ++entry_index)
	{
		std::cout << "Processing ChainNum_" << entry_index + 1
			<< " (run=" << entries[entry_index].runnum
			<< ", map_entry=" << entries[entry_index].map_entry << ")" << std::endl;
		if(!ProcessWaveEntry(entries[entry_index], entry_index + 1,
			*map_reader, *decoder, *output))
			had_error = true;
	}

	output->Write();
	output->Close();
	std::cout << "Wrote " << entries.size() << " canvases to " << output_path << std::endl;
	return had_error ? 1 : 0;
}
