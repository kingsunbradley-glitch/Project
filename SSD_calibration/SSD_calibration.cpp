// 编译命令:
// g++ SSD_calibration.cpp `root-config --cflags --libs` -lSpectrum -lRDataFrame -lTreePlayer -o SSD_calibration
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>
#include <memory> 
#include <stdexcept> 
#include <numeric> 
#include <thread> // 引入 thread 头文件

// RDataFrame 相关的头文件
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDFHelpers.hxx" // For AddProgressBar, RunGraphs

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TCut.h"
#include "TGraph.h"
#include "TF1.h"
#include "TAxis.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TStyle.h"

// --- 包含项目配置文件 ---
#include "config.h"

/* 函数：验证 config.h 中的参数是否合理 (保持不变)
void validate_config() {
    if (peak_search_sigma <= 0) {
        throw std::runtime_error("Error: peak_search_sigma must be positive.");
    }
    if (peak_search_threshold <= 0 || peak_search_threshold >= 1) {
        throw std::runtime_error("Error: peak_search_threshold must be between 0 and 1.");
    }
    if (MIN_HIST_ENTRIES <= 0) {
        throw std::runtime_error("Error: MIN_HIST_ENTRIES must be positive.");
    }
}
*/

ErrorCode calibrate_ssd_professional() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gSystem->mkdir(output_plot_dir, kTRUE);

    // --- 步骤 0: 启用多核并行 ---
    // 修正错误 1: 避免 GetImplicitMTRuntime() 的兼容性问题
    ROOT::EnableImplicitMT(num_threads);
    // 使用核数目
    std::cout << "\n--> Starting RDataFrame-based professional calibration process (using up to " 
              << num_threads << " threads)..." << std::endl;

    // --- 1. RDataFrame 准备 ---
    ROOT::RDataFrame df(tree_name, input_root_file);
    auto n_entries = df.Count();
    std::cout << "    Total entries in TTree: " << *n_entries << std::endl;
    if (*n_entries == 0) {
        std::cerr << "Error: No entries found in TTree." << std::endl;
        return FILE_NOT_FOUND;
    }

    // 定义一个 Lambda 函数，用于计算 Energy_Sum
    auto df_base = df.Define("Energy_Sum", [](const ROOT::RVec<double>& dssd_e, const ROOT::RVec<double>& ssd_e, UInt_t chain_mul) {
        ROOT::RVec<double> sum_e;
        for (UInt_t i = 0; i < chain_mul; ++i) {
             if (dssd_e[i] > 0 && ssd_e[i] > 0) {
                 sum_e.push_back(dssd_e[i] + ssd_e[i]);
             } else {
                 sum_e.push_back(-999.0); // 标记无效 hits
             }
        }
        return sum_e;
    }, {"DSSD_E", "SSD_E", "chain_mul"});

    // --- 2. 并行绘制 48 个能谱 (TH1D) ---
    std::vector<ROOT::RDF::RResultPtr<TH1D>> hists;
    for (int pos = 0; pos < NUM_SSD_POS; ++pos) {
        
        // --- 修正错误 3 和 4 的逻辑 ---
        // 步骤 1: 根据 SSD_Pos 筛选出对应的 Energy_Sum 和 SSD_Pos 数组元素
        auto df_pos_data = df_base.Define(TString::Format("Pos%d_Energy", pos), [pos](const ROOT::RVec<double>& energy_sum, const ROOT::RVec<double>& ssd_pos_vec, UInt_t chain_mul) {
            ROOT::RVec<double> pos_energy;
            for (UInt_t i = 0; i < chain_mul; ++i) {
                // 确保索引有效，且 Pos 匹配
                if (i < ssd_pos_vec.size() && TMath::Nint(ssd_pos_vec[i]) == pos && energy_sum[i] > 0) {
                    pos_energy.push_back(energy_sum[i]);
                }
            }
            return pos_energy;
        }, {"Energy_Sum", "SSD_Pos", "chain_mul"});

        // 步骤 2: 应用全局时间削减 (修正错误 2: global_event_cut 检查)
        // 注意：这里需要确保 global_event_cut 已经被正确包含在 config.h 并定义
       ROOT::RDF::RNode df_to_plot(df_pos_data);
        if (!std::string(global_event_cut).empty()) {
        df_to_plot = df_pos_data.Filter(global_event_cut, TString::Format("Global Cut for Pos %d", pos));
        }



        // 步骤 3: 仅对包含至少一个目标 Pos hit 的事件进行绘图
        hists.push_back(
    // 调用 .Data() 将 TString 转换为 const char*，让 RDataFrame 识别为表达式
    df_to_plot.Filter(TString::Format("Pos%d_Energy.size() > 0", pos).Data()) 
    .Histo1D(
        {TString::Format("h_sum_E_pos%d", pos),
         TString::Format("Energy Spectrum for Pos %d;Energy Sum (channels);Counts", pos),
         700, 3000, 10000},
        TString::Format("Pos%d_Energy", pos)
    )
);
    }
    
    // --- 3. 运行 RDataFrame (多核) ---
    std::cout << "\n    Executing RDataFrame tasks (plotting " << NUM_SSD_POS << " histograms)..." << std::endl;

    // 修正错误 5: 将 RResultPtr 转换为 RResultHandle 的 vector
    std::vector<ROOT::RDF::RResultHandle> tasks;
    tasks.reserve(hists.size());
    for (auto& h : hists) {
        tasks.push_back(h);
    }
    
    // 添加进度条并触发事件循环
    ROOT::RDF::Experimental::AddProgressBar(df); 
    ROOT::RDF::RunGraphs(tasks);
    std::cout << "    RDataFrame execution finished." << std::endl;


    // --- 4. 单核寻峰和拟合 ---
    std::map<int, double> k_params, b_params;
    std::map<int, int> confidence;
    auto canvas = std::make_unique<TCanvas>("canvas", "Calibration Canvas", 1200, 900);
    auto s = std::make_unique<TSpectrum>();

    std::cout << "\n--> Starting single-core fitting process (TSpectrum + TGraph)..." << std::endl;

    for (int pos = 0; pos < NUM_SSD_POS; ++pos) {
        TH1D* h_sum_E = hists[pos].GetPtr(); // 从 RResultPtr 获取指针
        
        std::cout << "    Processing Pos " << pos << " (Entries: " << h_sum_E->GetEntries() << ")..." << std::endl;

        if (h_sum_E->GetEntries() < MIN_HIST_ENTRIES) {
            std::cerr << "    Warning: Insufficient data for Pos " << pos << ". Skipping." << std::endl;
            k_params[pos] = 1.0; b_params[pos] = 0.0; confidence[pos] = -1; 
            continue;
        }

        // --- 步骤 4.1: 背景扣除 ---
        std::unique_ptr<TH1D> h_sub; 
        if (enable_background_subtraction) {
            canvas->Clear();
            TH1 *h_bg_raw = s->Background(h_sum_E, background_iterations, "same");
            auto h_bg = std::unique_ptr<TH1>(h_bg_raw); 
            
            h_sum_E->Draw();
            h_bg->SetLineColor(kRed);  
            h_bg->Draw("same"); 
            auto leg_bg = std::make_unique<TLegend>(0.7, 0.75, 0.9, 0.9);
            leg_bg->AddEntry(h_sum_E, "Original Spectrum", "l");
            leg_bg->AddEntry(h_bg.get(), "Estimated Background", "l");
            leg_bg->Draw();
            canvas->SaveAs(TString::Format("%s/pos_%d_1_background.png", output_plot_dir, pos));
            
            h_sub.reset((TH1D*)h_sum_E->Clone(TString::Format("h_sub_pos%d", pos)));
            h_sub->Add(h_bg.get(), -1);
        } else {
            h_sub.reset((TH1D*)h_sum_E->Clone(TString::Format("h_sub_pos%d", pos)));
        }


        // --- 步骤 4.2: 自动寻峰 ---
        Int_t n_found = s->Search(h_sub.get(), peak_search_sigma, "", peak_search_threshold);
        Double_t *found_peaks_x = s->GetPositionX();
        
        std::vector<double> cal_points, ref_points;
        std::vector<std::unique_ptr<TMarker>> used_markers, other_markers;
        std::vector<bool> is_peak_used(n_found, false);

        for(const auto& p_interest : peaks_of_interest) {
            double best_peak_x = -1;
            double max_peak_y = 0;
            int best_peak_idx = -1;
            for (int i = 0; i < n_found; ++i) {
                if (found_peaks_x[i] > p_interest.win_min && found_peaks_x[i] < p_interest.win_max) {
                    double current_peak_y = h_sub->GetBinContent(h_sub->FindBin(found_peaks_x[i]));
                    if (current_peak_y > max_peak_y) {
                        max_peak_y = current_peak_y;
                        best_peak_x = found_peaks_x[i];
                        best_peak_idx = i;
                    }
                }
            }
            if (best_peak_x > 0) {
                cal_points.push_back(best_peak_x);
                ref_points.push_back(p_interest.ref_E);
                if (best_peak_idx != -1) is_peak_used[best_peak_idx] = true;
            }
        }
        
        canvas->Clear();
        h_sub->SetTitle(TString::Format("Peak Search for Pos %d (Found %d peaks)", pos, n_found));
        h_sub->Draw();
        TMarker* used_marker_proxy = nullptr;
        TMarker* other_marker_proxy = nullptr;
        for (int i = 0; i < n_found; ++i) {
            double peak_x = found_peaks_x[i];
            double peak_y = h_sub->GetBinContent(h_sub->FindBin(peak_x));
            auto m = std::make_unique<TMarker>(peak_x, peak_y, 20);
            if (is_peak_used[i]) {
                m->SetMarkerColor(kRed); m->SetMarkerStyle(29); m->SetMarkerSize(2.0);
                m->Draw("same");
                if (!used_marker_proxy) used_marker_proxy = m.get();
                used_markers.push_back(std::move(m));
            } else {
                m->SetMarkerColor(kGreen); m->SetMarkerStyle(3); m->SetMarkerSize(1.5);
                m->Draw("same");
                if (!other_marker_proxy) other_marker_proxy = m.get();
                other_markers.push_back(std::move(m));
            }
        }

        auto leg_peaks = std::make_unique<TLegend>(0.7, 0.75, 0.9, 0.9);
        if (used_marker_proxy) leg_peaks->AddEntry(used_marker_proxy, "Peaks for Calibration", "p");
        if (other_marker_proxy) leg_peaks->AddEntry(other_marker_proxy, "Other Found Peaks", "p");
        leg_peaks->Draw();
        canvas->SaveAs(TString::Format("%s/pos_%d_2_peak_search.png", output_plot_dir, pos));


        // --- 步骤 4.3: 线性拟合 ---
        if (cal_points.size() < 2) {
            std::cerr << "    Warning: Fewer than 2 peaks found for Pos " << pos << ". Cannot perform fit. Using default parameters." << std::endl;
            k_params[pos] = 1.0; b_params[pos] = 0.0;
            confidence[pos] = cal_points.size();
        }
         else {
            // --- POS % 8 >= 3 拟合逻辑 ---
            if (pos % 8 >= 3 && cal_points.size() >= 3) {
                auto gr_temp = std::make_unique<TGraph>(cal_points.size(), cal_points.data(), ref_points.data());
                TFitResultPtr r_temp = gr_temp->Fit("pol1", "S Q");
                
                if (r_temp.Get()) {
                    double k_full = r_temp->Parameter(1);
                    double b_full = r_temp->Parameter(0);
                    std::vector<std::pair<double, int>> residuals; 
                    for (size_t i = 0; i < cal_points.size(); ++i) {
                        double E_cal = k_full * cal_points[i] + b_full; 
                        double residual = TMath::Abs(ref_points[i] - E_cal); 
                        residuals.push_back({residual, (int)i});
                    }
                    std::sort(residuals.begin(), residuals.end());
                    int idx1 = residuals[0].second;
                    int idx2 = residuals[1].second;
                    
                    std::vector<double> filtered_cal, filtered_ref;
                    filtered_cal.push_back(cal_points[idx1]);
                    filtered_cal.push_back(cal_points[idx2]);
                    filtered_ref.push_back(ref_points[idx1]);
                    filtered_ref.push_back(ref_points[idx2]);
                    
                    cal_points = std::move(filtered_cal);
                    ref_points = std::move(filtered_ref);
                    
                }
            }
            auto gr = std::make_unique<TGraph>(cal_points.size(), cal_points.data(), ref_points.data());
            TFitResultPtr r = gr->Fit("pol1", "S Q");
            if (r.Get()) {
                k_params[pos] = r->Parameter(1);
                b_params[pos] = r->Parameter(0);
                confidence[pos] = cal_points.size();
                canvas->Clear();
                gr->SetTitle(TString::Format("Linear Fit for Pos %d", pos));
                gr->GetXaxis()->SetTitle("Measured Peak Position (channels)");
                gr->GetYaxis()->SetTitle("Reference Energy (keV)");
                gr->SetMarkerStyle(20); gr->SetMarkerSize(1.5);
                gr->Draw("AP");
                auto pt = std::make_unique<TPaveText>(0.2, 0.65, 0.6, 0.85, "NDC");
                pt->SetFillColor(0); pt->SetBorderSize(1); pt->SetTextAlign(12);
                pt->AddText(TString::Format("Fit: E = k * Ch + b"));
                pt->AddText(TString::Format("k = %.5f #pm %.5f", r->Parameter(1), r->ParError(1)));
                pt->AddText(TString::Format("b = %.2f #pm %.2f keV", r->Parameter(0), r->ParError(0)));
                pt->AddText(TString::Format("#chi^{2} / NDF = %.2f / %d", r->Chi2(), r->Ndf()));
                pt->Draw();

                // --- 能量信息显示 ---
                auto pt_E = std::make_unique<TPaveText>(0.65, 0.2, 0.9, 0.45, "NDC");
                pt_E->SetFillColor(0); pt_E->SetBorderSize(1); pt_E->SetTextAlign(12);
                pt_E->AddText("Calibration Points:");
                pt_E->AddText("Ch | Ref E -> Cal E (keV)"); 
                
                double k_fit = r->Parameter(1);
                double b_fit = r->Parameter(0);

                for (size_t i = 0; i < cal_points.size(); ++i) {
                    double E_ref = ref_points[i];
                    double Channel = cal_points[i];
                    double E_cal = k_fit * Channel + b_fit;
                    pt_E->AddText(TString::Format("%.0f | %.2f -> %.2f", Channel, E_ref, E_cal));
                }
                pt_E->Draw();
                canvas->SaveAs(TString::Format("%s/pos_%d_3_linear_fit.png", output_plot_dir, pos));
            } else {
                 std::cerr << "    Error: Fit failed for Pos " << pos << ". Using default parameters." << std::endl;
                 k_params[pos] = 1.0; b_params[pos] = 0.0;
                 confidence[pos] = 0; 
            }
        }
    }
    // --- 步骤 5: 保存参数文件 ---
    std::cout << "\n--> Writing calibration parameters to " << calibration_param_file << "..." << std::endl;
    std::ofstream txt_out(calibration_param_file);
    if (!txt_out.is_open()) {
        std::cerr << "Error: Cannot open output parameter file: " << calibration_param_file << std::endl;
        return FILE_NOT_FOUND;
    }
    txt_out << "# Pos\tk\tb\tconfidence" << std::endl;
    for (int pos = 0; pos < NUM_SSD_POS; ++pos) {
        if(k_params.count(pos)) {
            txt_out << pos << "\t" << k_params[pos] << "\t" << b_params[pos] << "\t" << confidence[pos] << std::endl;
        }
    }
    txt_out.close();
    
    std::cout << "--> Calibration fitting finished successfully!" << std::endl;
    return SUCCESS;
}

// 在文件末尾添加 main 函数

int main(int argc, char* argv[]) {
    // 增加 try-catch 块以捕获配置或运行时的错误
    try {
        ErrorCode status = calibrate_ssd_professional();
        
        if (status == SUCCESS) {
            std::cout << "\nProgram finished successfully." << std::endl;
            return 0; // 返回 0 表示成功
        } else {
            std::cerr << "\nProgram finished with error code: " << status << std::endl;
            return 1; // 返回非 0 表示失败
        }
    } catch (const std::exception& e) {
        std::cerr << "\nRuntime error: " << e.what() << std::endl;
        return 2; // 返回 2 表示异常
    } catch (...) {
        std::cerr << "\nAn unknown error occurred." << std::endl;
        return 3;
    }
}