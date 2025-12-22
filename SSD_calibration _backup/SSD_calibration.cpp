// 编译命令:
// g++ SSD_calibration.cpp `root-config --cflags --libs` -lSpectrum -o SSD_calibration
// 基于自动寻峰的SSD内刻度程序
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <sstream>
#include <memory> // For std::unique_ptr
#include <stdexcept> // For std::runtime_error

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
// 所有路径、参数和结构体定义都从此文件获得
#include "config.h"

// 函数：验证 config.h 中的参数是否合理
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


ErrorCode calibrate_ssd_professional() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // --- 步骤 0: 验证配置参数 ---
    try {
      // validate_config();
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return CONFIG_ERROR;
    }

    gSystem->mkdir(output_plot_dir, kTRUE);

    // 使用智能指针管理ROOT对象，避免内存泄漏
    auto chain = std::make_unique<TChain>(tree_name);
    chain->Add(input_root_file);

    if (chain->GetEntries() == 0) {
        std::cerr << "Error: No entries found in TChain. Check input file path and tree name." << std::endl;
        return FILE_NOT_FOUND;
    }

    std::map<int, double> k_params, b_params;
    std::map<int, int> confidence;
    auto canvas = std::make_unique<TCanvas>("canvas", "Calibration Canvas", 1200, 900);

    std::cout << "\n--> Starting professional calibration process..." << std::endl;

    for (int pos = 0; pos < NUM_SSD_POS; ++pos) {
        std::cout << "    Processing Pos " << pos << "..." << std::endl;

        auto h_sum_E = std::make_unique<TH1D>(TString::Format("h_sum_E_pos%d", pos), TString::Format("Energy Spectrum for Pos %d;Energy Sum (channels);Counts", pos), 700, 3000, 10000);
        TCut cut = TCut(TString::Format("SSD_E > 0 && DSSD_E > 0 && TMath::Nint(SSD_Pos) == %d", pos));
        
        chain->Draw(TString::Format("DSSD_E + SSD_E >> %s", h_sum_E->GetName()), cut, "gooff");

        if (h_sum_E->GetEntries() < MIN_HIST_ENTRIES) {
            std::cerr << "    Warning: Insufficient data for Pos " << pos << " (Entries: " << h_sum_E->GetEntries() << "). Skipping." << std::endl;
            k_params[pos] = 1.0; b_params[pos] = 0.0; confidence[pos] = -1; // 标记为数据不足
            continue;
        }

        // --- 步骤 1: 背景扣除 ---
        auto s = std::make_unique<TSpectrum>();
        std::unique_ptr<TH1D> h_sub; 

        if (enable_background_subtraction) {
            canvas->Clear();

            TH1 *h_bg_raw = s->Background(h_sum_E.get(), background_iterations, "same");
            auto h_bg = std::unique_ptr<TH1>(h_bg_raw); // 接管裸指针

            
            h_sum_E->Draw();
            h_bg->SetLineColor(kRed);  // 设置背景线颜色
            h_bg->Draw("same");  // 绘制背景线
            auto leg_bg = std::make_unique<TLegend>(0.7, 0.75, 0.9, 0.9);
            leg_bg->AddEntry(h_sum_E.get(), "Original Spectrum", "l");
            leg_bg->AddEntry(h_bg.get(), "Estimated Background", "l");
            leg_bg->Draw();
            canvas->SaveAs(TString::Format("%s/pos_%d_1_background.png", output_plot_dir, pos));
            
            h_sub.reset((TH1D*)h_sum_E->Clone(TString::Format("h_sub_pos%d", pos)));
            h_sub->Add(h_bg.get(), -1);
        } else {
            h_sub.reset((TH1D*)h_sum_E->Clone(TString::Format("h_sub_pos%d", pos)));
        }

        // --- 步骤 2: 自动寻峰 ---
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

        // --- 步骤 3: 线性拟合 ---
        if (cal_points.size() < 2) {
            std::cerr << "    Warning: Fewer than 2 peaks found for Pos " << pos << ". Cannot perform fit. Using default parameters." << std::endl;
            k_params[pos] = 1.0; b_params[pos] = 0.0;
            confidence[pos] = cal_points.size();
        }
         else {
            // --- ADDED LOGIC FOR POS % 8 >= 3: 选择两个最接近整体拟合趋势的点 ---
            if (pos % 8 >= 3 && cal_points.size() >= 3) {
                std::cout << "    Special handling: Pos " << pos << " (%8 >= 3) has " << cal_points.size() << " points. Selecting 2 points closest to the overall fit line (based on residuals)." << std::endl;
                
                // 1. 拟合所有点，确定一个“理想”的线性关系
                auto gr_temp = std::make_unique<TGraph>(cal_points.size(), cal_points.data(), ref_points.data());
                TFitResultPtr r_temp = gr_temp->Fit("pol1", "S Q");
                
                if (r_temp.Get()) {
                    double k_full = r_temp->Parameter(1);
                    double b_full = r_temp->Parameter(0);

                    // 2. 计算每个点到这条拟合线的残差（Residuals）
                    std::vector<std::pair<double, int>> residuals; // pair<residual, original_index>
                    for (size_t i = 0; i < cal_points.size(); ++i) {
                        double E_cal = k_full * cal_points[i] + b_full; // 拟合线给出的能量
                        double residual = TMath::Abs(ref_points[i] - E_cal); // 实际能量与拟合能量的差
                        residuals.push_back({residual, (int)i});
                    }

                    // 3. 排序，选择残差最小的两个点
                    std::sort(residuals.begin(), residuals.end());
                    
                    // 选择前两个（残差最小）
                    int idx1 = residuals[0].second;
                    int idx2 = residuals[1].second;
                    
                    std::vector<double> filtered_cal, filtered_ref;
                    filtered_cal.push_back(cal_points[idx1]);
                    filtered_cal.push_back(cal_points[idx2]);
                    filtered_ref.push_back(ref_points[idx1]);
                    filtered_ref.push_back(ref_points[idx2]);
                    
                    cal_points = std::move(filtered_cal);
                    ref_points = std::move(filtered_ref);
                    
                    std::cout << "    Selected points (Ch, Ref E): (" << cal_points[0] << ", " << ref_points[0] << ") and (" << cal_points[1] << ", " << ref_points[1] << ")" << std::endl;
                    std::cout << "    Their residuals from the full-point fit: " << residuals[0].first << " and " << residuals[1].first << std::endl;

                } else {
                    std::cerr << "    Warning: Full-point fit failed for Pos " << pos << ". Using all " << cal_points.size() << " points." << std::endl;
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

                // --- ADDED LOGIC: 显示所有用于拟合的点的信息 ---

                auto pt_E = std::make_unique<TPaveText>(0.65, 0.2, 0.9, 0.45, "NDC");
                pt_E->SetFillColor(0); pt_E->SetBorderSize(1); pt_E->SetTextAlign(12);
                pt_E->AddText("Calibration Points:");
                pt_E->AddText("Ch | Ref E -> Cal E (keV)"); // 增加表头
                
                double k_fit = r->Parameter(1);
                double b_fit = r->Parameter(0);

                // 循环显示所有用于拟合的点的信息
                for (size_t i = 0; i < cal_points.size(); ++i) {
                    double E_ref = ref_points[i];
                    double Channel = cal_points[i];
                    double E_cal = k_fit * Channel + b_fit;
                    
                    // 格式：通道数(Ref E) -> Cal E keV
                    pt_E->AddText(TString::Format("%.0f | %.2f -> %.2f", Channel, E_ref, E_cal));
                }
                pt_E->Draw();
                // --- END ADDED LOGIC ---

                canvas->SaveAs(TString::Format("%s/pos_%d_3_linear_fit.png", output_plot_dir, pos));
                canvas->SaveAs(TString::Format("%s/pos_%d_3_linear_fit.png", output_plot_dir, pos));
            } else {
                 std::cerr << "    Error: Fit failed for Pos " << pos << ". Using default parameters." << std::endl;
                 k_params[pos] = 1.0; b_params[pos] = 0.0;
                 confidence[pos] = 0; // 标记为拟合失败
            }
        }
    }

    // --- 步骤 4: 保存参数文件 ---
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

int main() {
    return calibrate_ssd_professional();
}
