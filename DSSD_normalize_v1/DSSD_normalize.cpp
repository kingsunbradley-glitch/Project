#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <memory> 
#include <stdexcept> 
#include <algorithm>
#include <thread>

// RDataFrame 相关的头文件
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDFHelpers.hxx" 

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"

#include "config.h"

using namespace ROOT::VecOps;

// 归一化拟合的核心函数
ErrorCode DSSD_Normalize_Process() {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gSystem->mkdir(output_plot_dir, kTRUE);

    ROOT::EnableImplicitMT(num_threads);
    std::cout << "\n--> Starting DSSD Normalization Process (Step 1, using up to " 
              << num_threads << " threads)..." << std::endl;
    
    // 【优化】: 直接从预先筛选过的 sorted_root_file 读取数据
    std::cout << "    Input file (pre-filtered): " << sorted_root_file << std::endl;
    
    ROOT::RDataFrame df(tree_name, sorted_root_file);
    auto n_entries = df.Count();
    std::cout << "    Total entries in pre-filtered TTree: " << *n_entries << std::endl;
    if (*n_entries == 0) {
        std::cerr << "Error: No entries found in the pre-filtered TTree." << std::endl;
        std::cerr << "       Hint: Please run the Sorter program first." << std::endl;
        return INSUFFICIENT_DATA;
    }
    
    ROOT::RDF::Experimental::AddProgressBar(df);

    // 【优化】: 核心筛选步骤已由 Sorter 程序完成，此处不再需要全局 Filter
    // 仅定义后续分析所需的列
    auto df_defined = df.Define("ChX", "(int)DSSDX_Ch[0]")
                        .Define("EX", "DSSDX_E[0]")
                        .Define("ChY", "(int)DSSDY_Ch[0]")
                        .Define("EY", "DSSDY_E[0]");

    // --- 【修正部分 1】: 整合诊断文件和自动寻找逻辑 ---
    TString diagnostic_root_file = TString::Format("%s/%s_Normalize_Diagnostic.root", output_plot_dir, EXPERIMENT_ID);
    
    // 【修正部分 2】: 使用正确的 unique_ptr 构造方式
    std::unique_ptr<TFile> outFile(TFile::Open(diagnostic_root_file, "RECREATE"));
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot create diagnostic ROOT file: " << diagnostic_root_file << std::endl;
        return FILE_NOT_FOUND;
    }

    std::cout << "\n--> Finding optimal reference channels with highest statistics..." << std::endl;
    // 使用 df_defined 而不是 df_filtered
    auto h_x_hits = df_defined.Histo1D<int>({"h_x_hits", "X-side Hit Distribution;X Channel;Counts", NUM_DSSDX_POS, -0.5, NUM_DSSDX_POS - 0.5}, "ChX");
    auto h_y_hits = df_defined.Histo1D<int>({"h_y_hits", "Y-side Hit Distribution;Y Channel;Counts", NUM_DSSDY_POS, -0.5, NUM_DSSDY_POS - 0.5}, "ChY");

    int best_X_Normalize_Ch = h_x_hits->GetMaximumBin() - 1;
    int best_Y_Normalize_Ch = h_y_hits->GetMaximumBin() - 1;

    if (best_X_Normalize_Ch < 0 || h_x_hits->GetBinContent(best_X_Normalize_Ch + 1) < MIN_FIT_ENTRIES) {
        std::cerr << "Error: Insufficient data to determine a reliable X reference channel." << std::endl;
        return INSUFFICIENT_DATA;
    }
    if (best_Y_Normalize_Ch < 0 || h_y_hits->GetBinContent(best_Y_Normalize_Ch + 1) < MIN_FIT_ENTRIES) {
        std::cerr << "Error: Insufficient data to determine a reliable Y reference channel." << std::endl;
        return INSUFFICIENT_DATA;
    }

    std::cout << "    Optimal reference for Y-side normalization: X_Ch = " << best_X_Normalize_Ch << std::endl;
    std::cout << "    Optimal reference for X-side normalization: Y_Ch = " << best_Y_Normalize_Ch << std::endl;

    // --- 归一化操作 ---
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_x_results, h_y_results;

    for (int ch = 0; ch < NUM_DSSDX_POS; ++ch) {
        TString hname = TString::Format("h2_X%d", ch);
        TString htitle = TString::Format("X-side Ch %d vs Y-side Ch %d;DSSDX_E (X-axis);DSSDY_E (Y-axis)", ch, best_Y_Normalize_Ch);
        // 使用 df_defined，Filter仅用于此处的局部通道选择
        ROOT::RDF::RResultPtr<TH2D> h2_ptr = df_defined.Filter(TString::Format("ChX == %d && ChY == %d", ch, best_Y_Normalize_Ch).Data())
                             .Histo2D<TH2D>({hname, htitle, 1000, MIN_ENERGY_THRESHOLD, MAX_ENERGY_THRESHOLD, 1000, MIN_ENERGY_THRESHOLD, MAX_ENERGY_THRESHOLD}, "EX", "EY");
        h_x_results.push_back(std::move(h2_ptr));
    }

    for (int ch = 0; ch < NUM_DSSDY_POS; ++ch) {
        TString hname = TString::Format("h2_Y%d", ch);
        TString htitle = TString::Format("Y-side Ch %d vs X-side Ch %d;DSSDY_E (X-axis);DSSDX_E (Y-axis)", ch, best_X_Normalize_Ch);
        ROOT::RDF::RResultPtr<TH2D> h2_ptr = df_defined.Filter(TString::Format("ChY == %d && ChX == %d", ch, best_X_Normalize_Ch).Data())
                             .Histo2D<TH2D>({hname, htitle, 1000, MIN_ENERGY_THRESHOLD, MAX_ENERGY_THRESHOLD, 1000, MIN_ENERGY_THRESHOLD, MAX_ENERGY_THRESHOLD}, "EY", "EX");
        h_y_results.push_back(std::move(h2_ptr));
    }

    std::cout << "\n    Executing RDataFrame tasks (plotting 2D histograms)..." << std::endl;
    // 触发所有计算
    for (auto& h_ptr : h_x_results) h_ptr.GetValue();
    for (auto& h_ptr : h_y_results) h_ptr.GetValue();
    std::cout << "    RDataFrame execution finished." << std::endl;

    // --- 单核拟合和参数保存 --- (这部分代码无需改动)
    std::map<int, double> k_x_params, b_x_params;
    std::map<int, double> k_y_params, b_y_params;
    
    auto canvas = std::make_unique<TCanvas>("c_normalize", "Normalization Canvas", 1200, 600);
    canvas->Divide(2, 1);
    auto f1_linear = std::make_unique<TF1>("f1_linear", "pol1", MIN_ENERGY_THRESHOLD, MAX_ENERGY_THRESHOLD);
    
    std::cout << "\n--> Starting single-core linear fitting for normalization..." << std::endl;

    std::cout << "    Fitting X-side (relative to ChY=" << best_Y_Normalize_Ch << "):" << std::endl;
    for (int ch = 0; ch < NUM_DSSDX_POS; ++ch) {
        TH2D* h2 = h_x_results[ch].GetPtr();
        if (h2->GetEntries() < MIN_FIT_ENTRIES) {
            k_x_params[ch] = 1.0; b_x_params[ch] = 0.0;
            continue;
        }

        TProfile* prof = h2->ProfileX(TString::Format("prof_X%d", ch));
        canvas->cd(1); h2->Draw("COLZ");
        canvas->cd(2); prof->Draw();
        
        f1_linear->SetParameters(0.0, 1.0);
        TFitResultPtr r = prof->Fit(f1_linear.get(), "SQR");
        
        if (r.Get() && r->IsValid()) {
            k_x_params[ch] = r->Parameter(1);
            b_x_params[ch] = r->Parameter(0);
        } else {
            k_x_params[ch] = 1.0; b_x_params[ch] = 0.0;
        }
        
        canvas->SaveAs(TString::Format("%s/h2_X%d_fit.png", output_plot_dir, ch));
        outFile->cd();
        h2->Write();
        prof->Write();
        delete prof;
    }

    std::cout << "    Fitting Y-side (relative to ChX=" << best_X_Normalize_Ch << "):" << std::endl;
    for (int ch = 0; ch < NUM_DSSDY_POS; ++ch) {
        TH2D* h2 = h_y_results[ch].GetPtr();
        if (h2->GetEntries() < MIN_FIT_ENTRIES) {
            k_y_params[ch] = 1.0; b_y_params[ch] = 0.0;
            continue;
        }

        TProfile* prof = h2->ProfileX(TString::Format("prof_Y%d", ch));
        canvas->cd(1); h2->Draw("COLZ");
        canvas->cd(2); prof->Draw();

        f1_linear->SetParameters(0.0, 1.0);
        TFitResultPtr r = prof->Fit(f1_linear.get(), "SQR");
        
        if (r.Get() && r->IsValid()) {
            k_y_params[ch] = r->Parameter(1);
            b_y_params[ch] = r->Parameter(0);
        } else {
            k_y_params[ch] = 1.0; b_y_params[ch] = 0.0;
        }
        
        canvas->SaveAs(TString::Format("%s/h2_Y%d_fit.png", output_plot_dir, ch));
        outFile->cd();
        h2->Write();
        prof->Write();
        delete prof;
    }
    
    outFile->Close();
    
    std::cout << "\n--> Writing normalization parameters to " << normalize_param_file << "..." << std::endl;
    std::ofstream txt_out(normalize_param_file);
    if (!txt_out.is_open()) {
        std::cerr << "Error: Cannot open output parameter file: " << normalize_param_file << std::endl;
        return FILE_NOT_FOUND;
    }
    txt_out << "# DSSD Normalization Parameters: Channel k b" << std::endl;
    txt_out << "# X-side (0-" << NUM_DSSDX_POS - 1 << ") relative to ChY=" << best_Y_Normalize_Ch << std::endl;
    for (int ch = 0; ch < NUM_DSSDX_POS; ++ch) {
        txt_out << ch << "\t" << k_x_params[ch] << "\t" << b_x_params[ch] << std::endl;
    }
    txt_out << "# Y-side (0-" << NUM_DSSDY_POS - 1 << ") relative to ChX=" << best_X_Normalize_Ch << std::endl;
    for (int ch = 0; ch < NUM_DSSDY_POS; ++ch) {
        txt_out << ch << "\t" << k_y_params[ch] << "\t" << b_y_params[ch] << std::endl;
    }
    txt_out.close();
    
    std::cout << "--> DSSD Normalization finished successfully! Diagnostic ROOT: " << diagnostic_root_file << std::endl;
    return SUCCESS;
}

int main() {
    try {
        return DSSD_Normalize_Process();
    } catch (const std::exception& e) {
        std::cerr << "\nRuntime error: " << e.what() << std::endl;
        return 2; 
    } catch (...) {
        std::cerr << "\nAn unknown error occurred." << std::endl;
        return 3;
    }
}