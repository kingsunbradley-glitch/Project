#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TChain.h"
#include "TFile.h"
#include <vector>
#include <iostream>

void analyze_rdf() {
    ROOT::EnableImplicitMT(8);

    TChain chain("tr_map");
    chain.Add("/home/evalie2/Project/document/273Ds/inter_map/SS0320001*_map.root");
    ROOT::RDataFrame df(chain);

    auto n_entries = df.Count();
    std::cout << "使用 RDataFrame (8核心) 准备处理 " << *n_entries << " 个事件..." << std::endl;

    // --- 所有分析逻辑和直方图预定都保持不变 ---
    auto df_filtered = df.Filter("DSSDX_mul == 1 && DSSDY_mul == 1", "DSSD XY multiplicity == 1");
    auto df_dssd_vars = df_filtered.Define("Ediff", "DSSDX_E[0] - DSSDY_E[0]")
                                   .Define("Tdiff", "(double)DSSDX_Ts[0] - (double)DSSDY_Ts[0]");
    auto df_dssd_cut = df_dssd_vars.Filter("abs(Ediff) < 500", "DSSD Ediff cut");

    auto h_DSSD_Ediff = df_dssd_vars.Histo1D<double>({"DSSD_Ediff", "DSSD X-Y Energy Difference;E_x - E_y (keV);Counts", 200, -300, 300}, "Ediff");
    auto h_DSSD_tdiff = df_dssd_cut.Histo1D<double>({"DSSD_tdiff", "DSSD X-Y Time Difference;T_x - T_y (ns);Counts", 200, -300, 300}, "Tdiff");

    auto df_ssd = df_dssd_cut.Filter("SSD_mul == 1", "SSD multiplicity == 1");
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_SSD_tdiff;
    for (int i = 0; i < 5; ++i) {
        int ch_min = i * 8;
        int ch_max = i * 8 + 7;
        auto df_ssd_ch = df_ssd.Filter([=](const ROOT::RVec<double>& ch_vec) {
            return ch_vec[0] >= ch_min && ch_vec[0] <= ch_max;
        }, {"SSD_Ch"}, Form("SSD Ch %d-%d", ch_min, ch_max));
        h_SSD_tdiff.push_back(
            df_ssd_ch.Define("dssd_ssd_tdiff", "(double)DSSDX_Ts[0] - (double)SSD_Ts[0]")
                     .Histo1D<double>({Form("SSD_tdiff_%d_%d", ch_min, ch_max), Form("DSSD-SSD Time Difference (Ch %d-%d);T_DSSD - T_SSD (ns);Counts", ch_min, ch_max), 200, -300, 300}, "dssd_ssd_tdiff")
        );
    }
    
    auto df_mwpc = df_dssd_cut.Filter("MWPC_mul == 2 && MWPC_Ch[0] == 0 && MWPC_Ch[1] == 1", "Good MWPC event");
    auto h_MWPC_tdiff_0 = df_mwpc.Define("dssd_mwpc0_tdiff", "(double)DSSDX_Ts[0] - (double)MWPC_Ts[0]").Histo1D<double>({"MWPC_tdiff_0", "DSSD-MWPC Time Difference (Ch 0);T_DSSD - T_MWPC (ns);Counts", 200, -2000, 2000}, "dssd_mwpc0_tdiff");
    auto h_MWPC_tdiff_1 = df_mwpc.Define("dssd_mwpc1_tdiff", "(double)DSSDX_Ts[0] - (double)MWPC_Ts[1]").Histo1D<double>({"MWPC_tdiff_1", "DSSD-MWPC Time Difference (Ch 1);T_DSSD - T_MWPC (ns);Counts", 200, -2000, 2000}, "dssd_mwpc1_tdiff");

    auto df_veto = df_dssd_cut.Filter("Veto_mul == 1", "Veto multiplicity == 1");
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_Veto_tdiff;
    for (int i = 0; i < 3; ++i) {
        auto df_veto_ch = df_veto.Filter([=](const ROOT::RVec<double>& ch_vec) {
            return ch_vec[0] == i;
        }, {"Veto_Ch"}, Form("Veto Ch %d", i));
        h_Veto_tdiff.push_back(
            df_veto_ch.Define("dssd_veto_tdiff", "(double)DSSDX_Ts[0] - (double)Veto_Ts[0]")
                      .Histo1D<double>({Form("Veto_tdiff_%d", i), Form("DSSD-Veto Time Difference (Ch %d);T_DSSD - T_Veto (ns);Counts", i), 50, -200, 200}, "dssd_veto_tdiff")
        );
    }

    // 【修正点】使用最简单、最标准的方式触发带进度条的事件循环
    // -----------------------------------------------------------------
    std::cout << "所有计算任务已预定，开始执行事件循环..." << std::endl;
    // 只需访问任意一个结果，即可触发所有计算。GetValue() 是一个很明确的选择。
    // 这时，进度条会自动出现在终端。
    h_DSSD_Ediff.GetValue();
    // -----------------------------------------------------------------

    // --- 事件循环结束后，所有直方图都已填好，现在保存它们 ---
    TFile *outputFile = new TFile("histograms_rdf.root", "RECREATE");
    
    h_DSSD_Ediff->Write();
    h_DSSD_tdiff->Write();
    for(auto& h : h_SSD_tdiff) h->Write();
    h_MWPC_tdiff_0->Write();
    h_MWPC_tdiff_1->Write();
    for(auto& h : h_Veto_tdiff) h->Write();
    
    outputFile->Close();
    
    std::cout << "处理完成！所有直方图已保存至 histograms_rdf.root" << std::endl;
}

int main() {
    analyze_rdf();
    return 0;
}