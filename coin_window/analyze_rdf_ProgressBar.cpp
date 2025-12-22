//g++ analyze_rdf_ProgressBar.cpp $(root-config --cflags --libs) -o analyze_rdf_ProgressBar
//using version 
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TChain.h"
#include "TFile.h"
#include "ROOT/RVersion.hxx"
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
#include "ROOT/RDFHelpers.hxx" // 仅 ROOT>=6.22 才有 RunGraphs
#endif
#include <vector>
#include <iostream>

void analyze_rdf() {
    ROOT::EnableImplicitMT(12);

    TChain chain("tr_map");
    chain.Add("/home/evalie2/Project/document/273Ds/inter_map/SS032001*_map.root");
    ROOT::RDataFrame df(chain);

    auto n_entries = df.Count();
    std::cout << "使用 RDataFrame (12核心) 准备处理 " << *n_entries << " 个事件..." << std::endl;

    // -------------------- DSSD --------------------
    auto df_filtered = df.Filter("DSSDX_mul == 1 && DSSDY_mul == 1", "DSSD XY multiplicity == 1");
    auto df_dssd_vars = df_filtered.Define("Ediff", "DSSDX_E[0] - DSSDY_E[0]")
                                   .Define("Tdiff", "(double)DSSDX_Ts[0] - (double)DSSDY_Ts[0]");
    auto df_dssd_cut = df_dssd_vars.Filter("abs(Ediff) < 500", "DSSD Ediff cut");

    auto h_DSSD_Ediff = df_dssd_vars.Histo1D<double>(
        {"DSSD_Ediff", "DSSD X-Y Energy Difference;E_x - E_y (keV);Counts", 60, -300, 300}, "Ediff");
    auto h_DSSD_tdiff = df_dssd_cut.Histo1D<double>(
        {"DSSD_tdiff", "DSSD X-Y Time Difference;T_x - T_y (ns);Counts", 60, -300, 300}, "Tdiff");

    // -------------------- SSD --------------------
    auto df_ssd = df_dssd_cut.Filter("SSD_mul == 1", "SSD multiplicity == 1");
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_SSD_tdiff;
    for (int i = 0; i < 6; ++i) {
        int ch_min = i * 8;
        int ch_max = i * 8 + 7;
        auto df_ssd_ch = df_ssd.Filter([=](const ROOT::RVec<double>& ch_vec) {
            return ch_vec[0] >= ch_min && ch_vec[0] <= ch_max;
        }, {"SSD_Ch"}, Form("SSD Ch %d-%d", ch_min, ch_max));
        h_SSD_tdiff.push_back(
            df_ssd_ch.Define("dssd_ssd_tdiff", "(double)DSSDX_Ts[0] - (double)SSD_Ts[0]")
                     .Histo1D<double>({Form("SSD_tdiff_%d_%d", ch_min, ch_max),
                                       Form("DSSD-SSD Time Difference (Ch %d-%d);T_DSSD - T_SSD (ns);Counts", ch_min, ch_max),
                                       60, -300, 300},
                                      "dssd_ssd_tdiff"));
    }

    // -------------------- MWPC --------------------
    //auto df_mwpc = df_dssd_cut.Filter("MWPC_mul == 2 && MWPC_Ch[0] == 0 && MWPC_Ch[1] == 1 && DSSDX_E[0] > 4000", "Good MWPC event");
    auto df_mwpc = df_dssd_cut.Filter("MWPC_mul == 2 && MWPC_Ch[0] == 0 && MWPC_Ch[1] == 1 ", "Good MWPC event");
    auto h_MWPC_tdiff_0     = df_mwpc.Define("dssd_mwpc0_tdiff", "(double)DSSDX_Ts[0] - (double)MWPC_Ts[0]")
                                 .Histo1D<double>({"MWPC_tdiff_0", "DSSD-MWPC Time Difference (Ch 0);T_DSSD - T_MWPC (ns);Counts", 60, -300, 300}, "dssd_mwpc0_tdiff");
    auto h_MWPC_tdiff_1     = df_mwpc.Define("dssd_mwpc1_tdiff", "(double)DSSDX_Ts[0] - (double)MWPC_Ts[1]")
                                 .Histo1D<double>({"MWPC_tdiff_1", "DSSD-MWPC Time Difference (Ch 1);T_DSSD - T_MWPC (ns);Counts", 60, -300, 300}, "dssd_mwpc1_tdiff");
    
    auto h_MWPC_tdiff_MWPC  = df_mwpc.Define("mwpc0_mwpc1_tdiff", "((double)MWPC_Ts[0] - (double)MWPC_Ts[1])") // 1. 修正列名 2. 补全括号
                                        .Histo1D<double>({"MWPC_tdiff_MWPC", // 3. 修正直方图唯一名称
                                                        "MWPC_Ts[0]-MWPC_Ts[1] Time Difference;T_MWPC_Ts[0] - T_MWPC[1] (ns);Counts", 60, -300, 300}, 
                                                       "mwpc0_mwpc1_tdiff"); // 4. 匹配 Define 的列名


    
    // -------------------- VETO --------------------
    auto df_veto = df_dssd_cut.Filter("Veto_mul == 1", "Veto multiplicity == 1");
    std::vector<ROOT::RDF::RResultPtr<TH1D>> h_Veto_tdiff;
    for (int i = 0; i < 3; ++i) {
        auto df_veto_ch = df_veto.Filter([=](const ROOT::RVec<double>& ch_vec) {
            return ch_vec[0] == i;
        }, {"Veto_Ch"}, Form("Veto Ch %d", i));
        h_Veto_tdiff.push_back(
            df_veto_ch.Define("dssd_veto_tdiff", "(double)DSSDX_Ts[0] - (double)Veto_Ts[0]")
                      .Histo1D<double>({Form("Veto_tdiff_%d", i),
                                        Form("DSSD-Veto Time Difference (Ch %d);T_DSSD - T_Veto (ns);Counts", i),
                                        60, -300, 300},
                                       "dssd_veto_tdiff"));
    }

    // -------------------- 执行带进度条 --------------------
    std::cout << "\n所有计算任务已注册，准备执行 (带进度条)..." << std::endl;

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
    // ✅ ROOT >= 6.22：使用 RunGraphs() 并行运行所有任务
    std::cout << "检测到 ROOT 版本 >= 6.22，使用 ROOT::RDF::RunGraphs() 显示进度条。" << std::endl;
    std::vector<ROOT::RDF::RResultHandle> tasks = {
        //n_entries,  // 触发 Count() 显示总事件数
        h_DSSD_Ediff, h_DSSD_tdiff, h_MWPC_tdiff_0, h_MWPC_tdiff_1,h_MWPC_tdiff_MWPC
    };
    for (auto& h : h_SSD_tdiff) tasks.push_back(h);
    for (auto& h : h_Veto_tdiff) tasks.push_back(h);
    
    // 🔑 关键：触发事件循环并显示进度
    ROOT::RDF::Experimental::AddProgressBar(df);  // 如果可用
    ROOT::RDF::RunGraphs(tasks);
#else
    // ⚙️ 旧 ROOT：手动触发 Count() 显示进度条（同时触发所有任务）
    std::cout << "检测到 ROOT 版本 < 6.22，使用 Count() 显示进度条。" << std::endl;
    auto progress = df.Count();
    progress.GetValue();
#endif

    std::cout << "\n所有任务执行完毕，开始保存输出文件..." << std::endl;

    // -------------------- 保存输出 --------------------
    TFile *outputFile = new TFile("coin_window_with_energy_rdf.root", "RECREATE");
    h_DSSD_Ediff->Write();
    h_DSSD_tdiff->Write();
    for (auto& h : h_SSD_tdiff) h->Write();
    h_MWPC_tdiff_0->Write();
    h_MWPC_tdiff_1->Write();
    h_MWPC_tdiff_MWPC->Write();
    for (auto& h : h_Veto_tdiff) h->Write();
    outputFile->Close();

    std::cout << "✅ 处理完成！所有直方图已保存至 coin_window_with_energy_rdf.root" << std::endl;
}

int main() {
    analyze_rdf();
    return 0;
}
