
//g++ analyze_rdf_2D.cpp $(root-config --cflags --libs) -o analyze_rdf_2D
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TColor.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "ROOT/RVersion.hxx"
#include "ROOT/RDF/HistoModels.hxx"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#endif
#include <vector>
#include <iostream>

void analyze_rdf_2D() {
    gStyle->SetPalette(kViridis);
    ROOT::EnableImplicitMT(12);

    TChain chain("tr_map");
    chain.Add("/home/evalie2/Project/document/273Ds/inter_map/SS032001*_map.root");
    ROOT::RDataFrame df(chain);

    auto n_entries = df.Count();
    std::cout << "使用 RDataFrame (12核心) 准备处理 " << *n_entries << " 个事件..." << std::endl;

    // -------------------- DSSD (这部分不变) --------------------
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
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_SSD_E_vs_tdiff;
    for (int i = 0; i < 5; ++i) {
        int ch_min = i * 8;
        int ch_max = i * 8 + 7;
        // ✅ 修正: 将 RVec<int> 修改为 RVec<double> 并添加 (int) 转换
        auto df_ssd_ch = df_ssd.Filter([=](const ROOT::RVec<double>& ch_vec) {
            return (int)ch_vec[0] >= ch_min && (int)ch_vec[0] <= ch_max;
        }, {"SSD_Ch"}, Form("SSD Ch %d-%d", ch_min, ch_max));

        auto model_ssd =ROOT::RDF::TH2DModel(
            Form("SSD_E_vs_tdiff_%d_%d", ch_min, ch_max),
            Form("DSSD E vs Time Diff (SSD Ch %d-%d);T_DSSD - T_SSD (ns);DSSDX_E (keV)", ch_min, ch_max),
            60, -300, 300,
            200, 0, 15000);

        h_SSD_E_vs_tdiff.push_back(
            df_ssd_ch.Define("dssd_ssd_tdiff", "(double)DSSDX_Ts[0] - (double)SSD_Ts[0]")
                     .Histo2D(model_ssd, "dssd_ssd_tdiff", "DSSDX_E"));
    }

    // -------------------- MWPC --------------------
    auto df_mwpc = df_dssd_cut.Filter("MWPC_mul == 2 && MWPC_Ch[0] == 0 && MWPC_Ch[1] == 1", "Good MWPC event");
    
    auto model_mwpc0 = ROOT::RDF::TH2DModel("MWPC_E_vs_tdiff_0", "DSSD E vs Time Diff (MWPC Ch 0);T_DSSD - T_MWPC (ns);DSSDX_E (keV)", 100, -500, 500, 4000, 0, 40000);
    auto h_MWPC_E_vs_tdiff_0 = df_mwpc.Define("dssd_mwpc0_tdiff", "(double)DSSDX_Ts[0] - (double)MWPC_Ts[0]")
                                      .Histo2D(model_mwpc0, "dssd_mwpc0_tdiff", "DSSDX_E");

    auto model_mwpc1 = ROOT::RDF::TH2DModel("MWPC_E_vs_tdiff_1", "DSSD E vs Time Diff (MWPC Ch 1);T_DSSD - T_MWPC (ns);DSSDX_E (keV)", 100, -500, 500, 4000, 0, 40000);
    auto h_MWPC_E_vs_tdiff_1 = df_mwpc.Define("dssd_mwpc1_tdiff", "(double)DSSDX_Ts[0] - (double)MWPC_Ts[1]")
                                      .Histo2D(model_mwpc1, "dssd_mwpc1_tdiff", "DSSDX_E");
    auto model_mwpc = ROOT::RDF::TH2DModel("MWPC_E_vs_tdiff", "DSSD E vs Time Diff (MWPC[0]-MWPC[1]);T_MWPC[0] - T_MWPC[1] (ns);DSSDX_E (keV)", 100, -500, 500, 4000, 0, 40000);                                  
    auto h_MWPC_E_vs_tdiff   = df_mwpc.Define("mwpc0_mwpc1_tdiff", "(double)MWPC_Ts[0] - (double)MWPC_Ts[1]")
                                      .Histo2D(model_mwpc, "mwpc0_mwpc1_tdiff", "DSSDX_E");

    // -------------------- VETO --------------------
    auto df_veto = df_dssd_cut.Filter("Veto_mul == 1", "Veto multiplicity == 1");
    std::vector<ROOT::RDF::RResultPtr<TH2D>> h_Veto_E_vs_tdiff;
    for (int i = 0; i < 3; ++i) {
        // ✅ 修正: 将 RVec<int> 修改为 RVec<double> 并添加 (int) 转换
        auto df_veto_ch = df_veto.Filter([=](const ROOT::RVec<double>& ch_vec) {
            return (int)ch_vec[0] == i;
        }, {"Veto_Ch"}, Form("Veto Ch %d", i));
        
        auto model_veto = ROOT::RDF::TH2DModel(
            Form("Veto_E_vs_tdiff_%d", i),
            Form("DSSD E vs Time Diff (Veto Ch %d);T_DSSD - T_Veto (ns);DSSDX_E (keV)", i),
            100, -500, 500,
            4000, 0, 40000);

        h_Veto_E_vs_tdiff.push_back(
            df_veto_ch.Define("dssd_veto_tdiff", "(double)DSSDX_Ts[0] - (double)Veto_Ts[0]")
                      .Histo2D(model_veto, "dssd_veto_tdiff", "DSSDX_E"));
    }

    // -------------------- 执行带进度条 --------------------
    std::cout << "\n所有计算任务已注册，准备执行 (带进度条)..." << std::endl;

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
    std::cout << "检测到 ROOT 版本 >= 6.22，使用 ROOT::RDF::RunGraphs() 显示进度条。" << std::endl;
    std::vector<ROOT::RDF::RResultHandle> tasks = {
        h_DSSD_Ediff, h_DSSD_tdiff,
        h_MWPC_E_vs_tdiff_0, h_MWPC_E_vs_tdiff_1, h_MWPC_E_vs_tdiff
    };
    for (auto& h : h_SSD_E_vs_tdiff) tasks.push_back(h);
    for (auto& h : h_Veto_E_vs_tdiff) tasks.push_back(h);
    
    ROOT::RDF::Experimental::AddProgressBar(df);
    ROOT::RDF::RunGraphs(tasks);
#else
    std::cout << "检测到 ROOT 版本 < 6.22，不支持 RunGraphs。将通过触发最后一个直方图来运行所有任务。" << std::endl;
    h_Veto_E_vs_tdiff.back()->GetValue();
#endif

    std::cout << "\n所有任务执行完毕，开始保存输出文件..." << std::endl;

    // -------------------- 保存输出 --------------------
    TFile *outputFile = new TFile("coin_window_2D_rdf.root", "RECREATE");
    h_DSSD_Ediff->Write();
    h_DSSD_tdiff->Write();
    for (auto& h : h_SSD_E_vs_tdiff) h->Write();
    h_MWPC_E_vs_tdiff_0->Write();
    h_MWPC_E_vs_tdiff_1->Write();
    h_MWPC_E_vs_tdiff->Write();
    for (auto& h : h_Veto_E_vs_tdiff) h->Write();
    outputFile->Close();

    std::cout << "✅ 处理完成！所有直方图已保存至 coin_window_2D_rdf.root" << std::endl;
}

int main() {
    analyze_rdf_2D();
    return 0;
}