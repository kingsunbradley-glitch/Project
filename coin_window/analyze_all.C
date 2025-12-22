#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <vector>

// ===================================================================
//
//  高效整合分析宏 (内存版本): analyze_from_memory.C
//
//  功能:
//  - 直接使用已经存在于ROOT内存中的 TChain 或 TTree 对象进行分析
//  - 通过一次事件循环完成所有直方图的填充，极大提升运行效率
//
// ===================================================================

void analyze_all() {
    // --- 步骤 1: 从ROOT内存中获取数据对象 ---
    // ###############################################################
    // #  重要：这里将直接寻找名为 "tr_map" 的对象                  #
    // #  它可以是你手动创建的 TChain，也可以是打开文件得到的 TTree   #
    // ###############################################################
    
    const char* object_name = "tr_map"; // 假设你的 TChain/TTree 对象名为 "tr_map"
    
    // TTree 是 TChain 的基类，所以这个指针可以同时指向 TTree 或 TChain
    TTree *tr_map = (TTree*)gROOT->FindObject(object_name);

    // 检查是否成功获取对象
    if (!tr_map) {
        std::cerr << "错误: 在当前ROOT内存中未找到名为 '" << object_name << "' 的对象。" << std::endl;
        std::cerr << "请确保你已经成功创建了 TChain 或打开了 TFile 里的 TTree，"
                  << "并且它的对象名就是 '" << object_name << "'。" << std::endl;
        return;
    }
    std::cout << "成功从内存中获取对象 '" << object_name << "'，总事件数: " << tr_map->GetEntries() << std::endl;

    // --- 后续所有代码与之前的版本完全相同 ---

    // --- 步骤 2: 提前定义所有需要的直方图 ---
    std::cout << "正在初始化所有直方图..." << std::endl;
    TH1F *h_DSSD_Ediff = new TH1F("DSSD_Ediff", "DSSD X-Y Energy Difference;E_x - E_y (keV);Counts", 200, -500, 500);
    TH1F *h_DSSD_tdiff = new TH1F("DSSD_tdiff", "DSSD X-Y Time Difference;T_x - T_y (ns);Counts", 100, -500, 500);
    std::vector<TH1F*> h_SSD_tdiff;
    for (int i = 0; i < 5; ++i) {
        h_SSD_tdiff.push_back(new TH1F(Form("SSD_tdiff_%d_%d", i * 8, i * 8 + 7), Form("DSSD-SSD Time Difference (Ch %d-%d);T_DSSD - T_SSD (ns);Counts", i * 8, i * 8 + 7), 200, -2000, 2000));
    }
    TH1F *h_MWPC_tdiff_0 = new TH1F("MWPC_tdiff_0", "DSSD-MWPC Time Difference (Ch 0);T_DSSD - T_MWPC (ns);Counts", 200, -2000, 2000);
    TH1F *h_MWPC_tdiff_1 = new TH1F("MWPC_tdiff_1", "DSSD-MWPC Time Difference (Ch 1);T_DSSD - T_MWPC (ns);Counts", 200, -2000, 2000);
    std::vector<TH1F*> h_Veto_tdiff;
    for (int i = 0; i < 3; ++i) {
        h_Veto_tdiff.push_back(new TH1F(Form("Veto_tdiff_%d", i), Form("DSSD-Veto Time Difference (Ch %d);T_DSSD - T_Veto (ns);Counts", i), 100, -200, 200));
    }

    // --- 步骤 3: 设置Branch Address以高效读取数据 ---
    const int MAX_HITS = 1024;
    UShort_t DSSDX_mul, DSSDY_mul, SSD_mul;
    UShort_t MWPC_mul, Veto_mul;
    Double_t DSSDX_E[MAX_HITS], DSSDY_E[MAX_HITS];
    ULong64_t DSSDX_Ts[MAX_HITS], DSSDY_Ts[MAX_HITS];
    Double_t SSD_Ch[MAX_HITS];
    ULong64_t SSD_Ts[MAX_HITS];
    Double_t MWPC_Ch[MAX_HITS];
    ULong64_t MWPC_Ts[MAX_HITS];
    Double_t Veto_Ch[MAX_HITS];
    ULong64_t Veto_Ts[MAX_HITS];
    tr_map->SetBranchAddress("DSSDX_mul", &DSSDX_mul);
    tr_map->SetBranchAddress("DSSDY_mul", &DSSDY_mul);
    tr_map->SetBranchAddress("DSSDX_E", DSSDX_E);
    tr_map->SetBranchAddress("DSSDY_E", DSSDY_E);
    tr_map->SetBranchAddress("DSSDX_Ts", DSSDX_Ts);
    tr_map->SetBranchAddress("DSSDY_Ts", DSSDY_Ts);
    tr_map->SetBranchAddress("SSD_mul", &SSD_mul);
    tr_map->SetBranchAddress("SSD_Ch", SSD_Ch);
    tr_map->SetBranchAddress("SSD_Ts", SSD_Ts);
    tr_map->SetBranchAddress("MWPC_mul", &MWPC_mul);
    tr_map->SetBranchAddress("MWPC_Ch", MWPC_Ch);
    tr_map->SetBranchAddress("MWPC_Ts", MWPC_Ts);
    tr_map->SetBranchAddress("Veto_mul", &Veto_mul);
    tr_map->SetBranchAddress("Veto_Ch", Veto_Ch);
    tr_map->SetBranchAddress("Veto_Ts", Veto_Ts);

    // --- 步骤 4: 单次事件循环，填充所有直方图 ---
    Long64_t n_entries = tr_map->GetEntries();
    std::cout << "开始处理 " << n_entries << " 个事件..." << std::endl;
    for (Long64_t i = 0; i < n_entries; ++i) {
        tr_map->GetEntry(i);
        if (i % 1000000 == 0 && i > 0) {
            std::cout << "已处理 " << i << " / " << n_entries << " (" << (100.0 * i / n_entries) << "%)" << std::endl;
        }
        if (DSSDX_mul == 1 && DSSDY_mul == 1) {
            h_DSSD_Ediff->Fill(DSSDX_E[0] - DSSDY_E[0]);
            if (fabs(DSSDX_E[0] - DSSDY_E[0]) < 500) {
                h_DSSD_tdiff->Fill(DSSDX_Ts[0] - DSSDY_Ts[0]);
            }
            if (SSD_mul == 1) {
                int index = SSD_Ch[0] / 8;
                if (index >= 0 && index < h_SSD_tdiff.size()) {
                    h_SSD_tdiff[index]->Fill(DSSDX_Ts[0] - SSD_Ts[0]);
                }
            }
            if (MWPC_mul == 2 && MWPC_Ch[0] == 0 && MWPC_Ch[1] == 1) {
                h_MWPC_tdiff_0->Fill(DSSDX_Ts[0] - MWPC_Ts[0]);
                h_MWPC_tdiff_1->Fill(DSSDX_Ts[0] - MWPC_Ts[1]);
            }
            if (Veto_mul == 1) {
                int veto_ch = Veto_Ch[0];
                if (veto_ch >= 0 && veto_ch < h_Veto_tdiff.size()) {
                    h_Veto_tdiff[veto_ch]->Fill(DSSDX_Ts[0] - Veto_Ts[0]);
                }
            }
        }
    }
    std::cout << "事件处理完成！" << std::endl;

    // --- 步骤 5: 绘制所有直方图 ---
    std::cout << "正在绘制所有直方图..." << std::endl;
    TCanvas *c_dssd = new TCanvas("c_dssd", "DSSD Information", 1200, 600);
    c_dssd->Divide(2, 1);
    c_dssd->cd(1); h_DSSD_Ediff->Draw("h");
    c_dssd->cd(2); h_DSSD_tdiff->Draw("h");
    TCanvas *c_ssd = new TCanvas("c_ssd", "DSSD-SSD Time Difference", 1200, 800);
    c_ssd->Divide(3, 2);
    for (size_t i = 0; i < h_SSD_tdiff.size(); ++i) {
        c_ssd->cd(i + 1); h_SSD_tdiff[i]->Draw("h");
    }
    TCanvas *c_mwpc = new TCanvas("c_mwpc", "DSSD-MWPC Time Difference", 1200, 600);
    c_mwpc->Divide(2, 1);
    c_mwpc->cd(1); h_MWPC_tdiff_0->Draw("h");
    c_mwpc->cd(2); h_MWPC_tdiff_1->Draw("h");
    TCanvas *c_veto = new TCanvas("c_veto", "DSSD-Veto Time Difference", 1200, 500);
    c_veto->Divide(3, 1);
    for (size_t i = 0; i < h_Veto_tdiff.size(); ++i) {
        c_veto->cd(i + 1); h_Veto_tdiff[i]->Draw("h");
    }
    std::cout << "所有图像已生成。" << std::endl;
}