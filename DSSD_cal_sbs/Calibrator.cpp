#include "Calibrator.h"
#include "Config.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <algorithm>

using namespace std;
using namespace ROOT; // For RDataFrame, RVec

Calibrator::Calibrator() {
    // 启用多线程 (必须在构建 RDF 之前调用)
    if (NUM_THREADS > 0) {
        ROOT::EnableImplicitMT(NUM_THREADS);
        cout << "--> Multi-threading enabled with " << NUM_THREADS << " threads." << endl;
    } else {
        ROOT::EnableImplicitMT(); // 自动检测
        cout << "--> Multi-threading enabled with auto-detected threads." << endl;
    }

    h_resolution = new TH1F("h_resolution", "FWHM Resolution vs Channel;Channel;FWHM (keV)", TOTAL_CH, 0, TOTAL_CH);
}

Calibrator::~Calibrator() {
    // 内存清理
}

void Calibrator::InitialiseHistograms() {
    // 初始化主直方图容器
    vec_h_raw.resize(TOTAL_CH);
    for(int i=0; i<TOTAL_CH; ++i) {
        vec_h_raw[i] = new TH1F(TString::Format("h_raw_%d", i), 
                                TString::Format("Raw Spectrum Ch %d", i), 
                                HIST_BINS, HIST_MIN, HIST_MAX);
        // 关键：由于我们要手动合并，这里先不让它们关联到任何目录，避免多线程冲突
        vec_h_raw[i]->SetDirectory(0);
    }
}

// ================= RDataFrame 核心填充逻辑 =================
void Calibrator::FillHistogramsRDF(TChain* chain) {
    cout << "--> Setting up RDataFrame..." << endl;
    ROOT::RDataFrame df(*chain);

    // 1. 定义筛选条件 (Filter)
    // 筛选条件：DSSDX_mul==1 && DSSDY_mul==1 && MWPC_mul==0 && Veto_mul==0 && SSD_mul==0
    auto df_filter = df.Filter(" DSSDX_mul==1 && DSSDY_mul==1 && MWPC_mul == 0 && Veto_mul == 0 && SSD_mul == 0", "Event Cut");

    // 2. 准备线程局部存储 (Thread Local Storage)
    // RDataFrame 会并行运行，为了避免锁竞争，我们为每个 Slot (线程) 创建一套独立的直方图
    // 获取实际使用的 Slot 数量 (通常等于线程数)
    unsigned int nSlots = df.GetNSlots();
    cout << "--> Allocating histograms for " << nSlots << " slots..." << endl;

    // 结构：[SlotID][ChannelID] -> TH1F*
    std::vector<std::vector<TH1F*>> thread_hists(nSlots);
    
    for(unsigned int s = 0; s < nSlots; ++s) {
        thread_hists[s].resize(TOTAL_CH);
        for(int i = 0; i < TOTAL_CH; ++i) {
            // 创建线程私有的直方图，名字需要不同或不注册到gDirectory
            TString tmpName = TString::Format("h_raw_%d_slot_%u", i, s);
            thread_hists[s][i] = new TH1F(tmpName, "", HIST_BINS, HIST_MIN, HIST_MAX);
            thread_hists[s][i]->SetDirectory(0); // 关键：完全脱离 ROOT 全局管理，纯内存对象
        }
    }

    cout << "--> Starting Event Loop (Processing all files)..." << endl;

    // 3. 执行并行循环 (ForeachSlot)
    // 我们直接读取需要的列：Energy 和 Channel，以及 YH_mul (用于判断 YH 是否有效)
    // 注意：RDataFrame 读取数组列类型为 RVec<double>
    df_filter.ForeachSlot(
        [&thread_hists](unsigned int slot, 
                        const ROOT::VecOps::RVec<double>& x_e, const ROOT::VecOps::RVec<double>& x_ch,
                        const ROOT::VecOps::RVec<double>& y_e, const ROOT::VecOps::RVec<double>& y_ch,
                        const ROOT::VecOps::RVec<double>& yh_e, const ROOT::VecOps::RVec<double>& yh_ch,
                        UShort_t yh_mul) {
            
            // 此时已经是筛选过的事件：X_mul=1, Y_mul=1
            // 也就是 x_ch.size() >= 1, y_ch.size() >= 1

            // --- Fill X (0-127) ---
            int chX = (int)x_ch[0];
            if (chX >= DSSD_X_CH_LOW && chX <= DSSD_X_CH_HIGH) {
                thread_hists[slot][chX]->Fill(x_e[0]);
            }

            // --- Fill Y (128-175) ---
            int chY = (int)y_ch[0];
            // 映射：0-47 -> 128-175
            if (chY >= 0 && chY <= (DSSD_Y_CH_HIGH - DSSD_Y_CH_LOW)) {
                thread_hists[slot][DSSD_Y_CH_LOW + chY]->Fill(y_e[0]);
            }

            // --- Fill YH (176-223) ---
            // YH 需要额外判断 mul >= 1
            if (yh_mul >= 1 && yh_ch.size() > 0) {
                int chYH = (int)yh_ch[0];
                if (chYH >= 0 && chYH <= (DSSD_YH_CH_HIGH - DSSD_YH_CH_LOW)) {
                    thread_hists[slot][DSSD_YH_CH_LOW + chYH]->Fill(yh_e[0]);
                }
            }
        },
        {"DSSDX_E", "DSSDX_Ch", "DSSDY_E", "DSSDY_Ch", "DSSDYH_E", "DSSDYH_Ch", "DSSDYH_mul"}
    );

    // 4. 合并结果 (Reduce)
    cout << "--> Merging results from " << nSlots << " threads..." << endl;
    for(int i = 0; i < TOTAL_CH; ++i) {
        // 将所有 Slot 的直方图加到 vec_h_raw[i] 上
        for(unsigned int s = 0; s < nSlots; ++s) {
            vec_h_raw[i]->Add(thread_hists[s][i]);
            delete thread_hists[s][i]; // 合并后删除临时直方图，节省内存
        }
    }
    
    // 打印处理的事件总数报告
    auto count = df_filter.Count();
    cout << "--> Total valid events processed: " << *count << endl;
}


// ================= 以下是刻度逻辑 (基本未变) =================

TH1F* Calibrator::BackgroundSubtract(TH1F* h_in) {
    if(!ENABLE_BKG_SUB) return (TH1F*)h_in->Clone();
    
    TH1F* h_out = (TH1F*)h_in->Clone();
    TSpectrum *s = new TSpectrum();
    TH1 *h_bkg = s->Background(h_out, BKG_ITERATIONS, "same");
    h_out->Add(h_bkg, -1); 
    delete s;
    return h_out; 
}

bool Calibrator::FindSinglePeak(TH1F* h, double center, double window, double &mean, double &sigma) {
    h->GetXaxis()->SetRangeUser(center - window, center + window);
    int maxBin = h->GetMaximumBin();
    double maxE = h->GetBinCenter(maxBin);
    double maxContent = h->GetBinContent(maxBin);

    if(maxContent < 5) return false; 

    TF1 *fGaus = new TF1("fGaus", "gaus", maxE - window, maxE + window);
    fGaus->SetParameters(maxContent, maxE, 5.0); 
    fGaus->SetParLimits(1, center - window, center + window); 

    int status = h->Fit(fGaus, "RQN"); 
    
    if(status == 0) {
        mean = fGaus->GetParameter(1);
        sigma = fGaus->GetParameter(2);
        delete fGaus;
        return true;
    }
    delete fGaus;
    return false;
}

void Calibrator::ProcessChannel(int chID, ofstream &ofs) {
    TH1F* h_process = vec_h_raw[chID];
    
    // 1. 本底扣除
    TH1F* h_sub = BackgroundSubtract(h_process);

    // 2. 寻找参考峰
    vector<double> x_meas; 
    vector<double> y_true; 
    
    for(const auto& peak : REF_PEAKS) {
        double mean = 0, sigma = 0;
        if(FindSinglePeak(h_sub, peak.energy, peak.window, mean, sigma)) {
            x_meas.push_back(mean);
            y_true.push_back(peak.energy);
        }
    }

    // 3. 线性拟合
    double k = 1.0, b = 0.0, chi2 = 0.0;
    bool fit_success = false;

    if(x_meas.size() >= 2) {
        TGraph *gr = new TGraph(x_meas.size(), &x_meas[0], &y_true[0]);
        TF1 *pol1 = new TF1("pol1", "pol1");
        gr->Fit(pol1, "Q");
        b = pol1->GetParameter(0); 
        k = pol1->GetParameter(1); 
        chi2 = pol1->GetChisquare();
        fit_success = true;
        delete gr;
        delete pol1;
    } else {
        k = 1.0; b = 0.0; chi2 = 999.0;
    }

    // 4. 生成刻度后直方图 (诊断用)
    TString name = TString::Format("h%d", chID); 
    TH1F* h_cal = new TH1F(name, TString::Format("Calibrated Ch %d;Energy (keV);Counts", chID), 
                           HIST_BINS, k*HIST_MIN+b, k*HIST_MAX+b); 
    
    for(int i=1; i<=h_process->GetNbinsX(); ++i) {
        double raw_center = h_process->GetBinCenter(i);
        double content = h_process->GetBinContent(i);
        if(content > 0)
            h_cal->Fill(k * raw_center + b, content);
    }

    // 5. 诊断峰 FWHM
    double fwhm_max = 0.0;
    if(fit_success) {
        for(const auto& diag : DIAG_PEAKS) {
            double mean = 0, sigma = 0;
            if(FindSinglePeak(h_cal, diag.energy, diag.window, mean, sigma)) {
                double fwhm = 2.355 * sigma;
                if(fwhm > fwhm_max) fwhm_max = fwhm;
                
                TLine *l = new TLine(mean, 0, mean, h_cal->GetMaximum()*0.8);
                l->SetLineColor(kRed); l->SetLineStyle(2);
                h_cal->GetListOfFunctions()->Add(l);
            }
            TLine *l_ref = new TLine(diag.energy, 0, diag.energy, h_cal->GetMaximum());
            l_ref->SetLineColor(kGreen+2);
            h_cal->GetListOfFunctions()->Add(l_ref);
        }
    }

    // 6. 写入 Dat
    ofs << TString::Format("%-4d\t%10.4f\t%10.6f\t%10.4f\t%10.4f", chID, b, k, fwhm_max, chi2).Data() << endl;

    // 7. 分辨率图
    h_resolution->SetBinContent(chID+1, fwhm_max);

    // 存入容器
    vec_h_cal.push_back(h_cal);
    
    delete h_sub;
}

void Calibrator::Run() {
    InitialiseHistograms();

    // === 构建 TChain (核心修改) ===
    TChain* chain = new TChain(TREE_NAME);
    cout << "--> Adding files to TChain from run " << RUN_START << " to " << RUN_END << "..." << endl;
    
    int fileCount = 0;
    for(int i = RUN_START; i <= RUN_END; ++i) {
        TString filename = TString::Format(INPUT_DIR_PATTERN, i);
        // 使用 AccessPathName 检查文件是否存在 (返回 false 表示存在)
        if(gSystem->AccessPathName(filename) == false) {
            chain->Add(filename);
            fileCount++;
        }
    }
    
    if(fileCount == 0) {
        cerr << "Error: No files found matching pattern!" << endl;
        return;
    }
    cout << "--> Added " << fileCount << " files. Total TChain entries (estimated): " << chain->GetEntries() << endl;

    // 1. 使用 RDataFrame 填充直方图
    FillHistogramsRDF(chain);
    
    // TChain 任务结束，删除
    delete chain; 

    // 2. 准备输出
    f_out_root = new TFile(OUTPUT_ROOT_FILE, "RECREATE");
    ofstream ofs(OUTPUT_DAT_FILE);
    
    cout << "--> Performing Calibration Fits..." << endl;
    for(int i=0; i<TOTAL_CH; ++i) {
        ProcessChannel(i, ofs);
        if(i % 10 == 0) cout << "Calibrating channel " << i << " / " << TOTAL_CH << "\r" << flush;
    }
    cout << endl;

    // 3. 写入 ROOT
    f_out_root->cd();
    for(auto h : vec_h_cal) {
        if(h) h->Write();
    }
    h_resolution->Write();
    
    f_out_root->Close();
    ofs.close();

    cout << "Done! Results saved to " << OUTPUT_DAT_FILE << " and " << OUTPUT_ROOT_FILE << endl;
}