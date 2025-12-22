/*
 * @brief 使用 RDataFrame 多核并行分析 MWPC 时间，复刻 MWPCTs.c 的逻辑。
 * * 逻辑：
 * 1. 逐个文件循环（因为 MWPCTs.c 的逻辑是 per-file 的）。
 * 2. 在每个文件内部，使用 RDataFrame 多核并行处理事件。
 * 3. 筛选 MWPC 事件：(MWPC_E[0] > 0 || MWPC_E[1] > 0) && (MWPC_Ts[0] > 0 || MWPC_Ts[1] > 0)
 * 4. 定义 "second_bin"：(Long64_t)(MWPC_Ts[0] / 1.e9 + 0.5) - 1
 * 5. 统计三个值：
 * a) TotalEvents：通过筛选的总事件数 (不受阈值影响)。
 * b) Histo：将 "second_bin" 填入一个大直方图（20000 bins）。
 * c) (后处理) SecondsAboveThreshold：'n'，有多少个"秒"的计数 > threshold。
 * d) (后处理) EventsInGoodSeconds：只在(c)的"好秒"中累加的事件总和。
 * 6. 触发 RDataFrame 计算当前文件。
 * 7. 将 RunNum, TotalEvents, EventsInGoodSeconds, SecondsAboveThreshold 写入 txt 文件。
 * 8. 同时保存一个 rate_distribution.root 用于诊断阈值。
 * * 编译命令:
 * g++ analyze_MWPC_time_rdf.cpp $(root-config --cflags --libs) -o analyze_MWPC_time_rdf
 * * 运行命令:
 * ./analyze_MWPC_time_rdf_full
 */

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RVersion.hxx"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include <iostream>
#include <stdio.h> // 用于 FILE, fopen, fprintf, fclose
#include <string>
#include <vector>

// 仅 ROOT>=6.22 才有 RunGraphs
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 22, 0)
#include "ROOT/RDFHelpers.hxx"
#endif

// ----- 可配置参数 -----
// 文件名前缀 (来自用户输入)
const std::string FILENAME_PREFIX = "/home/evalie2/Project/document/273Ds/inter_map/SS032";
// 时间 bin 数量，来自 MWPCTs.c 的 const long NUM = 20000;
const int NUM_BINS = 20000;
// ----- 可配置参数 -----


int main() {
    // 1. 获取用户输入 (使用用户提供的硬编码值)
    int start_file, stop_file, threshold, num_cores;
    start_file = 0;
    stop_file = 180; 
    threshold = 0;   // 默认阈值
    num_cores = 12;   // 默认使用 12 核心

    printf("--- 配置参数 ---\n");
    printf("文件路径: %s\n", FILENAME_PREFIX.c_str());
    printf("文件编号: %d 到 %d\n", start_file, stop_file);
    printf("计数阈值: %d\n", threshold);
    printf("使用核心: %d\n", num_cores);
    printf("------------------\n");

    // 2. 启用 ROOT 隐式多线程
    ROOT::EnableImplicitMT(num_cores);
    std::cout << "已启用 " << num_cores << " 个核心进行处理。" << std::endl;

    // 3. 准备输出文件
    FILE *f_out = fopen("mwpc_counts_per_file.txt", "w");
    if (!f_out) {
        std::cerr << "错误：无法打开输出文件 mwpc_counts_per_file.txt" << std::endl;
        return 1;
    }
    // 写入表头 (包含所有三个统计量)
    fprintf(f_out, "# %-7s %-15s %-20s %-22s\n", "RunNum", "TotalEvents", "EventsInGoodSeconds", "SecondsAboveThreshold");

    // 4. 准备诊断文件
    TFile *f_diag = new TFile("rate_distribution.root", "RECREATE");
    // 这个直方图的X轴是“每秒计数”，Y轴是“有多少秒达到了这个计数”
    // 范围(200, 0, 200)假设你每秒的计数不会超过200，如果超过了请调大
    TH1D *h_rate_distribution = new TH1D("h_rate_distribution",
                                         "Distribution of MWPC Counts per Second;Counts/Second;Number of Seconds",
                                         200000000, 0, 200000000);


    std::cout << "\n开始处理文件，阈值 = " << threshold << " ..." << std::endl;

    // 5. 逐个文件循环
    for (int runnum = start_file; runnum <= stop_file; ++runnum) {
        TString str_map = TString::Format("%s%05d_map.root", FILENAME_PREFIX.c_str(), runnum);
        std::cout << "****************** 正在处理文件: " << str_map << " ******************" << std::endl;

        // 检查文件是否存在
        TFile *f_check = TFile::Open(str_map.Data());
        if (f_check == NULL || f_check->IsZombie()) {
            std::cerr << "错误：无法打开文件 " << str_map << "，跳过此文件。" << std::endl;
            if (f_check) f_check->Close();
            continue;
        }
        f_check->Close();

        // 为当前文件创建 RDataFrame
        ROOT::RDataFrame df("tr_map", str_map.Data());

        // 定义筛选条件 (来自 MWPCTs.c)
        auto df_filtered = df.Filter("(MWPC_E[0] > 0. || MWPC_E[1] > 0.) && (MWPC_Ts[0] > 0. || MWPC_Ts[1] > 0)", "MWPC valid event");

        // 定义 "second_bin" (来自 MWPCTs.c 的 t-1)
        auto df_binned = df_filtered.Define("second_bin", "(Long64_t)(MWPC_Ts[0] / 1.e9 + 0.5) - 1");

        // 注册计算任务
        // 任务1: 总有效事件数 (不受阈值影响)
        auto total_events = df_filtered.Count();

        // 任务2: "秒"谱直方图 (模拟 mwpctime 数组)
        // (修复 TH1DModel 错误：使用 {...} 语法)
        auto h_time_seconds = df_binned.Histo1D(
            {"h_time_seconds", "MWPC Events per Second", NUM_BINS, -0.5, (double)NUM_BINS - 0.5},
            "second_bin");


        // 触发当前文件的计算
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 22, 0)
        // std::cout << "  (使用 RunGraphs 触发计算...)" << std::endl;
        ROOT::RDF::RunGraphs({total_events, h_time_seconds});
#else
        // std::cout << "  (使用 GetValue 触发计算...)" << std::endl;
        // 在旧版 ROOT 中，我们必须手动触发
        total_events.GetValue();   // 触发 Count()
        h_time_seconds.GetValue(); // 触发 Histo1D()
#endif

        // --- 后处理 ---
        // 获取计算结果
        long long total_events_val = *total_events; // (A) 总事件数
        TH1D h_result = *h_time_seconds;            // mwpctime 数组

        // 遍历直方图，统计超阈值的 bin 和 事件
        long long seconds_above_threshold = 0;  // (B) "好秒" 总数
        long long events_in_good_seconds = 0;   // (C) "好秒" 内的事件总和

        for (int i = 1; i <= h_result.GetNbinsX(); ++i) { // TH1 bin 编号从 1 开始
            double rate_in_this_second = h_result.GetBinContent(i);
            
            // 填充诊断直方图
            h_rate_distribution->Fill(rate_in_this_second);

            if (rate_in_this_second > threshold) {
                seconds_above_threshold++; // (B) "好秒"+1
                events_in_good_seconds += (long long)rate_in_this_second; // (C) 累加“好”事件
            }
        }

        // 输出结果到控制台和文件
        printf("  > 文件: %s\n", str_map.Data());
        printf("    - (A) 总有效事件数 (不受阈值): %lld\n", total_events_val);
        printf("    - (C) 超阈值秒内的事件总和: %lld\n", events_in_good_seconds);
        printf("    - (B) 计数 > %d 的秒数 (n): %lld\n", threshold, seconds_above_threshold);

        fprintf(f_out, "  %-7d %-15lld %-20lld %-22lld\n", runnum, total_events_val, events_in_good_seconds, seconds_above_threshold);
    }

    // 6. 清理
    fclose(f_out);
    
    // 保存诊断文件
    f_diag->Write();
    f_diag->Close();

    std::cout << "\n✅ 处理完成！结果已保存至 mwpc_counts_per_file.txt" << std::endl;
    std::cout << "✅ 诊断直方图已保存至 rate_distribution.root" << std::endl;

    return 0;
}