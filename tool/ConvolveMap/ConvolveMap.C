#include "TH2.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include <iostream>

/**
 * @brief 对二维直方图进行高斯卷积以模拟探测器分辨率
 * * @param h_in 输入的原始 TH2 直方图 (Truth)
 * @param sigma_x X轴方向的分辨率 (物理单位，非Bin数)
 * @param sigma_y Y轴方向的分辨率 (物理单位，非Bin数)
 * @return TH2* 返回卷积后的新直方图 (Reco)，调用者负责 delete
 */
TH2* ApplyDetectorResolution(TH2* h_in, double sigma_x, double sigma_y) {
    if (!h_in) {
        std::cerr << "Error: Input histogram is null!" << std::endl;
        return nullptr;
    }

    // 1. 克隆直方图结构，并重命名，清空内容
    TH2* h_out = (TH2*)h_in->Clone(Form("%s_smeared", h_in->GetName()));
    h_out->SetTitle(Form("%s (Convolved #sigma_{x}=%.2f, #sigma_{y}=%.2f)", h_in->GetTitle(), sigma_x, sigma_y));
    h_out->Reset(); // 清空内容，我们将重新填充

    // 获取轴的指针
    TAxis* xAxis = h_in->GetXaxis();
    TAxis* yAxis = h_in->GetYaxis();
    int nBinsX = h_in->GetNbinsX();
    int nBinsY = h_in->GetNbinsY();

    // 2. 遍历原始直方图的每一个 Bin (源 Bin)
    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            
            double content = h_in->GetBinContent(i, j);
            if (content == 0) continue; // 跳过空 Bin 以节省时间

            double x_center = xAxis->GetBinCenter(i);
            double y_center = yAxis->GetBinCenter(j);

            // 3. 优化：只计算源 Bin 周围 +/- 3 sigma 范围内的目标 Bin
            // 这样避免了对整个直方图进行循环，极大提高速度
            int min_k = xAxis->FindBin(x_center - 3 * sigma_x);
            int max_k = xAxis->FindBin(x_center + 3 * sigma_x);
            int min_l = yAxis->FindBin(y_center - 3 * sigma_y);
            int max_l = yAxis->FindBin(y_center + 3 * sigma_y);

            // 边界检查
            if (min_k < 1) min_k = 1;
            if (max_k > nBinsX) max_k = nBinsX;
            if (min_l < 1) min_l = 1;
            if (max_l > nBinsY) max_l = nBinsY;

            // 临时存储这个点源弥散后的总权重，用于归一化
            double weight_sum = 0.0;
            
            // 第一次循环：计算权重并存储（如果追求极致速度可以合并，但分开写逻辑更清晰）
            // 这里我们直接在一个步骤中计算并填入，依靠高斯函数的归一化系数
            
            // 高斯归一化系数 A = 1 / (2 * pi * sx * sy)
            // 注意：因为我们是离散Bin，最好的方式是计算完所有权重后强制归一化到1，
            // 确保能量（Counts）守恒。
            
            // 这种局部卷积需要一个小的 buffer 或者直接计算
            // 为了保持代码简洁，我们使用"直接投射法"，并假设 bin width 远小于 sigma
            // 如果 bin width 很大，需要对高斯函数在 bin 内积分。这里假设 bin 中心点近似。
            
            std::vector<std::tuple<int, int, double>> targets;
            
            for (int k = min_k; k <= max_k; ++k) {
                for (int l = min_l; l <= max_l; ++l) {
                    double x_dest = xAxis->GetBinCenter(k);
                    double y_dest = yAxis->GetBinCenter(l);

                    // 计算高斯权重
                    double dx = x_dest - x_center;
                    double dy = y_dest - y_center;
                    double g_val = TMath::Gaus(dx, 0, sigma_x) * TMath::Gaus(dy, 0, sigma_y);
                    
                    targets.emplace_back(k, l, g_val);
                    weight_sum += g_val;
                }
            }

            // 4. 将源 Bin 的内容按权重分配给目标 Bins
            if (weight_sum > 0) {
                for (const auto& t : targets) {
                    int binX = std::get<0>(t);
                    int binY = std::get<1>(t);
                    double weight = std::get<2>(t);
                    
                    // 新内容 = 旧内容 * (该位置的权重 / 总权重)
                    // 这样保证了 sum(新内容) == 旧内容，即粒子数守恒
                    double add_val = content * (weight / weight_sum);
                    
                    // 累加到目标直方图
                    double current_val = h_out->GetBinContent(binX, binY);
                    h_out->SetBinContent(binX, binY, current_val + add_val);
                }
            }
        }
    }

    return h_out;
}

// ==========================================
// 测试函数：生成一个 dummy 直方图并运行上述算法
// ==========================================
void ConvolveMap() {
    // 设置绘图风格
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    // 1. 创建一个示例直方图 (Truth)
    TH2D* h_truth = new TH2D("h_truth", "Truth Distribution;X;Y", 100, -5, 5, 100, -5, 5);
    
    // 填充一些数据（例如两个清晰的圆环或点）
    TRandom3 rnd(0);
    for (int i = 0; i < 10000; ++i) {
        // 模拟一个中心在 (-2, -2) 的点源
        h_truth->Fill(rnd.Gaus(-2, 0.2), rnd.Gaus(-2, 0.2));
        // 模拟一个中心在 (2, 2) 的环形
        double theta = rnd.Uniform(0, TMath::TwoPi());
        double r = rnd.Gaus(1.5, 0.1);
        h_truth->Fill(2 + r*cos(theta), 2 + r*sin(theta));
    }

    // 2. 设置分辨率参数 (例如 X方向分辨差，Y方向分辨好)
    double sigmaX = 0.5; 
    double sigmaY = 0.2;

    // 3. 调用我们的卷积函数
    TH2* h_reco = ApplyDetectorResolution(h_truth, sigmaX, sigmaY);

    // 4. 绘图对比
    TCanvas* c1 = new TCanvas("c1", "Convolution Demo", 1200, 500);
    c1->Divide(2, 1);

    c1->cd(1);
    h_truth->Draw("COLZ");
    
    c1->cd(2);
    if (h_reco) h_reco->Draw("COLZ");

    std::cout << "Done! Original Integral: " << h_truth->Integral() 
              << ", Smeared Integral: " << h_reco->Integral() << std::endl;
}