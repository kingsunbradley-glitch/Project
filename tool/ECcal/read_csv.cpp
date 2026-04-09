#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cmath>

#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"   // for R² calculation

// 去除字符串首尾的所有空白字符
std::string trim(const std::string& str) {
    auto start = std::find_if_not(str.begin(), str.end(), [](unsigned char ch) {
        return std::isspace(ch);
    });
    auto end = std::find_if_not(str.rbegin(), str.rend(), [](unsigned char ch) {
        return std::isspace(ch);
    }).base();
    return (start < end) ? std::string(start, end) : "";
}

// 去掉 UTF-8 BOM
std::string removeBOM(const std::string& str) {
    if (str.size() >= 3 && str[0] == '\xEF' && str[1] == '\xBB' && str[2] == '\xBF') {
        return str.substr(3);
    }
    return str;
}

// 计算拟合优度 R^2
double ComputeR2(TGraph* gr, TF1* f) {
    double chi2 = f->GetChisquare();
    double sumY = 0.0;
    double sumY2 = 0.0;
    int n = gr->GetN();
    for (int i = 0; i < n; ++i) {
        double y = gr->GetPointY(i);
        sumY += y;
        sumY2 += y * y;
    }
    double tss = sumY2 - sumY * sumY / n;   // 总平方和
    if (tss <= 0) return 1.0;
    return 1.0 - chi2 / tss;
}

void read_csv() {
    // ================== 1. 读取 CSV 文件 ==================
    std::string filename = "ECdata.csv";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "can not open file: " << filename << std::endl;
        return;
    }

    std::vector<double> Q, V, Z;
    std::string line;
    int lineNum = 0;

    while (std::getline(file, line)) {
        lineNum++;
        if (line.find_first_not_of(" ,\t\r") == std::string::npos) continue;

        if (lineNum == 1) {
            line = removeBOM(line);
        }

        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> fields;

        while (std::getline(ss, token, ',')) {
            token = trim(token);
            fields.push_back(token);
        }

        if (fields.size() < 3) {
            std::cout << "Warning: line " << lineNum << " insufficient fields, skipped." << std::endl;
            continue;
        }

        try {
            double q = std::stod(fields[0]);
            double v = std::stod(fields[1]);
            double z = std::stod(fields[2]);
            Q.push_back(q);
            V.push_back(v);
            Z.push_back(z);
        } catch (const std::exception& e) {
            std::cout << "Warning: line " << lineNum << " conversion failed (" << e.what() << "), skipped." << std::endl;
            continue;
        }
    }
    file.close();

    std::cout << "Successfully read " << Q.size() << " rows of data." << std::endl;
    if (Q.size() < 2) {
        std::cout << "Error: insufficient data points." << std::endl;
        return;
    }

    // ================== 2. 初始线性拟合 (V -> Q) ==================
    TGraph* grLin = new TGraph(Q.size(), &V[0], &Q[0]);
    double vMin = *std::min_element(V.begin(), V.end());
    double vMax = *std::max_element(V.begin(), V.end());
    TF1* fLin = new TF1("linear", "[0] + [1]*x", vMin, vMax);
    grLin->Fit(fLin, "Q");
    double k0 = fLin->GetParameter(1);   // 斜率
    double b0 = fLin->GetParameter(0);   // 截距
    std::cout << "\nInitial linear fit: Q = " << k0 << " * V + " << b0 << std::endl;

    delete grLin;
    delete fLin;

    // ================== 3. 低分辨率sin拟合 ==================
    const int N = 20;
    double kmin = k0 - 0.1;
    double kmax = k0 + 0.1;
    double bmin = b0 - 1.0;
    double bmax = b0 + 1.0;

    std::vector<double> k_vals(N), b_vals(N);
    for (int i = 0; i < N; ++i) {
        k_vals[i] = kmin + (kmax - kmin) * i / (N - 1);
        b_vals[i] = bmin + (bmax - bmin) * i / (N - 1);
    }

    std::vector<std::vector<double>> r2(N, std::vector<double>(N, 0.0));

    const double pi16 = TMath::Pi() / 16.0;   // π/16

    for (int i = 0; i < N; ++i) {
        double k = k_vals[i];
        for (int j = 0; j < N; ++j) {
            double b = b_vals[j];

            // 计算 y_fit, y_cor, x_fit
            std::vector<double> y_cor(Q.size()), x_fit(Q.size());
            for (size_t m = 0; m < Q.size(); ++m) {
                double y_fit = k * V[m] + b;
                y_cor[m] = Q[m] - y_fit;
                x_fit[m] = Z[m] - y_fit;
            }

            // 正弦拟合: a * sin(π/16 * (x - phase))
            TGraph* grSin = new TGraph(Q.size(), &x_fit[0], &y_cor[0]);
            double xmin = *std::min_element(x_fit.begin(), x_fit.end());
            double xmax = *std::max_element(x_fit.begin(), x_fit.end());
            if (xmax - xmin < 1e-9) {
                r2[i][j] = 0.0;
                delete grSin;
                continue;
            }

            TF1* fSin = new TF1("sinFit", "[0]*sin(pi/16*(x-[1]))", xmin, xmax);
            fSin->SetParameter(0, 0.7);
            fSin->SetParameter(1, 0.03);
            fSin->SetParName(0, "a");
            fSin->SetParName(1, "phase");

            grSin->Fit(fSin, "Q");   // 静默拟合
            r2[i][j] = ComputeR2(grSin, fSin);

            delete grSin;
            delete fSin;
        }
    }

    // 寻找最大 R² 对应的索引
    double maxR2 = 0.0;
    int i_max = 0, j_max = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (r2[i][j] > maxR2) {
                maxR2 = r2[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }

    double k_best1 = k_vals[i_max];
    double b_best1 = b_vals[j_max];
    std::cout << "\nCoarse scan: best (k, b) = (" << k_best1 << ", " << b_best1 << "), R^2= " << maxR2 << std::endl;

    // ================== 4. 高分辨率sin拟合 ==================
    kmin = k_best1 - 0.01;
    kmax = k_best1 + 0.01;
    bmin = b_best1 - 0.1;
    bmax = b_best1 + 0.1;

    for (int i = 0; i < N; ++i) {
        k_vals[i] = kmin + (kmax - kmin) * i / (N - 1);
        b_vals[i] = bmin + (bmax - bmin) * i / (N - 1);
    }

    double bestR2 = 0.0;
    double best_k = 0.0, best_b = 0.0;
    double best_a = 0.0, best_phase = 0.0;
    std::vector<double> Q_corrected(Q.size());   // 修正后的 Q 值

    for (int i = 0; i < N; ++i) {
        double k = k_vals[i];
        for (int j = 0; j < N; ++j) {
            double b = b_vals[j];

            std::vector<double> y_cor(Q.size()), x_fit(Q.size());
            for (size_t m = 0; m < Q.size(); ++m) {
                double y_fit = k * V[m] + b;
                y_cor[m] = Q[m] - y_fit;
                x_fit[m] = Z[m] - y_fit;
            }

            TGraph* grSin = new TGraph(Q.size(), &x_fit[0], &y_cor[0]);
            double xmin = *std::min_element(x_fit.begin(), x_fit.end());
            double xmax = *std::max_element(x_fit.begin(), x_fit.end());
            if (xmax - xmin < 1e-9) {
                delete grSin;
                continue;
            }

            TF1* fSin = new TF1("sinFitFine", "[0]*sin(pi/16*(x-[1]))", xmin, xmax);
            fSin->SetParameter(0, 0.7);
            fSin->SetParameter(1, 0.03);
            grSin->Fit(fSin, "Q");
            double r2_cur = ComputeR2(grSin, fSin);

            if (r2_cur > bestR2) {
                bestR2 = r2_cur;
                best_k = k;
                best_b = b;
                best_a = fSin->GetParameter(0);
                best_phase = fSin->GetParameter(1);
                // 计算修正后的 Q: Q_corr = Q - a * sin(π/16*(V - phase))
                for (size_t m = 0; m < Q.size(); ++m) {
                    Q_corrected[m] = Q[m] - best_a * sin(pi16 * (V[m] - best_phase));
                }
            }

            delete grSin;
            delete fSin;
        }
    }

    // ================== 5. 输出最终结果 ==================
    std::cout << "\n========== Final Result ==========" << std::endl;
    std::cout << "Linear part: Q_0 = " << best_k << " * (v/v_0 * Z^(1/3) + " << best_b << std::endl;
    std::cout << "Sinusoidal part: a * sin(pi/16 * (Z - Q_0 - phase))" << std::endl;
    std::cout << "a     = " << best_a << std::endl;
    std::cout << "phase = " << best_phase << " rad" << std::endl;
    std::cout << "R^2    = " << bestR2 << std::endl;

      }