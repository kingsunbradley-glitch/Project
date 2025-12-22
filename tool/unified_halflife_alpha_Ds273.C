#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <TMath.h>

// 定义物理常数
const double r0 = 1.2249;      // fm
const double e2 = 1.43996;     // MeV*fm
const double d_centrifugal = 0.0669;


// 新增辅助函数：将半衰期（秒）格式化为易读的字符串
std::string format_half_life(double T_sec) {
    char buffer[50];
    if (T_sec >= 1.0) {
        snprintf(buffer, sizeof(buffer), "%.2f s", T_sec);
    } else if (T_sec >= 1e-3) {
        snprintf(buffer, sizeof(buffer), "%.2f ms", T_sec * 1e3);
    } else if (T_sec >= 1e-6) {
        snprintf(buffer, sizeof(buffer), "%.2f us", T_sec * 1e6);
    } else {
        snprintf(buffer, sizeof(buffer), "%.2f ns", T_sec * 1e9);
    }
    return std::string(buffer);
}

// 函数：计算角动量 l
int calculate_l(double jp, int pp, double jd, int pd) {
    double delta_j_float = std::abs(jp - jd);
    int delta_j = static_cast<int>(round(delta_j_float)); // Round to nearest integer for comparison
    
    bool parity_conserved = (pp == pd);

    if (delta_j % 2 == 0) { // even delta_j
        return parity_conserved ? delta_j : delta_j + 1;
    } else { // odd delta_j
        return parity_conserved ? delta_j + 1 : delta_j;
    }
}

// 函数：计算半衰期 log10(T_1/2)
double calculate_logT12(int Z, int Ap, double Qa, double jp, int pp, double jd, int pd) {
    // 1. 基本参数
    int Ad = Ap - 4;
    int Zd = Z - 2;
    int N = Ap - Z;
    int A_alpha = 4;

    // 2. 计算 F(Z)
    double F_Z = 28.274 * sqrt(Z) + 2920.347 / Z - 204.086;

    // 3. 计算 X
    double X = r0 * (pow(Ad, 1.0/3.0) + pow(A_alpha, 1.0/3.0)) * Qa / (2.0 * Zd * e2);

    // 4. 计算 arccos 项
    double X_term = acos(sqrt(X)) - sqrt(X * (1.0 - X));

    // 5. 计算 C(Z,N)
    double C_ZN = 0.0;
    if (Z >= 78 && Z <= 82 && N > 100 && N < 126) {
        C_ZN = 1.547 - 0.077 * (82 - Z) - 0.050 * (126 - N);
    } else if (Z > 82 && Z <= 90 && N >= 110 && N <= 126) {
        C_ZN = 1.397 - 0.116 * (Z - 82) - 0.061 * (126 - N);
    }

    // 6. 确定 h
    double h_blocking = 0.0;
    if ((Z % 2 != 0) && (N % 2 != 0)) { // odd-odd
        h_blocking = 0.4036;
    } else if ((Z % 2 == 0) && (N % 2 == 0)) { // even-even
        h_blocking = 0.0;
    } else { // odd-A
        h_blocking = 0.2018;
    }

    // 7. 计算 l
    int l = calculate_l(jp, pp, jd, pd);

    // 8. 组合所有项
    double term1 = F_Z * sqrt(static_cast<double>(Ad) / (Ap * Qa)) * X_term;
    double term_l = d_centrifugal * l * (l + 1);
    
    double logT12 = term1 - 20.446 + C_ZN + term_l + h_blocking;
    
    return logT12;
}

// 主执行函数
void calculate_Ds273() {
    std::cout << "======================================================================" << std::endl;
    std::cout << "Exhaustive Calculation for 273Ds Isomers using Xu et al. formula" << std::endl;
    std::cout << "======================================================================" << std::endl;
    
    // --- 定义 273Ds^a (长寿命) 的参数 ---
    int Z_a = 110, Ap_a = 273;
    double Qa_a = 11.09; // MeV
    double jp_a = 11.0/2.0, jd_a = 9.0/2.0;
    double T_exp_a = -1.52;

    std::cout << "\n--- Case 1: 273Ds^a (long-lived isomer) ---" << std::endl;
    std::cout << "Q_alpha = " << Qa_a << " MeV, Transition: " << jp_a << " -> " << jd_a << std::endl;
    std::cout << "Experimental log10(T_1/2) approx: " << T_exp_a << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "jp_pi -> jd_pi | Parity Conserved? | l | log10(T_1/2) [s]" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    for (int pp : {1, -1}) {
        for (int pd : {1, -1}) {
             // 在这里添加新代码
            // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            double logT = calculate_logT12(Z_a, Ap_a, Qa_a, jp_a, pp, jd_a, pd);
            double T_half_life = pow(10.0, logT); // 计算半衰期
            std::string T_str = format_half_life(T_half_life); // 格式化半衰期

            int l_val = calculate_l(jp_a, pp, jd_a, pd);
            std::string pp_str = (pp > 0) ? "+" : "-";
            std::string pd_str = (pd > 0) ? "+" : "-"; 
            printf("  %s -> %s     | %-18s | %d | %16.2f | %s\n", // 修改 printf 格式
                   pp_str.c_str(), pd_str.c_str(), (pp==pd ? "Yes" : "No"), l_val, logT, T_str.c_str());
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        }
    }
    
    // --- 定义 273Ds^b (短寿命) 的参数 ---
    int Z_b = 110, Ap_b = 273;
    double Qa_b = 11.27; // MeV
    double jp_b = 1.0/2.0, jd_b = 1.0/2.0;
    double T_exp_b = -3.74;

    std::cout << "\n--- Case 2: 273Ds^b (short-lived isomer/ground state) ---" << std::endl;
    std::cout << "Q_alpha = " << Qa_b << " MeV, Transition: " << jp_b << " -> " << jd_b << std::endl;
    std::cout << "Experimental log10(T_1/2) approx: " << T_exp_b << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "jp_pi -> jd_pi | Parity Conserved? | l | log10(T_1/2) [s]" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    
    for (int pp : {1, -1}) {
        for (int pd : {1, -1}) {
            // 在这里添加新代码
            // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
            double logT = calculate_logT12(Z_b, Ap_b, Qa_b, jp_b, pp, jd_b, pd);
            double T_half_life = pow(10.0, logT); // 计算半衰期
            std::string T_str = format_half_life(T_half_life); // 格式化半衰期

            int l_val = calculate_l(jp_b, pp, jd_b, pd);
            std::string pp_str = (pp > 0) ? "+" : "-";
            std::string pd_str = (pd > 0) ? "+" : "-";
            printf("  %s -> %s     | %-18s | %d | %16.2f | %s\n", // 修改 printf 格式
                   pp_str.c_str(), pd_str.c_str(), (pp==pd ? "Yes" : "No"), l_val, logT, T_str.c_str());
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        }
    }
    std::cout << "======================================================================" << std::endl;

}