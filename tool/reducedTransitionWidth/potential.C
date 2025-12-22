#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TAxis.h"
#include <iostream>

using namespace std;

// ================================================================
//  第一部分：辅助计算函数 (必须放在主函数上面)
// ================================================================

// 1. 核势 (Igo Potential)
Double_t Vn(Double_t *x, Double_t *Massnum)     
{
    // -1100 * exp(...)
    return -1100.*TMath::Exp(-((x[0]-1.17*TMath::Power(Massnum[0],1./3.))/0.574));
}

// 2. 库仑势 (改名为 Vcoul 以避免与 ROOT Vc 库冲突)
Double_t Vcoul(Double_t *x, Double_t *Chargenum)   
{
    return 2.*Chargenum[0]*1.4399644/x[0];
}

// 3. 离心势
Double_t Vcent(Double_t *x, Double_t *par)      
{
    Double_t Rmass = 4.0*par[0]/(4.+par[0])*931.494013;
    return par[1]*(par[1]+1.)*197.32696*197.32696/(2.*Rmass*x[0]*x[0]);
}

// 4. 总势能
Double_t Vtot(Double_t *x, Double_t *par)       
{
    // par[0]=A_daughter, par[1]=Z_daughter, par[2]=A_daughter (for Rmass), par[3]=L
    return Vn(x,&par[0]) + Vcoul(x,&par[1]) + Vcent(x,&par[3]); 
}

// 5. WKB 积分函数
Double_t Func(Double_t *x, Double_t *par)
{
    // par[4] 是总能量 Etot
    Double_t val = Vtot(x,par) - par[4];
    if (val < 0) return 0; 
    return TMath::Sqrt(val);
}

// ================================================================
//  第二部分：主函数 (改为带参数的形式)
// ================================================================

/**
 * 输入参数说明:
 * Am_in     : 母核质量数 (例如 273)
 * Zm_in     : 母核质子数 (例如 110)
 * Ealpha_in : 实验室测量的 Alpha 粒子能量 (MeV) (例如 11.09)
 * T12_in    : 半衰期 (秒) (例如 0.00017)
 * BR        : 分支比 (默认为 1.0，即 100%)
 */
void potential(Double_t Am_in, Double_t Zm_in, Double_t Ealpha_in, Double_t T12_in, Double_t BR = 1.0)
{
    // --- 1. 参数初始化 ---
    Double_t Am = Am_in;
    Double_t Zm = Zm_in;
    Double_t Ealpha = Ealpha_in;
    Double_t T_half_life = T12_in;
    Double_t Intensity = BR;
    
    Double_t L = 0.0; // 默认角动量 L=0 (基态到基态)
    Double_t A = Am - 4.0; // 子核质量数
    Double_t Z = Zm - 2.0; // 子核质子数

    // 误差设为0 (因为函数输入没包含误差项，简化显示)
    Double_t T_lefterror = 0.0;
    Double_t T_righterror = 0.0;
    Double_t IRerror = 0.0;
    Double_t ILerror = 0.0;

    // 物理常数
    Double_t Hplk = 4.13566727e-21;   // MeV*s (Planck Constant h)
    Double_t Rmass = 4.*A/(4.+A)*931.494013; // Reduced Mass
    Double_t eSqua = 1.4399644;       // e^2
    Double_t Hbarc = 197.32696;       // hbar*c
    
    // --- 2. 计算有效 Q 值 (Etot) ---
    // 电子屏蔽修正
    Double_t Escr = (65.3*TMath::Power(Zm,7./5.)-80.*TMath::Power(Zm,2./5.))*1.0E-6;
    
    // 总能量 = Lab能量 + 反冲能 + 屏蔽能
    // 原代码逻辑: Etot = Ealpha + 4./A*Ealpha + Escr
    Double_t Etot = Ealpha * (1.0 + 4.0/A) + Escr;

    cout << "==================================================" << endl;
    cout << " Input: Parent Z=" << Zm << " A=" << Am << endl;
    cout << " E_alpha(Lab) = " << Ealpha << " MeV" << endl;
    cout << " T1/2         = " << T_half_life << " s" << endl;
    cout << "--------------------------------------------------" << endl;
    cout << " [Calc] Q_eff (Etot) = " << Etot << " MeV" << endl;

    // --- 3. 寻找外转折点 Rout ---
    Double_t Rout = Z*eSqua/Etot + TMath::Sqrt(Z*Z*eSqua*eSqua + Etot*Hbarc*Hbarc*L*(L+1)/2./Rmass)/Etot;
    cout << " [Calc] Rout         = " << Rout << " fm" << endl;
    
    // --- 4. 寻找内转折点 Rin (数值迭代) ---
    Double_t array[5] = {A, Z, A, L, Etot}; // 传递给 Vtot 的参数
    Double_t rx, r, r0, r1;
    Double_t error = 1.0E-10;
    
    // 粗略搜索
    for(rx = Rout-0.2; rx > 0; rx = rx-0.2) {
        r = rx;
        if(Vtot(&r, array) < Etot) break;
    }
    // 精细搜索 (牛顿迭代)
    r0 = r; r1 = r+0.2;
    int iter = 0;
    while(1) {
        Double_t v1 = Vtot(&r1, array);
        Double_t v0 = Vtot(&r0, array);
        if (v1 == v0) break;
        r = r1 - (v1 - Etot)*(r1 - r0)/(v1 - v0);
        if (TMath::Abs(r - r1) < error) break;
        r0 = r1; r1 = r; iter++;
        if (iter > 1000) break; 
    }
    Double_t Rin = r;
    cout << " [Calc] Rin          = " << Rin << " fm" << endl;

    // --- 5. WKB 积分计算穿透率 P ---
    TF1 *fFunc = new TF1("Func", Func, Rin, Rout, 5);
    fFunc->SetParameters(A, Z, A, L, Etot);
    
    Double_t Integration = fFunc->Integral(Rin, Rout);
    Double_t Penetration = TMath::Exp(-2*TMath::Sqrt(2.*Rmass)*Integration/Hbarc);
    
    cout << " [Calc] Penetration P= " << Penetration << endl;

    // --- 6. 计算约化宽度 Delta^2 ---
    // Delta^2 = h * lambda / P  (注意：这里用的是 h 不是 hbar，对应原代码逻辑)
    // Lambda = ln2 / T1/2
    Double_t DeltSqua = Hplk * (TMath::Log(2)/T_half_life * Intensity) / Penetration * 1000.0; // 转换为 keV

    cout << "--------------------------------------------------" << endl;
    cout << " Reduced Width delta^2 = " << DeltSqua << " keV" << endl;
    cout << "==================================================" << endl;

    // --- 7. 画图 (可选) ---
    TCanvas *c1 = new TCanvas("c1", "Potentials", 800, 600);
    c1->Divide(2,2);
    
    // 总势能图
    c1->cd(1);
    TF1 *V = new TF1("Total_Potential", Vtot, 5, 30, 4);
    V->SetParameters(A, Z, A, L);
    V->SetTitle("Total Potential V(r);r [fm];V [MeV]");
    V->SetLineColor(kRed);
	V->GetYaxis()->SetRangeUser(-20.0, 50.0); // 强制 Y 轴只显示 -20 到 50 MeV
    V->GetXaxis()->SetRangeUser(6.0, 20.0);   // X 轴也缩短一点，聚焦势垒
    V->Draw();
    // 画一条能量线
    TLine *line = new TLine(5, Etot, 30, Etot);
    line->SetLineStyle(2);
    line->SetLineColor(kBlue);
    line->Draw();

    // 核势
    c1->cd(2);
    TF1 *fVn = new TF1("Nuclear_Vn", Vn, 5, 20, 1);
    fVn->SetParameter(0, A);
    fVn->SetTitle("Nuclear Potential (Igo);r [fm];V [MeV]");
    fVn->Draw();

    // 库仑势
    c1->cd(3);
    TF1 *fVc = new TF1("Coulomb_Vc", Vcoul, 5, 40, 1);
    fVc->SetParameter(0, Z);
    fVc->SetTitle("Coulomb Potential;r [fm];V [MeV]");
    fVc->Draw();

    // 离心势
    c1->cd(4);
    TF1 *fVcent = new TF1("Centrifugal", Vcent, 5, 20, 2);
    fVcent->SetParameters(A, L);
    fVcent->SetTitle("Centrifugal Potential;r [fm];V [MeV]");
    fVcent->Draw();
}