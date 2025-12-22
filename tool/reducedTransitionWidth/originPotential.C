#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TAxis.h"
#include <iostream>

using namespace std;

//originPotential(Double_t Am = 269, Double_t Zm = 110., Double_t Ealpha = 11.11, Double_t T_half_life = 179E-6, Double_t BR = 1.0)!20.2949 keV
//originPotential(Double_t Am = 270, Double_t Zm = 110., Double_t Ealpha = 11.1, Double_t T_half_life = 0.2E-3, Double_t BR = 0.5) !9.23117 keV
//originPotential(Double_t Am = 271, Double_t Zm = 110., Double_t Ealpha = 10.73, Double_t T_half_life = 1.63E-3, Double_t BR = 0.81)!10.9373 MeV
//originPotential(Double_t Am = 273, Double_t Zm = 110., Double_t Ealpha = 10.88, Double_t T_half_life = 17.34E-3, Double_t BR = 0.22)!originPotential(Double_t Am = 273, Double_t Zm = 110., Double_t Ealpha = 11.361, Double_t T_half_life = 0.16E-3, Double_t BR = 0.78)
//0.133811 keV  \  4.31912 keV
//originPotential(Double_t Am = 275, Double_t Zm = 110., Double_t Ealpha = 11.2, Double_t T_half_life = 0.43E-3, Double_t BR = 1.0)!4.31753 keV
//originPotential(Double_t Am = 276, Double_t Zm = 110., Double_t Ealpha = 10.75, Double_t T_half_life = 0.36E-3, Double_t BR = 0.43)!22.8526 keV

//originPotential(Double_t Am = 277, Double_t Zm = 110., Double_t Ealpha = 10.57, Double_t T_half_life = 6E-3, Double_t BR = 1.0)!8.31363 keV
//originPotential(Double_t Am = 279, Double_t Zm = 110., Double_t Ealpha = 9.71, Double_t T_half_life = 0.21E-3, Double_t BR = 0.11)!4221.24 keV
//originPotential(Double_t Am = 280, Double_t Zm = 110., Double_t Ealpha = 11.09, Double_t T_half_life = 0.17E-3, Double_t BR = 1.0)1.0SF!
//originPotential(Double_t Am = 281, Double_t Zm = 110., Double_t Ealpha = 8.73, Double_t T_half_life = 12.7, Double_t BR = 0.07)!39.2493 keV


//originPotential( 273,  110., 10.857, 17.34E-3,31.55E-3,6.83E-3, 1) 273Ds g.s. 0.687819 keV    -0.270923 keV    +1.25148 keV
//originPotential( 273,  110., 11.165, 0.16E-3,0.08e-3,0.05E-3, 1) 273Ds i.s.14.8718 keV    -4.64745 keV +7.43592 keV
//originPotential( 269,  108., 8.950, 6.49, 31.09,2.97, 0.83) 269Hs i.s.  49.5836 keV   -22.6908 keV +237.528 keV
//originPotential( 269,  108., 9.130, 9.66,4.41,2.74, 1) 269Hs g.s.11.2545 keV  -3.19227 keV    +5.13792 keV
//originPotential( 269,  108., 9.13, 10.12,5.03,3.03,0.9 ) 269Hs 1???   9.66862keV      -2.89485 keV    +4.80565 keV
//originPotential( 269,  108., 8.95, 10.12,5.03,3.03,0.1 ) 269Hs 1???   3.83111 keV     -1.14741 keV    +1.90718 keV
//originPotential( 265,  106., 8.69, 14.4,3.7,2.5, 0.8) 265Sg i.s.29.8284 keV -5.17854 keV +7.66424 keV
//originPotential( 269,  108., 8.88, 11.9,56.8,5.41, 0.17) 269Hs   
//originPotential( 265,  106., 8.84, 8.5,2.6,1.6, 0.91) 265Sg 
 //originPotential( 265,  106., 8.84, 8.5,2.6,1.6, 0.09) 265Sg
//originPotential( 261,  104., 8.28, 68,3,3, 1) 261Rf 
//originPotential( 261,  104., 8.51, 2.6 , 0.7 , 0.5 , 0.18) 261Rf
//originPotential( 265,  106., 8.69, 14.4,3.7,2.5, 0.2)




// ==========================================
// 第一步：把辅助函数移到最前面 (必须先定义，后使用)
// ==========================================

//--------- Alpha-nuclear potential (Nuclear Igo Potential)
Double_t Vn(Double_t *x, Double_t *Massnum)     // MeV
{
    Double_t Vn_val = -1100.*TMath::Exp(-((x[0]-1.17*TMath::Power(Massnum[0],1./3.))/0.574));
    return Vn_val;
}

//--------- Coulomb potential
// 【修改点】：名字改为 Vcoul，避免与 ROOT Vc 库冲突
Double_t Vcoul(Double_t *x, Double_t *Chargenum)   // MeV
{
    Double_t Vc_val = 2.*Chargenum[0]*1.4399644/x[0];
    return Vc_val;
}

//--------- Centrifugal potential
Double_t Vcent(Double_t *x, Double_t *par)      // MeV
{
    Double_t Rmass = 4.0*par[0]/(4.+par[0])*931.494013;
    Double_t Vcent_val = par[1]*(par[1]+1.)*197.32696*197.32696/(2.*Rmass*x[0]*x[0]);
    return Vcent_val;
}

//--------- Total potential
Double_t Vtot(Double_t *x, Double_t *par)       // MeV
{
    // 【修改点】：这里调用 Vcoul 而不是 Vc
    Double_t Vtot_val = Vn(x,&par[0]) + Vcoul(x,&par[1]) + Vcent(x,&par[2]);
    return Vtot_val;
}

//--------- Function for integration
Double_t Func(Double_t *x, Double_t *par)
{
    Double_t val = Vtot(x,par) - par[4];
    if (val < 0) return 0; // 安全检查
    Double_t Func_val = TMath::Sqrt(val);
    return Func_val;
}

// ==========================================
// 第二步：主函数放在最后
// ==========================================

void originPotential(Double_t Am = 273, Double_t Zm = 110., Double_t Ealpha = 11.09, Double_t T_half_life = 0.17E-3, Double_t T_righterror = 0.02E-3, Double_t T_lefterror = 0.02E-3,Double_t BR = 1.0) //right+ left-
{
    //--------- Input parameters
    Double_t A, Z, Escr, Etot;
    
    // 1. 修改母核信息
    //Double_t Am = 273, Zm = 110.;
    Double_t L = 0.; 
    
    A = Am - 4.; 
    Z = Zm - 2.;

    // 2. 修改 Alpha 粒子能量 (实验室坐标系, MeV)
    //Double_t Ealpha = 11.09;           
    
    // 3. 修改半衰期 (单位: 秒)
    //Double_t T_half_life = 0.17E-3;           

    // 4. 修改误差
    //Double_t T_lefterror = 0.02E-3;       
    //Double_t T_righterror = 0.02E-3;       

    // 5. 分支比
    Double_t Intensity = BR;     
    Double_t IRerror = 0.0;       
    Double_t ILerror = 0.0;      

    //---------------------------------------------------  

    Double_t Hplk = 4.13566727e-21;           // MeV*s
    Double_t Rmass = 4.*A/(4.+A)*931.494013;    // MeV/c2
    Double_t eSqua = 1.4399644;               // MeV*fm
    Double_t Hbarc = 197.32696;             // MeV*fm
    Double_t Rout, Rin;                    // fm
    Double_t array[4] = {A, Z, A, L};
    Double_t rx, r0, r1, r;
    Double_t error = 1.0E-10;
    Int_t i = 0;

    //--------- Total decay energy (计算放在画图前比较好，确保参数正确)
    Escr = (65.3*TMath::Power(Zm,7./5.)-80.*TMath::Power(Zm,2./5.))*1.0E-6;   
    
    // 计算有效 Q 值 (考虑反冲和屏蔽)
    // Etot = E_alpha * (Mass_Parent / Mass_Daughter) + E_screening
    Etot = Ealpha * (Am/A) + Escr;  
    // 原代码里的写法: Etot = Ealpha + 4./A*Ealpha + Escr; 也是对的，数学上等价

    cout << "----------------------------------" << endl;
    cout << "Etot (Q_eff) = " << Etot << " MeV" << endl;

    //--------- Potential plot
    TCanvas *c1 = new TCanvas("c1", "c1", 600, 0, 600, 500);
    c1->Divide(2,2);

    c1->cd(1);
        // 现在 Vtot 已经在上面定义过了，这里不会报错了
        TF1 *V = new TF1("Vn+Vc+Vcent", Vtot, 7, 50, 4);
        V->SetParameters(A, Z, A, L);
        V->SetTitle("Total Potential;r (fm);V (MeV)");
        V->GetYaxis()->SetRangeUser(-20.0, 50.0); // 强制 Y 轴只显示 -20 到 50 MeV
        V->GetXaxis()->SetRangeUser(6.0, 20.0);   // X 轴也缩短一点，聚焦势垒
        V->Draw();

    c1->cd(2);
        TF1 *fVn = new TF1("Vn", Vn, 7, 30, 1);
        fVn->SetParameter(0, A);
        fVn->SetTitle("Nuclear Potential;r (fm);V (MeV)");
        fVn->Draw();

    c1->cd(3);
        // 【修改点】：这里名字也改一下，避免歧义
        TF1 *fVc = new TF1("Vcoul", Vcoul, 7, 30, 1);
        fVc->SetParameter(0, Z);
        fVc->SetTitle("Coulomb Potential;r (fm);V (MeV)");
        fVc->Draw();

    c1->cd(4);
        TF1 *fVcent = new TF1("Vcent", Vcent, 7, 30, 2);
        fVcent->SetParameters(A, L);
        fVcent->SetTitle("Centrifugal Potential;r (fm);V (MeV)");
        fVcent->Draw();


    //--------- Outer turning point
    Rout = Z*eSqua/Etot + TMath::Sqrt(Z*Z*eSqua*eSqua + Etot*Hbarc*Hbarc*L*(L+1)/2./Rmass)/Etot;
    
   
    //--------- Inner turning point
    for(rx = Rout-0.2; rx > 0; rx = rx-0.2)
    {
        r = rx;
        if(Vtot(&r, array) < Etot) break;
    }
    r0 = r; r1 = r+0.2;
    while(1)             
    {
        // 牛顿迭代找交点
        Double_t v1 = Vtot(&r1, array);
        Double_t v0 = Vtot(&r0, array);
        if (v1 == v0) break; 
        r = r1 - (v1 - Etot)*(r1 - r0)/(v1 - v0);
        if( TMath::Abs(r - r1) < error ) break;
        r0 = r1;  r1 = r;  i++;
        if(i > 1000) break; // 防止死循环
    }
    Rin = r;
    cout << "Rin = " << Rin << " fm" << " (iterative number " << i << ")" << endl;
    cout << "Rout = " << Rout << " fm" << endl;
    //--------- Integration
    TF1 *fFunc = new TF1("Func", Func, Rin, Rout, 5);
    fFunc->SetParameters(A, Z, A, L, Etot);
    Double_t Integration = fFunc->Integral(Rin, Rout);
   
    //--------- Penetration factor
    Double_t Penetration = TMath::Exp(-2*TMath::Sqrt(2.*Rmass)*Integration/Hbarc);
    cout << "Penetration factor: P = " << Penetration << endl;
    
    //--------- Reduced alpha decay width
    // 公式：delta^2 = h * lambda / P
    // T_half = ln2 / lambda
    // Delta^2 (keV)
    Double_t DeltSqua = Hplk * (TMath::Log(2)/T_half_life * Intensity) / Penetration * 1000.;

    // 误差计算
    Double_t Left_error = 1000.*Hplk*TMath::Log(2)/Penetration * TMath::Power(TMath::Power(ILerror/T_half_life,2) + TMath::Power(Intensity*T_lefterror/T_half_life/T_half_life,2), 0.5);  
    
    Double_t Right_error = 1000.*Hplk*TMath::Log(2)/Penetration * TMath::Power(TMath::Power(IRerror/T_half_life,2) + TMath::Power(Intensity*T_righterror/T_half_life/T_half_life,2), 0.5);

    cout << "----------------------------------" << endl;
    cout << "Reduced alpha decay width: Delta^2 = " << DeltSqua << " keV" << endl;
    cout << "Left error = " << "-" <<Left_error << " keV" << endl;
    cout << "Right error =" << "+" << Right_error << " keV" << endl;
    cout << "----------------------------------" << endl;

    return ;
}