// HF_l.C
// ------------------------------------------------------------
// Estimate the centrifugal-barrier equivalent hindrance factor
//
//      HF_l = T_1/2(l) / T_1/2(l=0)
//           ~ P_0 / P_l
//           = exp[ 2 (G_l - G_0) ]
//
// using a simple Coulomb + centrifugal WKB barrier:
//
//      V(r) = 2 Z_d e^2 / r
//           + (hbar*c)^2 l(l+1) / (2 mu c^2 r^2)
//
// This is intended as a transparent estimate / cross-check.
// It is NOT a full Rasmussen exponential-nuclear-potential
// calculation.
//
// ROOT usage:
//   root -l
//   .L HF_l.C+
//   HF_l();
//   HF_l(265, 106, 8.69, 4);
//   HF_l(265, 106, 8.84, 4);
//
// If the third argument is already Q_alpha rather than E_alpha:
//   HF_l(265, 106, 8.82, 4, 1.20, true);
//
// Author: generated for ROOT/C++ analysis
// ------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

namespace HFLCalc {

// Constants in MeV, fm
constexpr double hbarc = 197.3269804;      // MeV fm
constexpr double e2    = 1.43996448;       // MeV fm
constexpr double amu   = 931.49410242;     // MeV/c^2
constexpr int    NSTEP = 200000;            // Simpson steps, must be even

double ReducedMassMeV(int Ad)
{
    // alpha + daughter reduced mass, using mass numbers
    // mu c^2 ~= [4*Ad/(4+Ad)] u c^2
    return (4.0 * Ad / (Ad + 4.0)) * amu;
}

double ContactRadius(int Ad, double r0)
{
    return r0 * (std::cbrt(static_cast<double>(Ad)) + std::cbrt(4.0));
}

double VCoulomb(double r, int Zd)
{
    return 2.0 * Zd * e2 / r;
}

double VCentrifugal(double r, int l, double mu)
{
    if (l == 0) return 0.0;
    return hbarc * hbarc * l * (l + 1.0) / (2.0 * mu * r * r);
}

double Veff(double r, int Zd, int l, double mu)
{
    return VCoulomb(r, Zd) + VCentrifugal(r, l, mu);
}

double OuterTurningPoint(double Q, int Zd, int l, double mu)
{
    // Solve:
    // Q = k/r + c/r^2
    // => Q r^2 - k r - c = 0
    //
    // Positive root:
    // r2 = [k + sqrt(k^2 + 4 Q c)] / (2 Q)

    const double k = 2.0 * Zd * e2;
    const double c = hbarc * hbarc * l * (l + 1.0) / (2.0 * mu);

    return (k + std::sqrt(k*k + 4.0*Q*c)) / (2.0*Q);
}

double Integrand(double r, double Q, int Zd, int l, double mu)
{
    const double dV = Veff(r, Zd, l, mu) - Q;
    if (dV <= 0.0) return 0.0;

    // Dimensionless WKB wave number:
    // kappa = sqrt[2 mu (V-Q)] / (hbar c)
    return std::sqrt(2.0 * mu * dV) / hbarc;
}

double WKBAction(double Q, int Ad, int Zd, int l, double r0)
{
    const double mu = ReducedMassMeV(Ad);
    const double r1 = ContactRadius(Ad, r0);
    const double r2 = OuterTurningPoint(Q, Zd, l, mu);

    if (r2 <= r1) {
        std::cerr << "ERROR: outer turning point <= inner radius.\n";
        std::cerr << "       r1 = " << r1 << " fm, r2 = " << r2 << " fm\n";
        return NAN;
    }

    int n = NSTEP;
    if (n % 2 != 0) ++n;

    const double h = (r2 - r1) / n;

    // Simpson integration
    double sum = Integrand(r1, Q, Zd, l, mu)
               + Integrand(r2, Q, Zd, l, mu);

    for (int i = 1; i < n; ++i) {
        const double r = r1 + i*h;
        sum += (i % 2 ? 4.0 : 2.0) * Integrand(r, Q, Zd, l, mu);
    }

    return sum * h / 3.0;
}

} // namespace HFLCalc


void HF_l(int Ap = 265,
          int Zp = 106,
          double Ealpha_or_Q = 8.69,
          int l = 4,
          double r0 = 1.20,
          bool inputIsQ = false)
{
    using namespace HFLCalc;

    if (Ap <= 4 || Zp <= 2) {
        std::cerr << "ERROR: invalid parent A or Z.\n";
        return;
    }
    if (l < 0) {
        std::cerr << "ERROR: l must be >= 0.\n";
        return;
    }
    if (Ealpha_or_Q <= 0.0) {
        std::cerr << "ERROR: E_alpha or Q_alpha must be positive.\n";
        return;
    }

    const int Ad = Ap - 4;
    const int Zd = Zp - 2;

    // Recoil correction:
    // Q_alpha ~= E_alpha * Ap/Ad
    const double Q = inputIsQ
                   ? Ealpha_or_Q
                   : Ealpha_or_Q * static_cast<double>(Ap) / static_cast<double>(Ad);

    const double mu = ReducedMassMeV(Ad);
    const double R  = ContactRadius(Ad, r0);

    const double r20 = OuterTurningPoint(Q, Zd, 0, mu);
    const double r2l = OuterTurningPoint(Q, Zd, l, mu);

    const double G0 = WKBAction(Q, Ad, Zd, 0, r0);
    const double Gl = WKBAction(Q, Ad, Zd, l, r0);

    if (!std::isfinite(G0) || !std::isfinite(Gl)) {
        std::cerr << "ERROR: WKB integration failed.\n";
        return;
    }

    const double deltaG = Gl - G0;

    // HF_l = exp[2(Gl-G0)]
    const double log10HF = 2.0 * deltaG / std::log(10.0);
    const double HFl = std::pow(10.0, log10HF);

    // P_l/P_0 = 1/HF_l
    const double PlOverP0 = 1.0 / HFl;

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n============================================================\n";
    std::cout << " Centrifugal equivalent hindrance factor HF_l\n";
    std::cout << "============================================================\n";
    std::cout << "Parent       : A = " << Ap << ", Z = " << Zp << "\n";
    std::cout << "Daughter     : A = " << Ad << ", Z = " << Zd << "\n";
    std::cout << "l            : " << l << "\n";
    std::cout << "r0           : " << r0 << " fm\n";

    if (inputIsQ) {
        std::cout << "Input Qalpha : " << Ealpha_or_Q << " MeV\n";
    } else {
        std::cout << "Input Ealpha : " << Ealpha_or_Q << " MeV\n";
        std::cout << "Qalpha       : " << Q
                  << " MeV  [Ealpha * Ap/Ad recoil correction]\n";
    }

    std::cout << "mu c^2       : " << mu << " MeV\n";
    std::cout << "R_inner      : " << R << " fm\n";
    std::cout << "r2(l=0)      : " << r20 << " fm\n";
    std::cout << "r2(l=" << l << ")      : " << r2l << " fm\n";

    std::cout << "\nWKB action:\n";
    std::cout << "G_0          : " << G0 << "\n";
    std::cout << "G_l          : " << Gl << "\n";
    std::cout << "Delta G      : " << deltaG << "\n";

    std::cout << "\nResult:\n";
    std::cout << "P_l / P_0    : " << PlOverP0 << "\n";
    std::cout << "log10(HF_l)  : " << log10HF << "\n";
    std::cout << "HF_l         : " << HFl << "\n";
    std::cout << "============================================================\n\n";

    std::cout << "Definition used:\n";
    std::cout << "HF_l = P_0/P_l = exp[2(G_l-G_0)]\n\n";
}
