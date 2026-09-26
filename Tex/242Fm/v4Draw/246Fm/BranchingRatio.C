#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "TEfficiency.h"

void BranchingRatio()
{
    int nModes;

    std::cout << "Input number of decay modes: ";
    std::cin >> nModes;

    if (nModes <= 0) {
        std::cout << "Error: number of modes must be > 0." << std::endl;
        return;
    }

    // Confidence levels
    std::vector<double> CLs = {
        0.683,
        0.84,
        0.90,
        0.954,
        0.9973
    };

    // ============================================================
    // Special case:
    // Only one decay mode observed.
    //
    // Calculate one-sided lower limit for this mode
    // and one-sided upper limit for all other possible modes.
    //
    // If N successes are observed in N trials:
    //
    //     b_lower = (1 - CL)^(1/N)
    //
    // therefore:
    //
    //     b_other_upper = 1 - b_lower
    // ============================================================
    if (nModes == 1) {

        int N;

        std::cout << "Input counts for Mode 1: ";
        std::cin >> N;

        if (N <= 0) {
            std::cout << "Error: counts must be > 0." << std::endl;
            return;
        }

        std::cout << "\nTotal counts = " << N << std::endl;
        std::cout << "One-sided confidence limits\n" << std::endl;

        std::cout << std::fixed << std::setprecision(6);

        std::cout
            << std::setw(12) << "CL"
            << std::setw(20) << "Mode1 lower"
            << std::setw(20) << "Other upper"
            << std::endl;

        for (double CL : CLs) {

            double bModeLower =
                std::pow(1.0 - CL, 1.0 / static_cast<double>(N));

            double bOtherUpper =
                1.0 - bModeLower;

            std::cout
                << std::setw(12) << CL
                << std::setw(20) << bModeLower
                << std::setw(20) << bOtherUpper
                << std::endl;

            std::cout
                << "CL = " << 100.0 * CL << "%"
                << " : Mode 1 > "
                << 100.0 * bModeLower << "%"
                << ", Other < "
                << 100.0 * bOtherUpper << "%"
                << std::endl;
        }

        return;
    }

    // ============================================================
    // General case: two or more observed decay modes
    // Exact two-sided Clopper-Pearson confidence intervals
    // ============================================================

    std::vector<int> counts(nModes);

    int Ntotal = 0;

    for (int i = 0; i < nModes; ++i) {

        std::cout
            << "Input counts for Mode "
            << i + 1
            << ": ";

        std::cin >> counts[i];

        if (counts[i] < 0) {
            std::cout
                << "Error: counts cannot be negative."
                << std::endl;
            return;
        }

        Ntotal += counts[i];
    }

    if (Ntotal <= 0) {
        std::cout
            << "Error: total counts must be > 0."
            << std::endl;
        return;
    }

    std::cout
        << "\nTotal counts = "
        << Ntotal
        << "\n"
        << std::endl;

    std::cout
        << std::fixed
        << std::setprecision(6);

    for (int i = 0; i < nModes; ++i) {

        int k = counts[i];

        double BR =
            static_cast<double>(k)
            / static_cast<double>(Ntotal);

        std::cout
            << "Decay Mode "
            << i + 1
            << "   Counts = "
            << k
            << "   BR = "
            << BR
            << std::endl;

        std::cout
            << std::setw(12) << "CL"
            << std::setw(15) << "BR"
            << std::setw(15) << "-Error"
            << std::setw(15) << "+Error"
            << std::setw(15) << "Lower"
            << std::setw(15) << "Upper"
            << std::endl;

        for (double CL : CLs) {

            double lower =
                TEfficiency::ClopperPearson(
                    Ntotal,
                    k,
                    CL,
                    false
                );

            double upper =
                TEfficiency::ClopperPearson(
                    Ntotal,
                    k,
                    CL,
                    true
                );

            double errLow  = BR - lower;
            double errHigh = upper - BR;

            std::cout
                << std::setw(12) << CL
                << std::setw(15) << BR
                << std::setw(15) << errLow
                << std::setw(15) << errHigh
                << std::setw(15) << lower
                << std::setw(15) << upper
                << std::endl;
        }

        std::cout << std::endl;
    }
}