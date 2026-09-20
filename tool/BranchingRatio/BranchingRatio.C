// BranchingRatio.C

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

#include "TEfficiency.h"

void BranchingRatio()
{
    int nModes;

    std::cout << "Number of observed decay modes: ";
    std::cin >> nModes;

    if (nModes < 1) {
        std::cout << "Error: number of decay modes must be >= 1." << std::endl;
        return;
    }

    const std::vector<double> CLs = {
        0.683,
        0.84,
        0.90,
        0.954,
        0.9973
    };

    // ============================================================
    // Special case:
    // Only ONE decay mode is observed.
    //
    // All N events belong to this decay mode:
    //
    // N_mode  = N
    // N_other = 0
    //
    // Use ONE-SIDED binomial confidence limits.
    // ============================================================

    if (nModes == 1) {

        int N;

        std::cout << "Counts of decay mode 1: ";
        std::cin >> N;

        if (N <= 0) {
            std::cout << "Error: counts must be > 0." << std::endl;
            return;
        }

        std::cout << "\n";
        std::cout << "===============================================================\n";
        std::cout << "Only one decay mode was observed.\n";
        std::cout << "Observed mode counts = " << N << "\n";
        std::cout << "Other modes counts   = 0\n";
        std::cout << "===============================================================\n\n";

        std::cout << std::fixed << std::setprecision(6);

        std::cout
            << std::setw(10) << "CL"
            << std::setw(20) << "Mode1 lower"
            << std::setw(20) << "Other upper"
            << "\n";

        std::cout
            << "--------------------------------------------------\n";

        for (double CL : CLs) {

            // Probability of observing zero "other" decays:
            //
            // (1 - b_other)^N = 1 - CL
            //
            // Therefore:
            //
            // b_other_upper = 1 - (1 - CL)^(1/N)
            //
            // and
            //
            // b_mode1_lower = (1 - CL)^(1/N)

            double bModeLower =
                std::pow(1.0 - CL, 1.0 / N);

            double bOtherUpper =
                1.0 - bModeLower;

            std::cout
                << std::setw(10) << CL
                << std::setw(20) << bModeLower
                << std::setw(20) << bOtherUpper
                << "\n";
        }

        std::cout << "\n";

        // Percentage form
        std::cout << "================ Percentage form ================\n\n";

        for (double CL : CLs) {

            double bModeLower =
                std::pow(1.0 - CL, 1.0 / N);

            double bOtherUpper =
                1.0 - bModeLower;

            std::cout
                << "CL = " << CL * 100.0 << "% :   "
                << "Mode 1 > " << bModeLower * 100.0 << "%,   "
                << "Other < " << bOtherUpper * 100.0 << "%"
                << "\n";
        }

        return;
    }


    // ============================================================
    // General case:
    // Two or more observed decay modes.
    //
    // Calculate the branching ratio of each mode and
    // two-sided Clopper-Pearson confidence intervals.
    // ============================================================

    std::vector<int> counts(nModes);

    int Ntotal = 0;

    std::cout << "\nInput counts for each decay mode:\n";

    for (int i = 0; i < nModes; ++i) {

        std::cout << "Mode " << i + 1 << ": ";
        std::cin >> counts[i];

        if (counts[i] < 0) {
            std::cout << "Error: counts cannot be negative." << std::endl;
            return;
        }

        Ntotal += counts[i];
    }

    if (Ntotal <= 0) {
        std::cout << "Error: total counts must be > 0." << std::endl;
        return;
    }

    std::cout << "\n";
    std::cout << "===============================================================\n";
    std::cout << "Total counts = " << Ntotal << "\n";
    std::cout << "===============================================================\n";

    std::cout << std::fixed << std::setprecision(6);

    for (int i = 0; i < nModes; ++i) {

        int k = counts[i];

        double BR =
            static_cast<double>(k) / static_cast<double>(Ntotal);

        std::cout << "\n";
        std::cout << "===============================================================\n";
        std::cout << "Decay Mode " << i + 1
                  << "   Counts = " << k
                  << "   BR = " << BR
                  << "\n";
        std::cout << "===============================================================\n";

        std::cout
            << std::setw(10) << "CL"
            << std::setw(14) << "BR"
            << std::setw(14) << "-Error"
            << std::setw(14) << "+Error"
            << std::setw(14) << "Lower"
            << std::setw(14) << "Upper"
            << "\n";

        std::cout
            << "----------------------------------------------------------------------------------\n";

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
                << std::setw(10) << CL
                << std::setw(14) << BR
                << std::setw(14) << errLow
                << std::setw(14) << errHigh
                << std::setw(14) << lower
                << std::setw(14) << upper
                << "\n";
        }
    }
}