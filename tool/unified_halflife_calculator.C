#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits> // Required for numeric_limits

// Enum to make parity representation clearer
enum Parity { UNKNOWN = 0, PLUS = 1, MINUS = -1 };

// --- Physical Constants and Formula Parameters ---
const double r0 = 1.2249;      // fm
const double e2 = 1.43996;     // MeV*fm
const double d_centrifugal = 0.0669;

// --- Helper Functions ---

// Converts a Parity enum to its string representation (+ or -)
std::string parity_to_string(Parity p) {
    if (p == PLUS) return "+";
    if (p == MINUS) return "-";
    return "?";
}

// Formats the half-life (in seconds) into a human-readable string (s, ms, us, ns)
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


// --- Core Calculation Functions ---

// Calculates the orbital angular momentum 'l' based on spin-parity selection rules
int calculate_l(double jp, Parity pp, double jd, Parity pd) {
    if (pp == UNKNOWN || pd == UNKNOWN) return -1; // Cannot calculate without parity
    double delta_j_float = std::abs(jp - jd);
    int delta_j = static_cast<int>(round(delta_j_float));
    
    bool parity_conserved = (pp == pd);

    if (delta_j % 2 == 0) { // even delta_j
        return parity_conserved ? delta_j : delta_j + 1;
    } else { // odd delta_j
        return parity_conserved ? delta_j + 1 : delta_j;
    }
}

// The main function implementing the Xu et al. formula
double calculate_logT12(int Z, int Ap, double Qa, double jp, Parity pp, double jd, Parity pd, int& out_l) {
    // 1. Basic parameters
    int Ad = Ap - 4;
    int Zd = Z - 2;
    int N = Ap - Z;

    // 2. Calculate l (and check for valid parities)
    out_l = calculate_l(jp, pp, jd, pd);
    if (out_l == -1) return std::numeric_limits<double>::quiet_NaN();

    // 3. F(Z) term
    double F_Z = 28.274 * sqrt(Z) + 2920.347 / Z - 204.086;

    // 4. X term and arccos term
    double X = r0 * (pow(Ad, 1.0/3.0) + pow(4.0, 1.0/3.0)) * Qa / (2.0 * Zd * e2);
    if (X > 1.0 || X < 0.0) {
        std::cerr << "Error: X is outside the valid range [0, 1]. X = " << X << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    double X_term = acos(sqrt(X)) - sqrt(X * (1.0 - X));

    // 5. Shell effect C(Z,N)
    double C_ZN = 0.0;
    if (Z >= 78 && Z <= 82 && N > 100 && N < 126) {
        C_ZN = 1.547 - 0.077 * (82 - Z) - 0.050 * (126 - N);
    } else if (Z > 82 && Z <= 90 && N >= 110 && N <= 126) {
        C_ZN = 1.397 - 0.116 * (Z - 82) - 0.061 * (126 - N);
    }

    // 6. Hindrance term h
    double h_blocking = 0.0;
    if ((Z % 2 != 0) && (N % 2 != 0)) h_blocking = 0.4036;      // odd-odd
    else if ((Z % 2 != 0) || (N % 2 != 0)) h_blocking = 0.2018; // odd-A

    // 7. Centrifugal term
    double term_l = d_centrifugal * out_l * (out_l + 1);
    
    // 8. Combine all terms
    double term1 = F_Z * sqrt(static_cast<double>(Ad) / (Ap * Qa)) * X_term;
    return term1 - 20.446 + C_ZN + term_l + h_blocking;
}

// --- User Interface and Execution Functions ---

// A single calculation and printout, used by other functions
void perform_and_print_calculation(int Z, int Ap, double Qa, double jp, Parity pp, double jd, Parity pd) {
    int l_val;
    double logT = calculate_logT12(Z, Ap, Qa, jp, pp, jd, pd, l_val);

    if (std::isnan(logT)) {
        printf("  %s -> %s     | %-18s | - | Calculation failed.\n",
               parity_to_string(pp).c_str(), parity_to_string(pd).c_str(), "N/A");
        return;
    }

    double T_half_life = pow(10.0, logT);
    std::string T_str = format_half_life(T_half_life);
    
    printf("  %s -> %s     | %-18s | %d | %16.2f | %s\n",
           parity_to_string(pp).c_str(), parity_to_string(pd).c_str(),
           (pp == pd ? "Yes" : "No"), l_val, logT, T_str.c_str());
}

// Example function for the specific 273Ds case, for demonstration
void calculate_Ds273_example() {
    std::cout << "\n=======================================================================" << std::endl;
    std::cout << "Running pre-defined example for 273Ds Isomers..." << std::endl;
    std::cout << "=======================================================================" << std::endl;
    
    // Case 1: 273Ds^a (long-lived)
    std::cout << "\n--- Case 1: 273Ds^a (long-lived isomer) ---" << std::endl;
    std::cout << "Q_alpha = 11.09 MeV, Transition: 5.5 -> 4.5" << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;
    std::cout << "jp_pi -> jd_pi | Parity Conserved? | l | log10(T_1/2) [s] |   T_1/2 " << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;
    for (Parity pp : {PLUS, MINUS}) {
        for (Parity pd : {PLUS, MINUS}) {
            perform_and_print_calculation(110, 273, 11.09, 5.5, pp, 4.5, pd);
        }
    }
    
    // Case 2: 273Ds^b (short-lived)
    std::cout << "\n--- Case 2: 273Ds^b (short-lived isomer/ground state) ---" << std::endl;
    std::cout << "Q_alpha = 11.27 MeV, Transition: 0.5 -> 0.5" << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;
    std::cout << "jp_pi -> jd_pi | Parity Conserved? | l | log10(T_1/2) [s] |   T_1/2 " << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;
    for (Parity pp : {PLUS, MINUS}) {
        for (Parity pd : {PLUS, MINUS}) {
            perform_and_print_calculation(110, 273, 11.27, 0.5, pp, 0.5, pd);
        }
    }
    std::cout << "=======================================================================\n" << std::endl;
}

// Main interactive function for general use
void calculate_interactive() {
    int Z, Ap;
    double Qa, jp, jd;
    char exhaustive_choice, p_char;

    std::cout << "\n--- Interactive Alpha Half-Life Calculator ---" << std::endl;
    std::cout << "Enter Parent Z: "; std::cin >> Z;
    std::cout << "Enter Parent A: "; std::cin >> Ap;
    std::cout << "Enter Q_alpha (MeV): "; std::cin >> Qa;
    std::cout << "Enter Parent Spin (e.g., 5.5 for 11/2): "; std::cin >> jp;
    std::cout << "Enter Daughter Spin (e.g., 4.5 for 9/2): "; std::cin >> jd;
    
    std::cout << "\nPerform exhaustive parity search for all 4 combinations? (y/n): ";
    std::cin >> exhaustive_choice;

    std::cout << "\n-------------------------------------------------------------------" << std::endl;
    std::cout << "jp_pi -> jd_pi | Parity Conserved? | l | log10(T_1/2) [s] |   T_1/2 " << std::endl;
    std::cout << "-------------------------------------------------------------------" << std::endl;

    if (exhaustive_choice == 'y' || exhaustive_choice == 'Y') {
        for (Parity pp : {PLUS, MINUS}) {
            for (Parity pd : {PLUS, MINUS}) {
                perform_and_print_calculation(Z, Ap, Qa, jp, pp, jd, pd);
            }
        }
    } else {
        Parity pp = UNKNOWN, pd = UNKNOWN;
        std::cout << "Enter Parent Parity (+ or -): ";
        std::cin >> p_char;
        if (p_char == '+') pp = PLUS; else if (p_char == '-') pp = MINUS;

        std::cout << "Enter Daughter Parity (+ or -): ";
        std::cin >> p_char;
        if (p_char == '+') pd = PLUS; else if (p_char == '-') pd = MINUS;

        if (pp != UNKNOWN && pd != UNKNOWN) {
            perform_and_print_calculation(Z, Ap, Qa, jp, pp, jd, pd);
        } else {
            std::cout << "Invalid parity input. Calculation aborted." << std::endl;
        }
    }
    std::cout << "-------------------------------------------------------------------\n" << std::endl;
}

// Entry point of the script
void unified_halflife_calculator() {
    while (true) {
        int choice;
        std::cout << "========= Unified Alpha Half-Life Calculator Menu =========\n";
        std::cout << "1. Calculate half-life for a specific nucleus (Interactive Mode)\n";
        std::cout << "2. Run the pre-defined 273Ds example\n";
        std::cout << "3. Exit\n";
        std::cout << "===========================================================\n";
        std::cout << "Enter your choice: ";
        std::cin >> choice;

        if (std::cin.fail()) {
            std::cout << "Invalid input. Please enter a number." << std::endl;
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        
        switch (choice) {
            case 1:
                calculate_interactive();
                break;
            case 2:
                calculate_Ds273_example();
                break;
            case 3:
                std::cout << "Exiting." << std::endl;
                return;
            default:
                std::cout << "Invalid choice. Please try again." << std::endl;
                break;
        }
    }
}