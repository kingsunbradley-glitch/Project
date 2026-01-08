import math
import sys

# Physical Constants
C_LIGHT = 299792458.0          # Speed of light in m/s
U_TO_MEV = 931.494102          # Energy equivalent of 1 u in MeV

def calculate_brho_comparison(m_u, e_kev, q_e):
    """
    Calculate and compare Relativistic vs Classical Magnetic Rigidity.
    """
    # 1. Prepare constants
    # Convert mass to rest energy (keV)
    rest_mass_kev = m_u * U_TO_MEV * 1000.0
    
    # 2. Relativistic Calculation
    # (pc)^2 = E^2 + 2*E*m
    pc_rel_kev = math.sqrt(e_kev**2 + 2 * e_kev * rest_mass_kev)
    brho_rel = (pc_rel_kev * 1000.0) / (C_LIGHT * q_e)
    
    # 3. Classical (Non-relativistic) Calculation
    # (pc)^2 = 2*E*m (Ignoring E^2 term)
    pc_cls_kev = math.sqrt(2 * e_kev * rest_mass_kev)
    brho_cls = (pc_cls_kev * 1000.0) / (C_LIGHT * q_e)
    
    # 4. Calculate Deviation Percentage
    diff_percent = ((brho_cls - brho_rel) / brho_rel) * 100
    
    return brho_rel, brho_cls, diff_percent

def main():
    print("========================================")
    print("   Magnetic Rigidity Calculator (B-rho)")
    print("========================================")
    print("Input Format: q(e), m(u), E(keV)")
    print("Example: 2, 4, 5804.8  (for Alpha particle from 244Cm decay)")
    print("Tip: Type 'q' or 'exit' to quit.")
    
    while True:
        # Get input
        raw_input = input("\nEnter q, m, E (or 'q' to quit) >>> ").strip()
        
        # 1. Check for exit command
        if raw_input.lower() in ['q', 'quit', 'exit']:
            print("Looking forward to seeing you again !")
            break
        
        # 2. Check for empty input
        if not raw_input:
            continue

        try:
            # 3. Parse Data
            # Replace commas with spaces and split
            parts = raw_input.replace(',', ' ').split()
            
            if len(parts) != 3:
                raise ValueError("Expected 3 values.")
            
            # --- Input Order: q, m, E ---
            q = float(parts[0])  # Charge
            m = float(parts[1])  # Mass
            E = float(parts[2])  # Energy
            
            if q == 0:
                raise ValueError("Charge 'q' cannot be zero.")
            
            # 4. Execute Calculation (Function expects: m, E, q)
            rel, cls, diff = calculate_brho_comparison(m, E, q)
            
            # 5. Output Results
            print("-" * 45)
            print(f"Particle: q={q} e, m={m} u, E={E} keV")
            print("-" * 45)
            # Table Header
            print(f"{'Model':<15} | {'Rigidity B-rho (T·m)':<20}")
            print("-" * 45)
            # Table Rows
            print(f"{'Relativistic':<15} | {rel:.6f} T·m")
            print(f"{'Classical':<15} | {cls:.6f} T·m")
            print("-" * 45)
            
            # Deviation message
            msg = f"Deviation: {diff:+.4f}%"
            if abs(diff) > 1.0:
                msg += " [Significant deviation! Use Relativistic]"
            else:
                msg += " [Good approximation]"
            print(msg)
            
        except ValueError as e:
            print(f"Input Error: {e} -> Please try again.")
        except Exception as e:
            print(f"Unknown Error: {e}")

if __name__ == "__main__":
    main()