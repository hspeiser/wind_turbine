import math
import numpy as np

def mph_to_mps(mph):
    return mph * 0.44704

#----- A more realistic CP function that depends on both λ and α -----
def cp_more_realistic(lambda_val, alpha_deg):
    """
    Returns an approximate CP that peaks near lambda ~ 5-6
    and alpha ~ 6 degrees.
    We'll define:
       CP_base(lambda) = 0.45 * (lambda/6) * exp(1 - lambda/6)
    Then a penalty for alpha deviation from alpha_opt=6:
       alpha_penalty(alpha) = 0.01 * (alpha - 6)^2
    So if alpha=6, penalty=0, we get the full CP_base.
    If alpha=8, penalty=0.04 => CP is 0.41 if base was 0.45, etc.
    We'll clamp CP >= 0 if penalty overshoots.

    Adjust constants as desired for realism.
    """
    # Baseline CP vs. lambda (peak ~0.45 at λ~6)
    lam_base = 6.0
    cp_base = 0.45 * (lambda_val/lam_base) * math.exp(1 - (lambda_val/lam_base))
    if cp_base < 0:
        cp_base = 0

    # alpha penalty
    alpha_opt = 6.0
    diff = alpha_deg - alpha_opt
    penalty = 0.01 * (diff**2)  # 0.01*(α-6)^2

    cp = cp_base - penalty
    if cp < 0:
        cp = 0.0

    return cp

def mechanical_power_25mph(R_inch, lambda_val, alpha_deg,
                           induction_a=0.3, wind_mph=25.0):
    """
    Compute mechanical power for a single design point,
    using an induction factor a => rotor sees (1-a)*25mph
    and a CP function that depends on both lambda & alpha.
    Returns (P_mech_W, cp_used).
    """
    # Effective wind speed
    V_eff_mps = mph_to_mps(wind_mph) * (1.0 - induction_a)

    # Convert R_inch -> meters
    R_m = R_inch * 0.0254
    area = math.pi * (R_m**2)
    rho = 1.225  # kg/m^3

    # Power in wind at effective speed
    P_wind = 0.5 * rho * area * (V_eff_mps**3)

    # CP depends on λ & α
    cp_val = cp_more_realistic(lambda_val, alpha_deg)

    # mechanical power
    P_mech = cp_val * P_wind
    return P_mech, cp_val

def twist_distribution_induction(R_inch, alpha_deg, lambda_val,
                                 induction_a=0.3,
                                 root_fraction=0.3, n_stations=6):
    """
    Compute twist distribution from root_fraction*R_inch to R_inch,
    ignoring swirl induction but factoring the induction for 'effective speed.'
    We'll treat the rotor as if it runs at λ = (tip speed)/( (1-a)*V ).
    So local inflow angle is still φ(r)=arctan(R/(λ*r)), i.e. ignoring swirl.

    Then local pitch = φ(r) - alpha.
    We return arrays for radius [m] and pitch [deg].
    """
    import numpy as np

    R_m = R_inch * 0.0254
    r_root = root_fraction * R_m
    r_vals = np.linspace(r_root, R_m, n_stations)

    twist_deg = []
    alpha_opt = alpha_deg
    for rr in r_vals:
        # local inflow angle in degrees
        phi_deg = math.degrees( math.atan( R_m / (lambda_val * rr ) ) )
        # pitch = φ - alpha
        beta = phi_deg - alpha_opt
        twist_deg.append(beta)

    return r_vals, twist_deg

if __name__ == "__main__":

    # We'll do a parameter sweep over alpha=4..8, λ=4..7, R=2..3 in
    alpha_list = [4,5,6,7,8]
    lambda_list = [4,5,6,7]
    R_inch_vals = np.arange(2.0, 3.001, 0.25)  # e.g. 2.0,2.25,2.5,2.75,3.0

    # Induction factor
    a = 0.30

    results = []
    for alpha_deg in alpha_list:
        for lam in lambda_list:
            for r_in in R_inch_vals:
                p_mech, cp_val = mechanical_power_25mph(r_in, lam, alpha_deg,
                                                        induction_a=a)
                results.append({
                   'alpha': alpha_deg,
                   'lambda': lam,
                   'R_in': r_in,
                   'P_mech_W': p_mech,
                   'cp': cp_val
                })

    results_sorted = sorted(results, key=lambda x: x['P_mech_W'], reverse=True)

    print("TOP 10 combos:")
    for i in range(10):
        r = results_sorted[i]
        print(f"{i+1}) alpha={r['alpha']:.1f}°, λ={r['lambda']}, R={r['R_in']:.2f}\" => "
              f"P_mech={r['P_mech_W']:.2f} W, CP={r['cp']:.3f}")

    best = results_sorted[0]
    alpha_star  = best['alpha']
    lambda_star = best['lambda']
    R_star_in   = best['R_in']
    P_star_W    = best['P_mech_W']

    print("\nBEST param set =>", best)

    rvals_m, twist_deg = twist_distribution_induction(R_star_in, alpha_star, lambda_star,
                                                      induction_a=a)
    print("\nTwist distribution from ~30% radius to tip:")
    for (rr, td) in zip(rvals_m, twist_deg):
        print(f"  r={rr*1000:.1f} mm, pitch={td:.2f}°")
