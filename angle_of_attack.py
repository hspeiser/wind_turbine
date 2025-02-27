import math
import numpy as np
import matplotlib.pyplot as plt

def blade_twist_distribution(R=0.0762,
                             alpha_opt_deg=6.0,
                             lambda_val=5.0,
                             root_fraction=0.3,
                             n_stations=6):
    """
    Compute a quick 'ideal' twist distribution from root_fraction*R to R
    ignoring induction.
    """
    r_root = root_fraction * R
    r_vals = np.linspace(r_root, R, n_stations)
    
    alpha_opt_rad = math.radians(alpha_opt_deg)

    beta_deg_list = []
    for r in r_vals:
        phi_rad = math.atan(R / (lambda_val*r))  # inflow angle
        phi_deg = math.degrees(phi_rad)
        beta_deg = phi_deg - alpha_opt_deg
        beta_deg_list.append(beta_deg)

    return r_vals, np.array(beta_deg_list)

# --- Example usage ---
if __name__ == "__main__":
    R_m = 0.0762          # 3 inches in meters
    alpha_opt_deg = 6.0   # best L/D angle of attack
    lambda_val     = 5.0  # chosen tip speed ratio
    root_fraction  = 0.3
    n_stations     = 6

    r_vals, twist_deg = blade_twist_distribution(R=R_m,
                                                 alpha_opt_deg=alpha_opt_deg,
                                                 lambda_val=lambda_val,
                                                 root_fraction=root_fraction,
                                                 n_stations=n_stations)

    # Print table
    print(" r [m]    Twist [deg]")
    print("---------------------")
    for rr, tt in zip(r_vals, twist_deg):
        print(f"{rr:7.4f}    {tt:7.3f}")

    # Plot
    plt.figure()
    plt.plot(r_vals*1000.0, twist_deg, 'ro-')
    plt.xlabel("Radius from hub [mm]")
    plt.ylabel("Twist angle [deg]")
    plt.title("Blade Twist Distribution")
    plt.grid(True)
    plt.show()
