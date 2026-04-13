import numpy as np

# Water properties from Incropera Table A.6 (saturated liquid)
# Columns: T(°C), rho(kg/m³), cp(J/kg·K), mu(N·s/m²), k(W/m·K), Pr
_prop_table = np.array([
    [0,    999.8,   4217,   1.792e-3,   0.561,   13.44],
    [5,    999.9,   4205,   1.519e-3,   0.571,   11.19],
    [10,   999.7,   4194,   1.307e-3,   0.580,    9.44],
    [15,   999.1,   4186,   1.138e-3,   0.589,    8.09],
    [20,   998.0,   4182,   1.002e-3,   0.598,    7.01],
    [25,   997.0,   4180,   0.891e-3,   0.607,    6.14],
    [30,   995.6,   4178,   0.798e-3,   0.615,    5.42],
    [35,   994.0,   4178,   0.720e-3,   0.623,    4.83],
    [40,   992.1,   4179,   0.653e-3,   0.631,    4.33],
    [45,   990.1,   4180,   0.596e-3,   0.637,    3.91],
    [50,   988.0,   4181,   0.547e-3,   0.644,    3.55],
])


def water_props(T_celsius):
    """Return (rho, cp, mu, k, Pr) at given temperature via linear interpolation."""
    T = _prop_table[:, 0]
    return tuple(np.interp(T_celsius, T, _prop_table[:, i]) for i in range(1, 6))


# ---- Inputs ----
T_h_i = float(input("Hot water inlet temp T_h,i (°C): "))
T_c_i = float(input("Cold water inlet temp T_c,i (°C): "))
V_dot = float(input("Flow rate (L/min): "))
V_dot_m3s = V_dot / 1000.0 / 60.0   # m³/s

# ---- Geometry ----
# Copper tube (helical coil)
d_o = 0.25 * 0.0254       # OD: 1/4" -> m
d_i = 0.19 * 0.0254       # ID: ~0.19" -> m
L_tube = 10.0 * 0.3048    # Total length: 10 ft -> m
k_copper = 385.0           # W/(m·K)

# Coil
D_coil_inner = 2.5 * 0.0254       # Inner coil diameter: 2.5" -> m
D_coil = D_coil_inner + d_o       # Centerline diameter
N_turns = 24

# Shell (PVC)
D_shell = 4.0 * 0.0254    # ID: 4" -> m
L_shell = 10.0 * 0.0254   # Length: 10" -> m

# Baffles at 1/3 and 2/3 positions
N_baffles = 2
baffle_spacing = L_shell / (N_baffles + 1)

pitch = L_shell / N_turns          # Average pitch per turn
A_o = np.pi * d_o * L_tube        # Outer surface area
A_i = np.pi * d_i * L_tube        # Inner surface area
A_tube_cs = np.pi / 4 * d_i**2    # Tube cross-section

# ---- Iterative solution (properties depend on unknown outlet temps) ----
T_h_o_guess = T_h_i - 5.0
T_c_o_guess = T_c_i + 5.0

for iteration in range(20):
    T_h_mean = (T_h_i + T_h_o_guess) / 2.0
    T_c_mean = (T_c_i + T_c_o_guess) / 2.0

    rho_h, cp_h, mu_h, k_h, Pr_h = water_props(T_h_mean)
    rho_c, cp_c, mu_c, k_c, Pr_c = water_props(T_c_mean)

    m_dot_h = rho_h * V_dot_m3s
    m_dot_c = rho_c * V_dot_m3s

    # ---- Tube-side h_i (hot water in coil) ----
    v_tube = V_dot_m3s / A_tube_cs
    Re_tube = rho_h * v_tube * d_i / mu_h
    De = Re_tube * np.sqrt(d_i / D_coil)    # Dean number

    if Re_tube < 2300:
        # Laminar: Schmidt / Manlapaz-Churchill correlations
        if De > 100:
            Nu_tube = 0.036 * De**0.5 * Pr_h**(1.0/3.0)
        elif De > 11.6:
            Nu_tube = (3.66**3 + 1.158 * (De * Pr_h**0.5)**1.5)**(1.0/3.0)
        else:
            Nu_tube = 3.66
    else:
        # Turbulent: Gnielinski (Eq 8.60) with coil enhancement (Schmidt)
        f_s = (0.790 * np.log(Re_tube) - 1.64)**(-2)
        Nu_straight = ((f_s / 8) * (Re_tube - 1000) * Pr_h /
                       (1 + 12.7 * (f_s / 8)**0.5 * (Pr_h**(2/3) - 1)))
        ratio = d_i / D_coil
        Nu_tube = Nu_straight * (1 + 3.6 * (1 - ratio) * ratio**0.8)

    h_i = Nu_tube * k_h / d_i

    # ---- Shell-side h_o (cold water, Kern method) ----
    turns_per_section = N_turns / (N_baffles + 1)
    A_cross = 0.05 * baffle_spacing * (D_shell - d_o)
    v_shell = 3 * V_dot_m3s / A_cross

    # Hydraulic diameter
    A_shell_cs = np.pi / 4 * D_shell**2 - turns_per_section * np.pi / 4 * d_o**2
    P_wetted = np.pi * D_shell + turns_per_section * np.pi * d_o
    D_h_shell = 4.0 * A_shell_cs / P_wetted

    Re_shell = rho_c * v_shell * D_h_shell / mu_c

    # Zukauskas correlation for flow over tube bank
    if Re_shell < 40:
        C_z, m_z = 0.75, 0.4
    elif Re_shell < 1000:
        C_z, m_z = 0.51, 0.5
    elif Re_shell < 200000:
        C_z, m_z = 0.27, 0.63
    else:
        C_z, m_z = 0.021, 0.84

    Nu_shell = C_z * Re_shell**m_z * Pr_c**0.36
    h_o = Nu_shell * k_c / D_h_shell

    # ---- Overall U based on outer area (Eq 11.1a, cylindrical wall Eq 3.33) ----
    U_o = 1.0 / (1/h_o + d_o * np.log(d_o/d_i) / (2*k_copper) + d_o/(h_i*d_i))

    # ---- epsilon-NTU, counterflow (Eq 11.29a) ----
    C_h = m_dot_h * cp_h
    C_c = m_dot_c * cp_c
    C_min = min(C_h, C_c)
    C_max = max(C_h, C_c)
    C_r = C_min / C_max

    NTU = U_o * A_o / C_min

    if abs(C_r - 1.0) < 1e-6:
        epsilon = NTU / (1 + NTU)
    else:
        epsilon = ((1 - np.exp(-NTU * (1 - C_r))) /
                   (1 - C_r * np.exp(-NTU * (1 - C_r))))

    q = epsilon * C_min * (T_h_i - T_c_i)
    T_h_o = T_h_i - q / C_h
    T_c_o = T_c_i + q / C_c

    if abs(T_h_o - T_h_o_guess) < 0.001 and abs(T_c_o - T_c_o_guess) < 0.001:
        break
    T_h_o_guess = T_h_o
    T_c_o_guess = T_c_o

# ---- Heat transfer rates ----
q_h = m_dot_h * cp_h * (T_h_i - T_h_o)
q_c = m_dot_c * cp_c * (T_c_o - T_c_i)
eta = q_c / q_h

# ---- Pressure drop: tube side (Darcy-Weisbach + coil correction) ----
if Re_tube < 2300:
    f_tube_straight = 64.0 / Re_tube
else:
    f_tube_straight = (0.790 * np.log(Re_tube) - 1.64)**(-2)

# Misra & Gupta coil friction enhancement
f_coil_ratio = 1 + 0.033 * (np.log10(De))**4 if De > 30 else 1.0
f_tube = f_tube_straight * f_coil_ratio
dP_hot = 0.5*f_tube * (L_tube / d_i) * (rho_h * v_tube**2 / 2)
dP_hot_kPa = dP_hot / 1000

# ---- Pressure drop: shell side (improved: friction + baffle losses) ----

# friction factor
f_shell = 64.0 / Re_shell if Re_shell < 2300 else (0.790 * np.log(Re_shell) - 1.64)**(-2)

# effective flow length (zig-zag due to baffles)
L_eff = (N_baffles + 1) * baffle_spacing

# frictional loss
dP_fric = f_shell * (L_eff / D_h_shell) * (rho_c * v_shell**2 / 2)

# baffle-induced minor losses
K_b = 1.0   # 推荐 0.5–2，可以调
dP_baffle = N_baffles * K_b * (rho_c * v_shell**2 / 2)

# entrance/exit + leak + turning losses（关键！！）
K_extra = 2.0   # 推荐 1–5，根据实验调
dP_extra = K_extra * (rho_c * v_shell**2 / 2)

# total shell-side pressure drop
dP_cold = dP_fric + dP_baffle + dP_extra
dP_cold_kPa = dP_cold / 1000

# ---- Output ----
print("=" * 60)
print("  HEAT EXCHANGER PERFORMANCE MODEL")
print("=" * 60)
print()
print("--- Operating Conditions ---")
print(f"  T_h,i = {T_h_i:.2f} °C    T_c,i = {T_c_i:.2f} °C")
print(f"  Flow rate = {V_dot:.1f} L/min (both sides)")
print()
print("--- Geometry ---")
print(f"  Tube: OD={d_o*1000:.2f}mm  ID={d_i*1000:.2f}mm  L={L_tube:.3f}m")
print(f"  Coil: D_cl={D_coil*1000:.1f}mm  {N_turns} turns  pitch={pitch*1000:.2f}mm")
print(f"  Shell: ID={D_shell*1000:.1f}mm  L={L_shell*1000:.1f}mm  {N_baffles} baffles")
print(f"  A_o = {A_o:.4f} m²")
print()
print("--- Intermediate ---")
print(f"  Re_tube={Re_tube:.0f}  De={De:.1f}  Nu_tube={Nu_tube:.2f}  h_i={h_i:.1f} W/(m²·K)")
print(f"  Re_shell={Re_shell:.0f}  Nu_shell={Nu_shell:.2f}  h_o={h_o:.1f} W/(m²·K)")
print(f"  U_o={U_o:.1f} W/(m²·K)  C_r={C_r:.4f}")
print()
print("=" * 60)
print("  RESULTS")
print("=" * 60)
print(f"  T_h,o  = {T_h_o:.2f} °C")
print(f"  T_c,o  = {T_c_o:.2f} °C")
print(f"  ΔP_hot  = {dP_hot_kPa:.6f} kPa")
print(f"  ΔP_cold = {dP_cold_kPa:.6f} kPa")
print(f"  q_h→   = {q_h:.2f} W  (= F1)")
print(f"  q→c    = {q_c:.2f} W")
print(f"  η      = {eta:.4f}")
print(f"  U      = {U_o:.2f} W/(m²·K)")
print(f"  ε      = {epsilon:.4f}")
print(f"  NTU    = {NTU:.4f}")
print("=" * 60)
