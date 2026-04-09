"""
Notebook Method Model (epsilon-NTU, Eq 11.30a shell-and-tube)
Based on the lecture notebook approach, applied to our heat exchanger geometry.

Differences from our main model (heat_exchanger_model.py):
 - Effectiveness: Eq 11.30a (1 shell pass) instead of Eq 11.29a (counterflow)
 - Shell-side h_o: 0.6 correction factor for non-ideal shell flow
 - Pressure drop: includes minor losses (K_minor)
"""

import numpy as np

# ---- Water properties ----
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
    T = _prop_table[:, 0]
    return tuple(np.interp(T_celsius, T, _prop_table[:, i]) for i in range(1, 6))


def effectiveness_shell_and_tube(NTU, Cr):
    """Eq 11.30a: one shell pass, any even number of tube passes."""
    if abs(1 - Cr) < 1e-12:
        return NTU / (1 + NTU)
    sqrt_term = np.sqrt(1 + Cr**2)
    exp_term = np.exp(-NTU * sqrt_term)
    return 2.0 / (1 + Cr + sqrt_term * (1 + exp_term) / (1 - exp_term))


def effectiveness_counterflow(NTU, Cr):
    """Eq 11.29a: pure counterflow."""
    if abs(Cr - 1.0) < 1e-6:
        return NTU / (1 + NTU)
    return (1 - np.exp(-NTU * (1 - Cr))) / (1 - Cr * np.exp(-NTU * (1 - Cr)))


# ---- Inputs ----
T_h_i = 30.0
T_c_i = 5.0
V_dot = 2.0  # L/min
V_dot_m3s = V_dot / 1000.0 / 60.0

# ---- Geometry ----
d_o = 0.25 * 0.0254
d_i = 0.19 * 0.0254
L_tube = 10.0 * 0.3048
k_copper = 385.0
D_coil_inner = 2.5 * 0.0254
D_coil = D_coil_inner + d_o
N_turns = 24
D_shell = 4.0 * 0.0254
L_shell = 10.0 * 0.0254
N_baffles = 2
baffle_spacing = L_shell / (N_baffles + 1)
pitch = L_shell / N_turns
A_o = np.pi * d_o * L_tube
A_i = np.pi * d_i * L_tube
A_tube_cs = np.pi / 4 * d_i**2

# ---- Iterative solution ----
T_h_o_guess = T_h_i - 5.0
T_c_o_guess = T_c_i + 5.0

# Store results for both methods
results = {}

for iteration in range(20):
    T_h_mean = (T_h_i + T_h_o_guess) / 2.0
    T_c_mean = (T_c_i + T_c_o_guess) / 2.0

    rho_h, cp_h, mu_h, k_h, Pr_h = water_props(T_h_mean)
    rho_c, cp_c, mu_c, k_c, Pr_c = water_props(T_c_mean)

    m_dot_h = rho_h * V_dot_m3s
    m_dot_c = rho_c * V_dot_m3s

    # ---- Tube-side h_i (same as our model: Gnielinski + Dean number) ----
    v_tube = V_dot_m3s / A_tube_cs
    Re_tube = rho_h * v_tube * d_i / mu_h
    De = Re_tube * np.sqrt(d_i / D_coil)

    if Re_tube < 2300:
        if De > 100:
            Nu_tube = 0.036 * De**0.5 * Pr_h**(1.0/3.0)
        elif De > 11.6:
            Nu_tube = (3.66**3 + 1.158 * (De * Pr_h**0.5)**1.5)**(1.0/3.0)
        else:
            Nu_tube = 3.66
    else:
        f_s = (0.790 * np.log(Re_tube) - 1.64)**(-2)
        Nu_straight = ((f_s / 8) * (Re_tube - 1000) * Pr_h /
                       (1 + 12.7 * (f_s / 8)**0.5 * (Pr_h**(2/3) - 1)))
        ratio = d_i / D_coil
        Nu_tube = Nu_straight * (1 + 3.6 * (1 - ratio) * ratio**0.8)

    h_i = Nu_tube * k_h / d_i

    # ---- Shell-side h_o (Zukauskas + 0.6 correction from notebook) ----
    turns_per_section = N_turns / (N_baffles + 1)
    A_cross = baffle_spacing * (D_shell - d_o)
    v_shell = V_dot_m3s / A_cross
    A_shell_cs = np.pi / 4 * D_shell**2 - turns_per_section * np.pi / 4 * d_o**2
    P_wetted = np.pi * D_shell + turns_per_section * np.pi * d_o
    D_h_shell = 4.0 * A_shell_cs / P_wetted
    Re_shell = rho_c * v_shell * D_h_shell / mu_c

    if Re_shell < 40:
        C_z, m_z = 0.75, 0.4
    elif Re_shell < 1000:
        C_z, m_z = 0.51, 0.5
    elif Re_shell < 200000:
        C_z, m_z = 0.27, 0.63
    else:
        C_z, m_z = 0.021, 0.84

    Nu_shell = C_z * Re_shell**m_z * Pr_c**0.36
    h_o_ideal = Nu_shell * k_c / D_h_shell
    h_o = 0.6 * h_o_ideal   # notebook correction factor

    # ---- U_o ----
    U_o = 1.0 / (1/h_o + d_o * np.log(d_o/d_i) / (2*k_copper) + d_o/(h_i*d_i))

    # ---- epsilon-NTU ----
    C_h = m_dot_h * cp_h
    C_c = m_dot_c * cp_c
    C_min = min(C_h, C_c)
    C_max = max(C_h, C_c)
    C_r = C_min / C_max
    NTU = U_o * A_o / C_min

    # Notebook method: Eq 11.30a (shell-and-tube)
    eps_st = effectiveness_shell_and_tube(NTU, C_r)
    # Our method: Eq 11.29a (counterflow)
    eps_cf = effectiveness_counterflow(NTU, C_r)

    # Use shell-and-tube effectiveness for this model
    epsilon = eps_st
    q = epsilon * C_min * (T_h_i - T_c_i)
    T_h_o = T_h_i - q / C_h
    T_c_o = T_c_i + q / C_c

    if abs(T_h_o - T_h_o_guess) < 0.001 and abs(T_c_o - T_c_o_guess) < 0.001:
        break
    T_h_o_guess = T_h_o
    T_c_o_guess = T_c_o

# ---- Pressure drop: tube side (with minor losses like notebook) ----
if Re_tube < 2300:
    f_tube_straight = 64.0 / Re_tube
else:
    f_tube_straight = (0.790 * np.log(Re_tube) - 1.64)**(-2)

f_coil_ratio = 1 + 0.033 * (np.log10(De))**4 if De > 30 else 1.0
f_tube = f_tube_straight * f_coil_ratio
K_minor = 5.0  # entry/exit + bends (from notebook)
dP_hot = f_tube * (L_tube / d_i) * (rho_h * v_tube**2 / 2) + K_minor * (rho_h * v_tube**2 / 2)
dP_hot_kPa = dP_hot / 1000

# ---- Pressure drop: shell side ----
f_shell = 64.0 / Re_shell if Re_shell < 2300 else (0.790 * np.log(Re_shell) - 1.64)**(-2)
dP_cold = f_shell * (D_shell / D_h_shell) * (N_baffles + 1) * (rho_c * v_shell**2 / 2)
dP_cold_kPa = dP_cold / 1000

# Also compute counterflow result for comparison
q_cf = eps_cf * C_min * (T_h_i - T_c_i)
T_h_o_cf = T_h_i - q_cf / C_h
T_c_o_cf = T_c_i + q_cf / C_c

# ---- Output ----
print("=" * 65)
print("  COMPARISON: Notebook Method vs Our Model")
print("=" * 65)
print()
print(f"  {'':30s} {'Notebook':>12s} {'Our Model':>12s}")
print(f"  {'':30s} {'(11.30a)':>12s} {'(11.29a)':>12s}")
print(f"  {'-'*54}")
print(f"  {'h_o (with 0.6 correction)':30s} {h_o:>12.1f}  {'—':>10s}")
print(f"  {'h_o (no correction)':30s} {'—':>12s}  {h_o_ideal:>10.1f}")
print(f"  {'U_o [W/(m²·K)]':30s} {U_o:>12.2f}  {'—':>10s}")
print(f"  {'NTU':30s} {NTU:>12.4f}  {NTU:>10.4f}")
print(f"  {'C_r':30s} {C_r:>12.4f}  {C_r:>10.4f}")
print(f"  {'epsilon':30s} {eps_st:>12.4f}  {eps_cf:>10.4f}")
print(f"  {'T_h,o [°C]':30s} {T_h_o:>12.2f}  {T_h_o_cf:>10.2f}")
print(f"  {'T_c,o [°C]':30s} {T_c_o:>12.2f}  {T_c_o_cf:>10.2f}")
print(f"  {'q [W]':30s} {q:>12.2f}  {q_cf:>10.2f}")
print(f"  {'ΔP_hot [kPa] (with K=5)':30s} {dP_hot_kPa:>12.3f}  {'—':>10s}")
print()

# Also print what our original model gives (no 0.6 correction, Eq 11.29a)
U_o_orig = 1.0 / (1/h_o_ideal + d_o * np.log(d_o/d_i) / (2*k_copper) + d_o/(h_i*d_i))
NTU_orig = U_o_orig * A_o / C_min
eps_orig = effectiveness_counterflow(NTU_orig, C_r)
q_orig = eps_orig * C_min * (T_h_i - T_c_i)

print(f"  For reference — our original model (no 0.6, Eq 11.29a):")
print(f"    U_o = {U_o_orig:.2f}, eps = {eps_orig:.4f}, q = {q_orig:.2f} W")
print("=" * 65)
