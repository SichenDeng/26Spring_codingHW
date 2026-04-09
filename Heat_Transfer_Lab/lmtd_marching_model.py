"""
LMTD Finite-Difference Marching Model (counterflow)
Reproduces the MATLAB HeatExchangerGUI demo in Python.

Divides the heat exchanger into N segments. For each segment, uses
bisection to solve the local LMTD energy balance. Plots T_h and T_c
along the exchanger length.
"""

import numpy as np
import matplotlib.pyplot as plt

# ---- Water properties (linear interpolation from Incropera Table A.6) ----
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


def counterflow_lmtd_march(mhdot, mcdot, ch, cc, Thi, Tci, U, L, D, N=10000):
    """
    Counterflow heat exchanger solved by LMTD marching (same as MATLAB demo).

    Parameters
    ----------
    mhdot, mcdot : mass flow rates (kg/s)
    ch, cc       : specific heats (J/kg·K)
    Thi, Tci     : inlet temperatures (°C)
    U            : overall heat transfer coefficient (W/m²·K)
    L            : exchanger length (m)
    D            : mean tube diameter for heat transfer area = π·D·L (m)
    N            : number of segments

    Returns
    -------
    Th, Tc : temperature arrays along length (N+1 points)
    Th_out, Tc_out, epsilon
    """
    dx = L / N
    Th = np.zeros(N + 1)
    Tc = np.zeros(N + 1)
    Th[0] = Thi
    Tc[N] = Tci  # counterflow: cold enters at x=L

    if mhdot * ch == mcdot * cc:
        # Special case: C_h == C_c → ΔT constant along length
        Tco = (Tci + U * np.pi * D * L * Thi / (mcdot * cc)) / \
              (1 + U * np.pi * D * L / (mcdot * cc))
        Tho = Thi + mcdot * cc * (Tci - Tco) / (mhdot * ch)
        Tc[0] = Tco
        Th[N] = Tho
        for i in range(1, N):
            Tc[i] = Tco + U * np.pi * D * (i * dx) * (Tco - Thi) / (mcdot * cc)
            Th[i] = Thi + mcdot * cc * (Tc[i] - Tco) / (mhdot * ch)
    else:
        # Step 1: find overall Tco by bisection on full-length LMTD balance
        Tco1, Tco2 = Tci, Thi
        Fold = 0.0
        for _ in range(200):
            Tco = (Tco1 + Tco2) / 2.0
            Tho = Thi + (mcdot * cc) * (Tci - Tco) / (mhdot * ch)
            DT1 = Thi - Tco
            DT2 = Tho - Tci
            if abs(DT1 - DT2) < 1e-12:
                Tlm = DT1
            else:
                Tlm = (DT1 - DT2) / np.log(DT1 / DT2)
            F = mcdot * cc * (Tco - Tci) - U * np.pi * D * L * Tlm
            if F == Fold:
                break
            if F < 0:
                Tco1 = Tco
            else:
                Tco2 = Tco
            Fold = F

        Tc[0] = Tco
        Th[N] = Tho

        # Step 2: march through segments
        if mhdot * ch < mcdot * cc:
            # March from hot inlet (i=0) to hot outlet (i=N)
            _Thi = Thi
            _Tco = Tco
            for i in range(N):
                Tci1 = Tc[N]
                Tci2 = _Tco
                Fold = 0.0
                for _ in range(200):
                    _Tci = (Tci1 + Tci2) / 2.0
                    _Tho = _Thi + (mcdot * cc) * (_Tci - _Tco) / (mhdot * ch)
                    DT1 = _Thi - _Tco
                    DT2 = _Tho - _Tci
                    if abs(DT1 - DT2) < 1e-12:
                        Tlm = DT1
                    else:
                        if DT1 / DT2 <= 0:
                            break
                        Tlm = (DT1 - DT2) / np.log(DT1 / DT2)
                    F = mcdot * cc * (_Tco - _Tci) - U * np.pi * D * dx * Tlm
                    if F == Fold or np.isnan(F):
                        break
                    if F > 0:
                        Tci1 = _Tci
                    else:
                        Tci2 = _Tci
                    Fold = F
                Th[i + 1] = _Tho
                Tc[i + 1] = _Tci
                _Thi = _Tho
                _Tco = _Tci
        else:
            # March from cold inlet (i=N) to cold outlet (i=0)
            _Tho = Tho
            _Tci = Tci
            for i in range(N, 0, -1):
                Thi1 = _Tho
                Thi2 = Th[0]
                Fold = 0.0
                for _ in range(200):
                    _Thi = (Thi1 + Thi2) / 2.0
                    _Tco = _Tci + (mhdot * ch) * (_Thi - _Tho) / (mcdot * cc)
                    DT1 = _Thi - _Tco
                    DT2 = _Tho - _Tci
                    if abs(DT1 - DT2) < 1e-12:
                        Tlm = DT1
                    else:
                        if DT1 / DT2 <= 0:
                            break
                        Tlm = (DT1 - DT2) / np.log(DT1 / DT2)
                    F = mhdot * ch * (_Thi - _Tho) - U * np.pi * D * dx * Tlm
                    if F == Fold or np.isnan(F):
                        break
                    if F < 0:
                        Thi1 = _Thi
                    else:
                        Thi2 = _Thi
                    Fold = F
                Th[i - 1] = _Thi
                Tc[i - 1] = _Tco
                _Tho = _Thi
                _Tci = _Tco

    Th_out = Th[-1]
    Tc_out = Tc[0]

    # Effectiveness
    qh = (mhdot * ch) * (Th[0] - Th_out)
    Cmin_side = mhdot * ch if mhdot * ch < mcdot * cc else mcdot * cc
    qmax = Cmin_side * (Th[0] - Tc[N])
    epsilon = qh / qmax

    return Th, Tc, Th_out, Tc_out, epsilon


# ========== Use our heat exchanger parameters ==========
T_h_i = 30.0
T_c_i = 5.0
V_dot = 2.0  # L/min
V_dot_m3s = V_dot / 1000.0 / 60.0

# Geometry
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

# Mean tube diameter for area calculation (same as epsilon-NTU model: A_o = π·d_o·L)
D_mean = d_o

# Evaluate properties at estimated mean temperatures
T_h_mean = (T_h_i + (T_h_i - 5)) / 2.0  # ~27.5
T_c_mean = (T_c_i + (T_c_i + 5)) / 2.0  # ~7.5

rho_h, cp_h, mu_h, k_h, Pr_h = water_props(T_h_mean)
rho_c, cp_c, mu_c, k_c, Pr_c = water_props(T_c_mean)

m_dot_h = rho_h * V_dot_m3s
m_dot_c = rho_c * V_dot_m3s

# Compute U_o the same way as our main model (copy the logic)
A_tube_cs = np.pi / 4 * d_i**2
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
h_o = Nu_shell * k_c / D_h_shell

U_o = 1.0 / (1/h_o + d_o * np.log(d_o/d_i) / (2*k_copper) + d_o/(h_i*d_i))

# ---- Run LMTD marching model ----
Th, Tc, Th_out, Tc_out, epsilon = counterflow_lmtd_march(
    m_dot_h, m_dot_c, cp_h, cp_c, T_h_i, T_c_i, U_o, L_tube, D_mean
)

q_h = m_dot_h * cp_h * (T_h_i - Th_out)
q_c = m_dot_c * cp_c * (Tc_out - T_c_i)
C_h = m_dot_h * cp_h
C_c = m_dot_c * cp_c
C_min = min(C_h, C_c)
C_max = max(C_h, C_c)
NTU = U_o * np.pi * D_mean * L_tube / C_min

# ---- Output ----
print("=" * 60)
print("  LMTD MARCHING MODEL (counterflow)")
print("=" * 60)
print(f"  U_o    = {U_o:.2f} W/(m²·K)")
print(f"  T_h,o  = {Th_out:.2f} °C")
print(f"  T_c,o  = {Tc_out:.2f} °C")
print(f"  q_h    = {q_h:.2f} W")
print(f"  q_c    = {q_c:.2f} W")
print(f"  η      = {q_c/q_h:.4f}")
print(f"  ε      = {epsilon:.4f}")
print(f"  NTU    = {NTU:.4f}")
print("=" * 60)

# ---- Plot temperature profiles ----
x = np.linspace(0, L_tube, len(Th))
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x, Th, 'r', linewidth=2, label='$T_h$ (hot)')
ax.plot(x, Tc, 'b', linewidth=2, label='$T_c$ (cold)')
ax.set_xlabel('Length (m)', fontsize=14)
ax.set_ylabel('Temperature (°C)', fontsize=14)
ax.set_title('Counterflow Heat Exchanger — LMTD Marching')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.tick_params(labelsize=12)
plt.tight_layout()
plt.savefig('lmtd_temperature_profile.png', dpi=150)
plt.show()
print("Saved: lmtd_temperature_profile.png")
