"""
cp1_base.py — 530.334 Heat Transfer, CP1
Satellite thermal model: 1D radial steady-state solver.
"""

import numpy as np
from scipy.optimize import brentq

# Geometry
R0 = 0.30       # inner cavity radius [m]
T1 = 0.002      # Layer 1 Al liner thickness [m]
T3 = 0.003      # Layer 3 CFRP shell thickness [m]
T2_MIN = 0.001  # aerogel thickness bounds [m]
T2_MAX = 0.040

def get_radii(t2):
    r1 = R0 + T1
    r2 = r1 + t2
    r3 = r2 + T3
    return R0, r1, r2, r3

# Temperature-dependent conductivities [W/m-K]
def k1(T):
    return 160.0 * (1.0 - 2.0e-4 * (T - 300.0))

def k2(T):
    dT = T - 300.0
    return 0.035 + 2.0e-4 * dT + 1.0e-7 * dT**2

def k3(T):
    return 6.0 * (1.0 + 1.0e-4 * (T - 300.0))

# FD marching: integrate dT/dr = -Q/(4*pi*r^2*k(T)) using RK2
def _march_layer(r_start, r_end, T_start, k_func, Q, N=20):
    r_arr = np.linspace(r_start, r_end, N)
    dr = r_arr[1] - r_arr[0]
    T_arr = np.zeros(N)
    T_arr[0] = T_start

    for i in range(N - 1):
        r, T = r_arr[i], T_arr[i]
        dTdr1 = -Q / (4.0 * np.pi * r**2 * k_func(T))
        r_m = r + 0.5 * dr
        T_m = T + 0.5 * dr * dTdr1
        dTdr2 = -Q / (4.0 * np.pi * r_m**2 * k_func(T_m))
        T_arr[i + 1] = T + dr * dTdr2

    return r_arr, T_arr

# Contact resistances [m^2*K/W]
RC12_PP = 2.0e-4  # Al / Aerogel at r1
RC23_PP = 1.0e-4  # Aerogel / CFRP at r2

Q_INT_DEFAULT = 50.0  # electronics waste heat [W]

# Inner boundary: natural convection
KG        = 0.030
NU_GAS    = 1.7e-5
ALPHA_GAS = 2.4e-5
BETA      = 1.0 / 300.0
LC        = R0
G_EFF     = 0.02

def h_inner(Tg, Tw):
    Ra = G_EFF * BETA * abs(Tg - Tw) * LC**3 / (NU_GAS * ALPHA_GAS)
    Nu = 2.0 + 0.55 * Ra**0.25
    return Nu * KG / LC

# Outer boundary: radiation
SIGMA      = 5.670e-8
T_SPACE    = 3.0
G_SOL      = 9100.0
G_SOL_FAC  = 0.25
G_MERC     = 1200.0
G_MERC_FAC = 0.25

def alpha_s(T):
    return 0.10 + 1.5e-4 * (T - 300.0)

def eps_s(T):
    return 0.92 - 1.0e-4 * (T - 300.0)

def Q_rad_net(Ts, r3):
    A = 4.0 * np.pi * r3**2
    a, e = alpha_s(Ts), eps_s(Ts)
    return A * (e * SIGMA * (Ts**4 - T_SPACE**4)
                - a * G_SOL_FAC * G_SOL
                - e * G_MERC_FAC * G_MERC)

# Core solver: march from outer surface inward
def solve(t2, Q_int=Q_INT_DEFAULT, N=20):
    r0, r1, r2, r3 = get_radii(t2)
    Q = Q_int

    # Outer surface temp from radiation balance
    Ts = brentq(lambda T: Q_rad_net(T, r3) - Q, 100.0, 1500.0)

    # Layer 3 (CFRP) inward
    r_L3, T_L3 = _march_layer(r3, r2, Ts, k3, Q, N)
    T2p = T_L3[-1]

    # Contact resistance at r2
    dTc23 = Q * RC23_PP / (4.0 * np.pi * r2**2)
    T2m = T2p + dTc23

    # Layer 2 (Aerogel) inward
    r_L2, T_L2 = _march_layer(r2, r1, T2m, k2, Q, N)
    T1p = T_L2[-1]

    # Contact resistance at r1
    dTc12 = Q * RC12_PP / (4.0 * np.pi * r1**2)
    T1m = T1p + dTc12

    # Layer 1 (Al) inward
    r_L1, T_L1 = _march_layer(r1, r0, T1m, k1, Q, N)
    Tw = T_L1[-1]

    # Cavity gas temp from convection balance
    Tg = brentq(
        lambda T: h_inner(T, Tw) * 4.0 * np.pi * r0**2 * (T - Tw) - Q_int,
        Tw + 1e-3, Tw + 2000.0
    )

    # Build outward profile (reverse each layer)
    r_profile = np.concatenate([r_L1[::-1], r_L2[::-1], r_L3[::-1]])
    T_profile = np.concatenate([T_L1[::-1], T_L2[::-1], T_L3[::-1]])

    return {
        'Tg': Tg, 'Tw': Tw,
        'T1m': T1m, 'T1p': T1p,
        'T2m': T2m, 'T2p': T2p,
        'Ts': Ts, 'Q': Q,
        'dTc12': dTc12, 'dTc23': dTc23,
        'r0': r0, 'r1': r1, 'r2': r2, 'r3': r3,
        't2': t2, 'Q_int': Q_int,
        'r_profile': r_profile, 'T_profile': T_profile, 'N': N,
    }

def temperature_profile(sol):
    return sol['r_profile'], sol['T_profile']

if __name__ == '__main__':
    sol = solve(t2=0.010)
    print("Baseline (t2=10mm, Q_int=50W)")
    print(f"Tg = {sol['Tg']:.3f} K, {'PASS' if sol['Tg'] <= 360 else 'FAIL'}")
    print(f"Tw = {sol['Tw']:.3f} K")
    print(f"Ts = {sol['Ts']:.3f} K, {'PASS' if sol['Ts'] <= 450 else 'FAIL'}")
    print(f"Q = {sol['Q']:.1f} W")
    print(f"dTc12 = {sol['dTc12']:.4f} K")
    print(f"dTc23 = {sol['dTc23']:.4f} K")
