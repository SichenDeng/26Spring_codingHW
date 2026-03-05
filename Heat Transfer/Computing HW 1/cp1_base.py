"""
cp1_base.py
───────────────────────────────────────────────────────────────────────────────
530.334 Heat Transfer — Spring 2026  |  Computer Project 1

Satellite Thermal Design: Mercury Orbit
(Steady-state, 1D radial, nonlinear)

Numerical approach — FD marching + brentq:
  At steady state, Q = Q_int is constant through all solid layers
  (no internal heat sources in the solid; solar/Mercury IR act only at
  the outer surface). The 1D spherical conduction equation then reduces to:

      dT/dr = -Q_int / (4π r² k(T))

  This is a simple 1st-order ODE integrated numerically using the
  RK2 (midpoint) finite-difference method with N nodes per layer.

  Solution algorithm (outside-in):
    1. Solve outer radiation BC for Ts  (brentq — nonlinear in Ts)
    2. March dT/dr inward through Layer 3  →  T2p  at r2
    3. Apply contact resistance jump  →  T2m
    4. March inward through Layer 2   →  T1p  at r1
    5. Apply contact resistance jump  →  T1m
    6. March inward through Layer 1   →  Tw   at r0
    7. Solve convection equation for Tg  (brentq — nonlinear in Tg)

  Nonlinear k(T) is handled naturally — k is evaluated at the current
  temperature at each FD step; no linearization or outer iteration needed.
───────────────────────────────────────────────────────────────────────────────
"""

import numpy as np
from scipy.optimize import brentq


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1  —  GEOMETRY
# ═══════════════════════════════════════════════════════════════════════════════

R0 = 0.30        # inner cavity radius [m]
T1 = 0.002       # Layer 1 (Al liner) thickness           [m]
T3 = 0.003       # Layer 3 (CFRP outer shell) thickness   [m]

T2_MIN = 0.001   # aerogel thickness design lower bound   [m]
T2_MAX = 0.040   # aerogel thickness design upper bound   [m]


def get_radii(t2):
    """
    Return shell radii (r0, r1, r2, r3) for aerogel thickness t2 [m].

      r0 : inner cavity surface  (fixed)
      r1 = r0 + t1  →  Al/Aerogel interface
      r2 = r1 + t2  →  Aerogel/CFRP interface
      r3 = r2 + t3  →  outer surface
    """
    r1 = R0 + T1
    r2 = r1 + t2
    r3 = r2 + T3
    return R0, r1, r2, r3


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2  —  TEMPERATURE-DEPENDENT CONDUCTIVITIES  k(T)  [W/m·K]
# ═══════════════════════════════════════════════════════════════════════════════

def k1(T):
    """Layer 1 — Al alloy liner.  Linearly decreasing with T."""
    return 160.0 * (1.0 - 2.0e-4 * (T - 300.0))


def k2(T):
    """Layer 2 — Aerogel insulation.  Low-k, polynomial in (T − 300)."""
    dT = T - 300.0
    return 0.035 + 2.0e-4 * dT + 1.0e-7 * dT**2


def k3(T):
    """Layer 3 — CFRP outer shell.  Linearly increasing with T."""
    return 6.0 * (1.0 + 1.0e-4 * (T - 300.0))


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3  —  FD LAYER MARCHING  (RK2 / midpoint method)
#
#  Governing ODE inside each solid layer:
#      dT/dr = -Q / (4π r² k(T))
#
#  Integrated inward (r_start > r_end) with RK2 using N equally-spaced nodes.
#  k(T) is evaluated at the current (updated) temperature each step —
#  this is how nonlinear k(T) is handled without any outer iteration.
# ═══════════════════════════════════════════════════════════════════════════════

def _march_layer(r_start, r_end, T_start, k_func, Q, N=20):
    """
    Integrate  dT/dr = -Q / (4π r² k(T))  from r_start to r_end using RK2.

    r_start > r_end  →  marching inward; dr < 0; T increases. ✓

    Parameters
    ----------
    r_start, r_end : layer boundaries  [m]  (r_start is the outer/starting end)
    T_start        : temperature at r_start  [K]
    k_func         : conductivity function k(T)  [W/m·K]
    Q              : constant heat flow (positive = outward)  [W]
    N              : number of FD nodes including both endpoints (default 20)

    Returns
    -------
    r_arr : np.ndarray, shape (N,)  —  radii from r_start to r_end  [m]
    T_arr : np.ndarray, shape (N,)  —  temperatures at each node    [K]
    """
    r_arr = np.linspace(r_start, r_end, N)
    dr    = r_arr[1] - r_arr[0]   # negative when marching inward
    T_arr = np.zeros(N)
    T_arr[0] = T_start

    for i in range(N - 1):
        r = r_arr[i]
        T = T_arr[i]

        # Stage 1: slope at current node
        dTdr1 = -Q / (4.0 * np.pi * r**2 * k_func(T))

        # Stage 2: slope at midpoint (RK2 midpoint rule)
        r_m   = r + 0.5 * dr
        T_m   = T + 0.5 * dr * dTdr1
        dTdr2 = -Q / (4.0 * np.pi * r_m**2 * k_func(T_m))

        T_arr[i + 1] = T + dr * dTdr2

    return r_arr, T_arr


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4  —  CONTACT RESISTANCES  (area-specific, R'' [m²K/W])
#
#  Temperature jump at interface:  ΔT = Q · R''_c / (4π r²)
#  The inner (hotter) side is always higher by this amount.
# ═══════════════════════════════════════════════════════════════════════════════

RC12_PP = 2.0e-4   # at r = r1  (Al / Aerogel interface)   [m²K/W]
RC23_PP = 1.0e-4   # at r = r2  (Aerogel / CFRP interface)  [m²K/W]


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5  —  INTERNAL HEAT SOURCE
# ═══════════════════════════════════════════════════════════════════════════════

Q_INT_DEFAULT = 50.0   # electronics steady waste heat  [W]


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6  —  INNER BOUNDARY: NATURAL CONVECTION
#
#  Enclosure correlation:
#      Ra = g_eff β |Tg − Tw| Lc³ / (ν α)
#      Nu = 2 + 0.55 Ra^(1/4)
#      h_i = Nu · kg / Lc
# ═══════════════════════════════════════════════════════════════════════════════

KG        = 0.030         # cavity gas thermal conductivity  [W/m·K]
NU_GAS    = 1.7e-5        # kinematic viscosity              [m²/s]
ALPHA_GAS = 2.4e-5        # thermal diffusivity              [m²/s]
BETA      = 1.0 / 300.0   # thermal expansion coefficient   [1/K]
LC        = R0            # characteristic length = r0       [m]
G_EFF     = 0.02          # effective disturbance gravity    [m/s²]


def h_inner(Tg, Tw):
    """
    Inner-wall natural convection coefficient  h_i  [W/m²·K].
    Uses |Tg − Tw| so it works regardless of which side is hotter.
    """
    Ra = G_EFF * BETA * abs(Tg - Tw) * LC**3 / (NU_GAS * ALPHA_GAS)
    Nu = 2.0 + 0.55 * Ra**0.25
    return Nu * KG / LC


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7  —  OUTER BOUNDARY: RADIATION
#
#  Outer surface simultaneously:
#    • absorbs solar (shortwave)  → use α(T)
#    • absorbs Mercury IR         → use ε(T)  [thermal IR; Kirchhoff: α_IR = ε]
#    • emits to deep space        → use ε(T)
#
#  Net outward heat flow [W]:
#    Q_rad = 4π r3² [ ε σ (Ts⁴ − T_space⁴) − α g_sol G_sol − ε g_m G_m ]
#
#  At steady state:  Q_rad = Q_int  (solar/IR do not penetrate the cavity).
# ═══════════════════════════════════════════════════════════════════════════════

SIGMA      = 5.670e-8   # Stefan-Boltzmann constant           [W/m²·K⁴]
T_SPACE    = 3.0        # deep-space background temperature    [K]
G_SOL      = 9100.0     # solar irradiance at Mercury orbit    [W/m²]
G_SOL_FAC  = 0.25       # solar geometric factor
G_MERC     = 1200.0     # Mercury IR equivalent irradiance     [W/m²]
G_MERC_FAC = 0.25       # Mercury IR geometric factor


def alpha_s(T):
    """Outer surface solar absorptivity (shortwave)."""
    return 0.10 + 1.5e-4 * (T - 300.0)


def eps_s(T):
    """Outer surface thermal emissivity."""
    return 0.92 - 1.0e-4 * (T - 300.0)


def Q_rad_net(Ts, r3):
    """
    Net heat dissipated from outer surface by radiation  [W].
    Positive = net outward = must be balanced by inward conduction.
    """
    A = 4.0 * np.pi * r3**2
    a = alpha_s(Ts)
    e = eps_s(Ts)
    return A * (e * SIGMA * (Ts**4 - T_SPACE**4)
                - a * G_SOL_FAC  * G_SOL
                - e * G_MERC_FAC * G_MERC)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8  —  CORE STEADY-STATE SOLVER
# ═══════════════════════════════════════════════════════════════════════════════

def solve(t2, Q_int=Q_INT_DEFAULT, N=20):
    """
    Solve the complete steady-state 1D radial thermal model.

    Algorithm (outside-in, FD marching with RK2):
      1. Q = Q_int is constant through the shell at steady state.
      2. Solve outer radiation BC for Ts  (brentq).
      3. March dT/dr inward through Layer 3  →  T2p at r2.
         Apply contact resistance  →  T2m.
      4. March inward through Layer 2  →  T1p at r1.
         Apply contact resistance  →  T1m.
      5. March inward through Layer 1  →  Tw at r0.
      6. Solve convection equation for Tg  (brentq).

    Parameters
    ----------
    t2    : aerogel insulation thickness  [m]
    Q_int : electronics waste heat        [W]  (default 50 W)
    N     : FD nodes per layer            (default 20)

    Returns
    -------
    dict with keys:
      Tg, Tw          : cavity gas / inner-wall temperatures        [K]
      T1m, T1p        : temp at r1 from Al side / Aerogel side      [K]
      T2m, T2p        : temp at r2 from Aerogel side / CFRP side    [K]
      Ts              : outer surface temperature                    [K]
      Q               : heat flow through shell                      [W]
      dTc12, dTc23    : contact-resistance temperature jumps        [K]
      r0..r3          : shell radii                                  [m]
      r_profile       : radii array for T(r) plot, r0 → r3          [m]
      T_profile       : temperatures at each profile point          [K]
                        (interface r values are duplicated to show jumps)
    """
    r0, r1, r2, r3 = get_radii(t2)
    Q = Q_int   # constant heat flow = electronics source (steady state)

    # ── Step 1: outer surface temperature from radiation balance ──────────────
    Ts = brentq(lambda T: Q_rad_net(T, r3) - Q, 100.0, 1500.0)

    # ── Step 2: march inward through Layer 3 (CFRP) ──────────────────────────
    r_L3, T_L3 = _march_layer(r3, r2, Ts, k3, Q, N)
    T2p = T_L3[-1]   # temp at r2, CFRP side

    # Contact resistance at r2 (Aerogel side is hotter)
    dTc23 = Q * RC23_PP / (4.0 * np.pi * r2**2)
    T2m   = T2p + dTc23

    # ── Step 3: march inward through Layer 2 (Aerogel) ───────────────────────
    r_L2, T_L2 = _march_layer(r2, r1, T2m, k2, Q, N)
    T1p = T_L2[-1]   # temp at r1, Aerogel side

    # Contact resistance at r1 (Al side is hotter)
    dTc12 = Q * RC12_PP / (4.0 * np.pi * r1**2)
    T1m   = T1p + dTc12

    # ── Step 4: march inward through Layer 1 (Al) ────────────────────────────
    r_L1, T_L1 = _march_layer(r1, r0, T1m, k1, Q, N)
    Tw = T_L1[-1]   # inner wall temperature T(r0)

    # ── Step 5: cavity gas temperature from convection balance ────────────────
    # Q_int = h_i(Tg, Tw) · 4π r0² · (Tg − Tw)   →  solve for Tg
    Tg = brentq(
        lambda T: h_inner(T, Tw) * 4.0 * np.pi * r0**2 * (T - Tw) - Q_int,
        Tw + 1e-3, Tw + 2000.0
    )

    # ── Build outward-ordered profile (r0 → r3) ───────────────────────────────
    # Each layer was marched inward → reverse to outward order.
    # At r1: L1 ends at T1m, L2 starts at T1p  (contact jump visible in plot)
    # At r2: L2 ends at T2m, L3 starts at T2p  (contact jump visible in plot)
    r_profile = np.concatenate([r_L1[::-1], r_L2[::-1], r_L3[::-1]])
    T_profile = np.concatenate([T_L1[::-1], T_L2[::-1], T_L3[::-1]])

    return {
        'Tg': Tg,  'Tw': Tw,
        'T1m': T1m, 'T1p': T1p,
        'T2m': T2m, 'T2p': T2p,
        'Ts': Ts,  'Q': Q,
        'dTc12': dTc12, 'dTc23': dTc23,
        'r0': r0, 'r1': r1, 'r2': r2, 'r3': r3,
        't2': t2, 'Q_int': Q_int,
        'r_profile': r_profile,
        'T_profile': T_profile,
        'N': N,
    }


def temperature_profile(sol):
    """Return the pre-computed (r, T) profile arrays from a solve() result."""
    return sol['r_profile'], sol['T_profile']


# ═══════════════════════════════════════════════════════════════════════════════
# Quick sanity check — run with:  python cp1_base.py
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    sol = solve(t2=0.010)   # baseline: t2 = 10 mm

    print("=" * 50)
    print("Baseline check  (t2 = 10 mm, Q_int = 50 W)")
    print("=" * 50)
    print(f"  Tg    = {sol['Tg']:.3f} K   (limit ≤ 360 K)  "
          f"{'PASS ✓' if sol['Tg'] <= 360 else 'FAIL ✗'}")
    print(f"  Tw    = {sol['Tw']:.3f} K")
    print(f"  Ts    = {sol['Ts']:.3f} K   (limit ≤ 450 K)  "
          f"{'PASS ✓' if sol['Ts'] <= 450 else 'FAIL ✗'}")
    print(f"  Q     = {sol['Q']:.1f} W")
    print(f"  ΔTc12 = {sol['dTc12']:.4f} K   (T(r1⁺) − T(r1⁻))")
    print(f"  ΔTc23 = {sol['dTc23']:.4f} K   (T(r2⁺) − T(r2⁻))")
    print(f"  Radii: r0={sol['r0']:.4f}, r1={sol['r1']:.4f}, "
          f"r2={sol['r2']:.4f}, r3={sol['r3']:.4f} m")
