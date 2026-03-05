"""
cp1_tasks.py
───────────────────────────────────────────────────────────────────────────────
530.334 Heat Transfer — Spring 2026  |  Computer Project 1

Tasks B through E — numerical results and plots.
All physics/solver logic lives in cp1_base.py; this file handles
post-processing, output formatting, and figure generation.

Run individual tasks:
    python3 cp1_tasks.py          ← runs all tasks in sequence
    python3 cp1_tasks.py B        ← runs only Task B
───────────────────────────────────────────────────────────────────────────────
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import cp1_base as base


# ─────────────────────────────────────────────────────────────────────────────
# Shared plot style
# ─────────────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 12,
    'figure.dpi': 120,
    'lines.linewidth': 1.8,
})


# ═══════════════════════════════════════════════════════════════════════════════
# TASK B  —  Baseline Case  (t₂ = 10 mm)
# ═══════════════════════════════════════════════════════════════════════════════

def task_B():
    print("\n" + "═" * 60)
    print("TASK B — Baseline Case  (t₂ = 10 mm, Q_int = 50 W)")
    print("═" * 60)

    sol = base.solve(t2=0.010)

    # ── Key temperatures ──
    print(f"\n  Tg          = {sol['Tg']:.3f} K    (limit ≤ 360 K)  "
          f"{'PASS ✓' if sol['Tg'] <= 360 else 'FAIL ✗'}")
    print(f"  Tw = T(r0)  = {sol['Tw']:.3f} K")
    print(f"  Ts = T(r3)  = {sol['Ts']:.3f} K    (limit ≤ 450 K)  "
          f"{'PASS ✓' if sol['Ts'] <= 450 else 'FAIL ✗'}")
    print(f"  Q           = {sol['Q']:.2f} W")

    # ── Interface temperature jumps ──
    # Convention: ΔTc = T(inner side) − T(outer side)  > 0  (drop in heat-flow dir.)
    print(f"\n  ΔTc,12 = T(r1⁻) − T(r1⁺) = {sol['dTc12']:.5f} K")
    print(f"  ΔTc,23 = T(r2⁻) − T(r2⁺) = {sol['dTc23']:.5f} K")

    # ── Additional interface temperatures for reference ──
    print(f"\n  T(r1⁻)  [Al   side] = {sol['T1m']:.3f} K")
    print(f"  T(r1⁺)  [Aero side] = {sol['T1p']:.3f} K")
    print(f"  T(r2⁻)  [Aero side] = {sol['T2m']:.3f} K")
    print(f"  T(r2⁺)  [CFRP side] = {sol['T2p']:.3f} K")

    # ── Pass/fail summary ──
    both_pass = (sol['Tg'] <= 360) and (sol['Ts'] <= 450)
    print(f"\n  Design requirements: {'BOTH SATISFIED ✓' if both_pass else 'NOT SATISFIED ✗'}")

    return sol


# ═══════════════════════════════════════════════════════════════════════════════
# TASK C  —  Maximum Insulation Thickness  (sweep t₂ ∈ [1, 40] mm)
# ═══════════════════════════════════════════════════════════════════════════════

def task_C():
    print("\n" + "═" * 60)
    print("TASK C — Maximum Insulation Thickness t₂")
    print("═" * 60)

    t2_arr = np.linspace(base.T2_MIN, base.T2_MAX, 200)   # [m]
    Tg_arr = np.zeros_like(t2_arr)
    Ts_arr = np.zeros_like(t2_arr)

    for i, t2 in enumerate(t2_arr):
        sol = base.solve(t2)
        Tg_arr[i] = sol['Tg']
        Ts_arr[i] = sol['Ts']

    # ── Find t2_max: largest t2 where Tg ≤ 360 K and Ts ≤ 450 K ──
    # In practice Ts ≪ 450 K across the full range; Tg is the binding constraint.
    feasible = (Tg_arr <= 360.0) & (Ts_arr <= 450.0)

    if feasible.any():
        # largest index in the feasible region
        t2_max = t2_arr[feasible][-1]
        print(f"\n  Maximum allowable t₂ = {t2_max * 1e3:.2f} mm")
        print(f"    (Tg = {Tg_arr[feasible][-1]:.3f} K ≤ 360 K,  "
              f"Ts = {Ts_arr[feasible][-1]:.3f} K ≤ 450 K)")
    else:
        t2_max = None
        print("\n  No feasible t₂ found in [1, 40] mm !")

    # ── Plot 1: Tg(t₂) ──
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.suptitle("Task C — Temperature vs. Aerogel Thickness $t_2$", fontsize=13)

    ax = axes[0]
    ax.plot(t2_arr * 1e3, Tg_arr, color='tab:red', label=r'$T_g(t_2)$')
    ax.axhline(360, color='k', linestyle='--', linewidth=1.2, label='360 K limit')
    if t2_max is not None:
        ax.axvline(t2_max * 1e3, color='tab:orange', linestyle=':', linewidth=1.5,
                   label=f'$t_{{2,max}}$ = {t2_max*1e3:.1f} mm')
    ax.set_xlabel('$t_2$ (mm)')
    ax.set_ylabel('$T_g$ (K)')
    ax.set_title('Cavity gas temperature $T_g$')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # ── Plot 2: Ts(t₂) ──
    ax = axes[1]
    ax.plot(t2_arr * 1e3, Ts_arr, color='tab:blue', label=r'$T_s(t_2)$')
    ax.axhline(450, color='k', linestyle='--', linewidth=1.2, label='450 K limit')
    if t2_max is not None:
        ax.axvline(t2_max * 1e3, color='tab:orange', linestyle=':', linewidth=1.5,
                   label=f'$t_{{2,max}}$ = {t2_max*1e3:.1f} mm')
    ax.set_xlabel('$t_2$ (mm)')
    ax.set_ylabel('$T_s$ (K)')
    ax.set_title('Outer surface temperature $T_s$')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('task_C.png', bbox_inches='tight')
    plt.show()
    print("  → Saved: task_C.png")

    return t2_max


# ═══════════════════════════════════════════════════════════════════════════════
# TASK D  —  Radial Temperature Distribution  (use t₂_max from Task C)
# ═══════════════════════════════════════════════════════════════════════════════

def task_D(t2_max):
    print("\n" + "═" * 60)
    print(f"TASK D — Radial Temperature Distribution  (t₂ = {t2_max*1e3:.2f} mm)")
    print("═" * 60)

    sol = base.solve(t2_max)
    r_arr, T_arr = base.temperature_profile(sol, n_per_layer=200)

    r0, r1, r2, r3 = sol['r0'], sol['r1'], sol['r2'], sol['r3']
    n = 200   # n_per_layer used above

    # ── Report interface values ──
    # End of Layer 1 (r = r1, Al side)  and start of Layer 2 (r = r1, Aerogel side)
    T_r1_minus = T_arr[n - 1]    # last point of Layer 1
    T_r1_plus  = T_arr[n]        # first point of Layer 2
    T_r2_minus = T_arr[2*n - 1]  # last point of Layer 2
    T_r2_plus  = T_arr[2*n]      # first point of Layer 3

    print(f"\n  Interface at r1 = {r1*1e3:.2f} mm:")
    print(f"    T(r1⁻) [Al   side] = {T_r1_minus:.4f} K")
    print(f"    T(r1⁺) [Aero side] = {T_r1_plus:.4f} K")
    print(f"    Jump ΔTc,12        = {T_r1_minus - T_r1_plus:.5f} K")

    print(f"\n  Interface at r2 = {r2*1e3:.2f} mm:")
    print(f"    T(r2⁻) [Aero side] = {T_r2_minus:.4f} K")
    print(f"    T(r2⁺) [CFRP side] = {T_r2_plus:.4f} K")
    print(f"    Jump ΔTc,23        = {T_r2_minus - T_r2_plus:.5f} K")

    # ── Plot T(r) ──
    fig, ax = plt.subplots(figsize=(9, 5))

    r_mm = r_arr * 1e3  # convert to mm for readability
    r_arr_L1 = r_arr[:n];      T_arr_L1 = T_arr[:n]
    r_arr_L2 = r_arr[n:2*n];   T_arr_L2 = T_arr[n:2*n]
    r_arr_L3 = r_arr[2*n:];    T_arr_L3 = T_arr[2*n:]

    # Plot each layer with a distinct color
    ax.plot(r_arr_L1 * 1e3, T_arr_L1, color='tab:blue',   label='Layer 1: Al liner')
    ax.plot(r_arr_L2 * 1e3, T_arr_L2, color='tab:orange', label='Layer 2: Aerogel')
    ax.plot(r_arr_L3 * 1e3, T_arr_L3, color='tab:green',  label='Layer 3: CFRP')

    # Mark contact-resistance jumps with vertical dashed lines
    ax.axvline(r1 * 1e3, color='gray', linestyle='--', linewidth=0.9, alpha=0.7)
    ax.axvline(r2 * 1e3, color='gray', linestyle='--', linewidth=0.9, alpha=0.7)

    # Annotate the jumps
    ax.annotate(f'ΔTc,12 = {sol["dTc12"]:.4f} K',
                xy=(r1 * 1e3, (T_r1_minus + T_r1_plus) / 2),
                xytext=(r1 * 1e3 - 1.5, (T_r1_minus + T_r1_plus) / 2 + 0.5),
                fontsize=9, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    ax.annotate(f'ΔTc,23 = {sol["dTc23"]:.4f} K',
                xy=(r2 * 1e3, (T_r2_minus + T_r2_plus) / 2),
                xytext=(r2 * 1e3 + 0.3, (T_r2_minus + T_r2_plus) / 2 + 0.5),
                fontsize=9, color='gray',
                arrowprops=dict(arrowstyle='->', color='gray', lw=0.8))

    # Shade layers for visual clarity
    ax.axvspan(r0 * 1e3, r1 * 1e3, alpha=0.07, color='tab:blue')
    ax.axvspan(r1 * 1e3, r2 * 1e3, alpha=0.07, color='tab:orange')
    ax.axvspan(r2 * 1e3, r3 * 1e3, alpha=0.07, color='tab:green')

    ax.set_xlabel('Radius $r$ (mm)')
    ax.set_ylabel('Temperature $T$ (K)')
    ax.set_title(f'Task D — Radial Temperature Distribution  '
                 f'($t_2$ = {t2_max*1e3:.1f} mm, $Q_{{int}}$ = {sol["Q_int"]:.0f} W)')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('task_D.png', bbox_inches='tight')
    plt.show()
    print("  → Saved: task_D.png")

    return sol


# ═══════════════════════════════════════════════════════════════════════════════
# TASK E  —  Electronics Power Sensitivity  (fix t₂ = t₂_max, sweep Q_int)
# ═══════════════════════════════════════════════════════════════════════════════

def task_E(t2_max):
    print("\n" + "═" * 60)
    print(f"TASK E — Electronics Power Sensitivity  (t₂ = {t2_max*1e3:.2f} mm)")
    print("═" * 60)

    Q_arr  = np.linspace(20.0, 200.0, 200)   # [W]
    Tg_arr = np.zeros_like(Q_arr)
    Ts_arr = np.zeros_like(Q_arr)

    for i, Q_int in enumerate(Q_arr):
        sol = base.solve(t2_max, Q_int=Q_int)
        Tg_arr[i] = sol['Tg']
        Ts_arr[i] = sol['Ts']

    # ── Find Q_int,max such that Tg = 360 K ──
    # Tg is monotonically increasing with Q_int; find crossing by interpolation.
    idx = np.searchsorted(Tg_arr, 360.0)

    if 0 < idx < len(Q_arr):
        # Linear interpolation between idx-1 and idx
        Q_lo, Q_hi = Q_arr[idx - 1], Q_arr[idx]
        Tg_lo, Tg_hi = Tg_arr[idx - 1], Tg_arr[idx]
        Q_int_max = Q_lo + (360.0 - Tg_lo) / (Tg_hi - Tg_lo) * (Q_hi - Q_lo)
        print(f"\n  Q_int,max (Tg = 360 K) = {Q_int_max:.2f} W")
    else:
        Q_int_max = None
        print("\n  Tg does not reach 360 K in the sweep range [20, 200] W")

    # ── Plot: Tg vs Q_int  and  Ts vs Q_int ──
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.suptitle(f"Task E — Sensitivity to $Q_{{int}}$  "
                 f"($t_2$ = {t2_max*1e3:.1f} mm)", fontsize=13)

    ax = axes[0]
    ax.plot(Q_arr, Tg_arr, color='tab:red', label=r'$T_g(Q_{int})$')
    ax.axhline(360, color='k', linestyle='--', linewidth=1.2, label='360 K limit')
    if Q_int_max is not None:
        ax.axvline(Q_int_max, color='tab:orange', linestyle=':', linewidth=1.5,
                   label=f'$Q_{{int,max}}$ = {Q_int_max:.1f} W')
    ax.set_xlabel('$Q_{int}$ (W)')
    ax.set_ylabel('$T_g$ (K)')
    ax.set_title('Cavity gas temperature $T_g$')
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(Q_arr, Ts_arr, color='tab:blue', label=r'$T_s(Q_{int})$')
    ax.axhline(450, color='k', linestyle='--', linewidth=1.2, label='450 K limit')
    ax.set_xlabel('$Q_{int}$ (W)')
    ax.set_ylabel('$T_s$ (K)')
    ax.set_title('Outer surface temperature $T_s$')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('task_E.png', bbox_inches='tight')
    plt.show()
    print("  → Saved: task_E.png")

    return Q_int_max


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN — run all (or selected) tasks
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':

    # Which tasks to run (default: all)
    run = sys.argv[1].upper() if len(sys.argv) > 1 else 'BCDE'

    if 'B' in run:
        task_B()

    t2_max = None

    if 'C' in run:
        t2_max = task_C()

    if 'D' in run:
        if t2_max is None:
            # Compute t2_max if Task C was skipped
            print("\n  [Task D] Task C not run; computing t2_max now...")
            t2_arr = np.linspace(base.T2_MIN, base.T2_MAX, 200)
            for t2 in t2_arr:
                sol = base.solve(t2)
                if sol['Tg'] > 360.0:
                    break
                t2_max = t2
        task_D(t2_max)

    if 'E' in run:
        if t2_max is None:
            print("\n  [Task E] Task C not run; computing t2_max now...")
            t2_arr = np.linspace(base.T2_MIN, base.T2_MAX, 200)
            for t2 in t2_arr:
                sol = base.solve(t2)
                if sol['Tg'] > 360.0:
                    break
                t2_max = t2
        task_E(t2_max)

    print("\nDone.\n")
