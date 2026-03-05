"""
cp1_tasks.py — Tasks B-E for CP1.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import cp1_base as base


def task_B():
    print("\nTASK B -- Baseline (t2=10mm, Q_int=50W)")
    sol = base.solve(t2=0.010)

    print(f"Tg = {sol['Tg']:.3f} K, {'PASS' if sol['Tg'] <= 360 else 'FAIL'}")
    print(f"Tw = {sol['Tw']:.3f} K")
    print(f"Ts = {sol['Ts']:.3f} K, {'PASS' if sol['Ts'] <= 450 else 'FAIL'}")
    print(f"Q = {sol['Q']:.2f} W")
    print(f"dTc12 = {sol['dTc12']:.5f} K")
    print(f"dTc23 = {sol['dTc23']:.5f} K")

    both_pass = (sol['Tg'] <= 360) and (sol['Ts'] <= 450)
    print(f"Design: {'BOTH PASS' if both_pass else 'NOT SATISFIED'}")
    return sol


def task_C():
    print("\nTASK C -- Max insulation thickness")
    t2_arr = np.linspace(base.T2_MIN, base.T2_MAX, 200)
    Tg_arr = np.zeros_like(t2_arr)
    Ts_arr = np.zeros_like(t2_arr)

    for i, t2 in enumerate(t2_arr):
        sol = base.solve(t2)
        Tg_arr[i] = sol['Tg']
        Ts_arr[i] = sol['Ts']

    feasible = (Tg_arr <= 360.0) & (Ts_arr <= 450.0)
    t2_max = t2_arr[feasible][-1] if feasible.any() else None

    if t2_max is not None:
        print(f"t2_max = {t2_max*1e3:.2f} mm")
    else:
        print("No feasible t2 found")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t2_arr * 1e3, Tg_arr, 'r-', label='Tg (cavity gas)')
    ax.plot(t2_arr * 1e3, Ts_arr, 'b-', label='Ts (outer surface)')
    ax.axhline(360, color='r', linestyle='--', label='Tg limit (360 K)')
    ax.axhline(450, color='b', linestyle='--', label='Ts limit (450 K)')
    ax.set_xlabel('Aerogel Thickness t2 (mm)')
    ax.set_ylabel('Temperature (K)')
    ax.set_title('Task C: Temperature vs Aerogel Thickness')
    ax.legend()
    plt.tight_layout()
    plt.savefig('task_C.png')
    plt.show()
    return t2_max


def task_D(t2_max):
    print(f"\nTASK D -- T(r) profile (t2={t2_max*1e3:.2f}mm)")
    sol = base.solve(t2_max, N=100)
    r_arr, T_arr = base.temperature_profile(sol)
    N = sol['N']

    print(f"dTc12 = {sol['dTc12']:.5f} K")
    print(f"dTc23 = {sol['dTc23']:.5f} K")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(r_arr[:N] * 1e3, T_arr[:N], 'b-', label='Layer 1: Al liner')
    ax.plot(r_arr[N:2*N] * 1e3, T_arr[N:2*N], 'r-', label='Layer 2: Aerogel')
    ax.plot(r_arr[2*N:] * 1e3, T_arr[2*N:], 'g-', label='Layer 3: CFRP')
    ax.axvline(sol['r1'] * 1e3, color='gray', linestyle='--', linewidth=0.8)
    ax.axvline(sol['r2'] * 1e3, color='gray', linestyle='--', linewidth=0.8)
    ax.set_xlabel('Radius r (mm)')
    ax.set_ylabel('Temperature T (K)')
    ax.set_title(f'Task D: Radial Temperature Distribution (t2 = {t2_max*1e3:.1f} mm)')
    ax.legend()
    plt.tight_layout()
    plt.savefig('task_D.png')
    plt.show()
    return sol


def task_E(t2_max):
    print(f"\nTASK E -- Sensitivity (t2={t2_max*1e3:.2f}mm)")
    Q_arr = np.linspace(20.0, 200.0, 200)
    Tg_arr = np.zeros_like(Q_arr)

    for i, Q_int in enumerate(Q_arr):
        sol = base.solve(t2_max, Q_int=Q_int)
        Tg_arr[i] = sol['Tg']

    idx = np.searchsorted(Tg_arr, 360.0)
    if 0 < idx < len(Q_arr):
        frac = (360.0 - Tg_arr[idx-1]) / (Tg_arr[idx] - Tg_arr[idx-1])
        Q_int_max = Q_arr[idx-1] + frac * (Q_arr[idx] - Q_arr[idx-1])
        print(f"Q_int_max (Tg=360K) = {Q_int_max:.2f} W")
    else:
        Q_int_max = None
        print("Tg does not reach 360K in range")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(Q_arr, Tg_arr, 'r-', label='Tg (cavity gas)')
    ax.axhline(360, color='k', linestyle='--', label='Tg limit (360 K)')
    ax.set_xlabel('Electronics Power Q_int (W)')
    ax.set_ylabel('Cavity Gas Temperature Tg (K)')
    ax.set_title(f'Task E: Tg Sensitivity to Q_int (t2 = {t2_max*1e3:.1f} mm)')
    ax.legend()
    plt.tight_layout()
    plt.savefig('task_E.png')
    plt.show()
    return Q_int_max


def _get_t2_max():
    t2_arr = np.linspace(base.T2_MIN, base.T2_MAX, 200)
    t2_max = t2_arr[0]
    for t2 in t2_arr:
        sol = base.solve(t2)
        if sol['Tg'] > 360.0:
            break
        t2_max = t2
    return t2_max


if __name__ == '__main__':
    run = sys.argv[1].upper() if len(sys.argv) > 1 else 'BCDE'

    if 'B' in run:
        task_B()

    t2_max = None
    if 'C' in run:
        t2_max = task_C()

    if 'D' in run:
        if t2_max is None:
            t2_max = _get_t2_max()
        task_D(t2_max)

    if 'E' in run:
        if t2_max is None:
            t2_max = _get_t2_max()
        task_E(t2_max)

    print("\nDone.")
