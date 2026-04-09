"""
part3.py — Determine the maximum allowable coolant temperature T_w such that
the panel temperature never exceeds 35 °C.

Linearity argument
------------------
Let T' = T - T_w.  Substituting into the governing equations:

  * Interior nodes: the Laplacian is linear, so ∇²T' = 0 (same as ∇²T = 0).
  * Top surface BC: -k_P ∂T'/∂y = q_pp  (unchanged, T_w cancels).
  * Insulated faces: ∂T'/∂n = 0  (unchanged).
  * Channel wall BC: -k_A ∂T'/∂n = h(T - T_w) = h·T'  → homogeneous Robin BC.
  * PV/Al interface: continuity of T' and normal flux (unchanged in form).

Because every boundary condition for T' is independent of T_w, the shifted
field T' is the same for all values of T_w.  In particular, T_max' = T_max - T_w
equals a constant C that depends only on the material properties and heat flux.

Therefore:  T_max(T_w) = C + T_w   (strictly linear in T_w).

Setting T_max = 35 °C:
    35 = C + T_w_max
    T_w_max = 35 - C = 35 - (T_max_20 - 20) = 55 - T_max_20

where T_max_20 is the maximum temperature from Part 2 (T_w = 20 °C).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve

from params import build_system, reconstruct_field

T_LIMIT = 35.0    # °C — maximum allowable panel temperature

# ---------------------------------------------------------------------------
# Compute T_max for T_w = 20 °C  (reproduces Part 2 result)
# ---------------------------------------------------------------------------
T_w_ref = 20.0
A, b = build_system(T_w_ref)
T_vec = spsolve(A, b)
T_field = reconstruct_field(T_vec)
T_max_20 = float(np.nanmax(T_field))

# Constant offset: T_max - T_w is material/flux-dependent, T_w-independent
C = T_max_20 - T_w_ref

T_w_max = T_LIMIT - C   # = 150 - T_max_20

print(f"T_max at T_w = {T_w_ref:.1f} °C  : {T_max_20:.4f} °C")
print(f"Linear offset C = T_max - T_w : {C:.4f} °C")
print(f"Maximum allowable T_w         : {T_w_max:.4f} °C")

# ---------------------------------------------------------------------------
# Numerical verification: solve for six values of T_w and check linearity
# ---------------------------------------------------------------------------
T_w_vals  = np.linspace(0.0, 50.0, 6)
T_max_vals = []

for Tw in T_w_vals:
    A_i, b_i = build_system(Tw)
    T_v = spsolve(A_i, b_i)
    T_max_vals.append(float(np.nanmax(reconstruct_field(T_v))))

T_max_vals = np.array(T_max_vals)

# ---------------------------------------------------------------------------
# Plot T_max vs T_w
# ---------------------------------------------------------------------------
T_w_line   = np.array([T_w_vals[0], T_w_vals[-1]])
T_max_line = C + T_w_line   # predicted linear trend

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(T_w_vals, T_max_vals, 'o', color='steelblue',
        markersize=8, label='Numerical solutions')
ax.plot(T_w_line, T_max_line, '--', color='steelblue', linewidth=1.2,
        label='Linear fit: $T_{max} = C + T_w$')
ax.axhline(T_LIMIT, color='red', linewidth=1.5,
           label=f'Panel limit = {T_LIMIT:.0f} °C')
ax.axvline(T_w_max, color='green', linestyle=':', linewidth=1.5,
           label=f'$T_{{w,max}}$ = {T_w_max:.2f} °C')

ax.set_xlabel('Coolant temperature $T_w$ (°C)')
ax.set_ylabel('Maximum panel temperature $T_{max}$ (°C)')
ax.set_title('$T_{max}$ vs. $T_w$ — Linearity verification')
ax.legend(fontsize=9)
ax.grid(True, linestyle=':', alpha=0.5)

fig.tight_layout()
fig.savefig(
    '/Users/sichendeng/JHU/26Spring/26Spring_codingHW/Heat_Transfer/Computing HW 2/Tmax_vs_Tw.png',
    dpi=150,
)
plt.show()
print("Figure saved: Tmax_vs_Tw.png")
