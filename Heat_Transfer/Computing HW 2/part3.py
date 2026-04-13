import numpy as np
import matplotlib.pyplot as plt

from params import (Nx, Ny, dx, build_node_map, solve_system,
                    CH_X_STARTS, CH_WIDTH, CH_J_TOP, CH_J_BOT)

T_PANEL_LIMIT = 35.0

NodeMap = build_node_map()
solid_mask = NodeMap != 99

# Reference solve
T_w_ref = 20.0
T_ref = solve_system(T_w_ref)
T_ref_solid = np.where(solid_mask, T_ref, np.nan)
T_max_ref = float(np.nanmax(T_ref_solid))
T_min_ref = float(np.nanmin(T_ref_solid))

# Linear offset
C = T_max_ref - T_w_ref
T_w_max = T_PANEL_LIMIT - C

# Shifted field at Tw = Tw_max
T_new = T_ref_solid + (T_w_max - T_w_ref)
T_min_new = T_min_ref + (T_w_max - T_w_ref)
T_max_new = T_PANEL_LIMIT

print("--- Problem 3 ---")
print(f"At Tw = {T_w_ref:.1f} C  : Tmax = {T_max_ref:.4f} C")
print(f"Offset C = Tmax - Tw     : {C:.4f} C")
print(f"Max allowable Tw         : {T_w_max:.4f} C")

# Hot-spot location (same indices as Part 2)
jmax, imax = np.unravel_index(np.nanargmax(T_new), T_new.shape)
x_max, y_max = imax * dx, jmax * dx

# Shared colorbar range with Part 2
vmin = min(T_min_ref, T_min_new)
vmax = max(T_max_ref, T_max_new)

# Plot
X = np.linspace(0.0, (Nx - 1) * dx, Nx)
Y = np.linspace(0.0, (Ny - 1) * dx, Ny)
Y_SCALE = 4.0
levels = np.linspace(vmin, vmax, 80)

fig, ax = plt.subplots(figsize=(13, 5.0))
cf = ax.contourf(X, Y, T_new, levels=levels, cmap='jet',
                 vmin=vmin, vmax=vmax, extend='both')

for x0 in CH_X_STARTS:
    ax.add_patch(plt.Rectangle((x0, CH_J_TOP * dx),
                               CH_WIDTH, (CH_J_BOT - CH_J_TOP) * dx,
                               fill=False, edgecolor='white',
                               linewidth=0.8, linestyle='--'))
ax.axhline(10 * dx, color='white', linewidth=0.6, linestyle=':')

ax.set_aspect(Y_SCALE)
ax.invert_yaxis()
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title(f'Temperature distribution  ($T_w$ = {T_w_max:.2f} °C,  y stretched {Y_SCALE:g}x)')

cbar = fig.colorbar(cf, ax=ax, pad=0.015, fraction=0.03)
cbar.set_label('T (°C)')

ax.plot(x_max, y_max, 'kx', markersize=10, markeredgewidth=2)
ax.annotate(f'$T_{{max}}$ = {T_PANEL_LIMIT:.2f} °C',
            xy=(x_max, y_max),
            xytext=(x_max - 0.09, y_max + 0.008),
            fontsize=10,
            arrowprops=dict(arrowstyle='->', color='black', lw=1))

fig.tight_layout()
plt.show()
