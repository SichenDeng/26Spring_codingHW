import numpy as np
import matplotlib.pyplot as plt

from params import (Nx, Ny, dx, build_node_map, solve_system,
                    CH_X_STARTS, CH_WIDTH, CH_J_TOP, CH_J_BOT)

T_PANEL_LIMIT = 35.0

T_w = 20.0
T = solve_system(T_w)

NodeMap = build_node_map()
T_solid = np.where(NodeMap == 99, np.nan, T)

# Hot-spot
T_max = np.nanmax(T_solid)
T_min = np.nanmin(T_solid)
jmax, imax = np.unravel_index(np.nanargmax(T_solid), T_solid.shape)
x_max, y_max = imax * dx, jmax * dx

# Tw,max from the linear offset
Tw_max = T_PANEL_LIMIT - (T_max - T_w)

# Shared colorbar range covering both Part 2 and Part 3 fields
vmin = min(T_min, T_min + (Tw_max - T_w))
vmax = max(T_max, T_max + (Tw_max - T_w))

print("--- Problem 2 ---")
print(f"Tw    = {T_w:.2f} C")
print(f"Tmax  = {T_max:.4f} C   at  (x, y) = ({x_max:.4f}, {y_max:.4f}) m")

# Plot
X = np.linspace(0.0, (Nx - 1) * dx, Nx)
Y = np.linspace(0.0, (Ny - 1) * dx, Ny)
Y_SCALE = 4.0
levels = np.linspace(vmin, vmax, 80)

fig, ax = plt.subplots(figsize=(13, 5.0))
cf = ax.contourf(X, Y, T_solid, levels=levels, cmap='jet',
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
ax.set_title(f'Temperature distribution  ($T_w$ = {T_w:.0f} °C,  y stretched {Y_SCALE:g}x)')

cbar = fig.colorbar(cf, ax=ax, pad=0.015, fraction=0.03)
cbar.set_label('T (°C)')

ax.plot(x_max, y_max, 'kx', markersize=10, markeredgewidth=2)
ax.annotate(f'$T_{{max}}$ = {T_max:.2f} °C',
            xy=(x_max, y_max),
            xytext=(x_max - 0.09, y_max + 0.008),
            fontsize=10,
            arrowprops=dict(arrowstyle='->', color='black', lw=1))

fig.tight_layout()
plt.show()
