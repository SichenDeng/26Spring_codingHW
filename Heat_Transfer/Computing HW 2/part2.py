"""
part2.py — Solve the 2D steady-state heat conduction problem for T_w = 20 °C
and visualise the temperature distribution.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.sparse.linalg import spsolve

from params import (
    NX, NY, D, J_INT, J_CH_TOP, J_CH_BOT, CH_IL, CH_IR,
    x_nodes, y_nodes, N_NODES,
    build_system, reconstruct_field,
)

# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
T_w = 20.0   # °C — coolant temperature (Part 2)

A, b = build_system(T_w)
T_vec = spsolve(A, b)
T_field = reconstruct_field(T_vec)

# Location of maximum temperature
T_max = np.nanmax(T_field)
jmax, imax = np.unravel_index(np.nanargmax(T_field), T_field.shape)
x_max = x_nodes[imax]
y_max = y_nodes[jmax]

print(f"Maximum temperature : {T_max:.4f} °C")
print(f"Location            : x = {x_max:.4f} m,  y = {y_max:.4f} m")

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(12, 4))

X, Y = np.meshgrid(x_nodes, y_nodes)

# Fill NaN (channel voids) with interpolated values for smooth contour plotting
T_plot = T_field.copy()
from scipy.ndimage import generic_filter
mask = np.isnan(T_plot)
if mask.any():
    # Fill void NaNs with nearest valid neighbor average for contour rendering
    from scipy.interpolate import NearestNDInterpolator
    valid = ~mask
    coords = np.array(np.where(valid)).T
    values = T_plot[valid]
    interp = NearestNDInterpolator(coords, values)
    void_coords = np.array(np.where(mask)).T
    T_plot[mask] = interp(void_coords)

# Filled contour plot with contour lines to show curvature
levels = np.linspace(np.nanmin(T_field), np.nanmax(T_field), 30)
cf = ax.contourf(X, Y, T_plot, levels=levels, cmap='hot')
# Add contour lines to make the curvature visible
cs = ax.contour(X, Y, T_plot, levels=levels, colors='k', linewidths=0.3, alpha=0.5)
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label('Temperature (°C)')

# Invert y-axis so j=0 (top surface) appears at the top of the plot
ax.invert_yaxis()

# PV / Al interface — cyan dashed horizontal line at y = J_INT * D
y_int = J_INT * D
ax.axhline(y=y_int, color='cyan', linestyle='--', linewidth=1.2,
           label=f'PV/Al interface  (y = {y_int*100:.0f} cm)')

# Water channel regions — blue filled rectangles
y_ch_top = J_CH_TOP * D
y_ch_bot = J_CH_BOT * D
ch_height = y_ch_bot - y_ch_top
first_ch  = True
for iL, iR in zip(CH_IL, CH_IR):
    x_ch_left  = iL * D
    ch_width   = (iR - iL) * D
    label = 'Water channels' if first_ch else None
    rect = patches.Rectangle(
        (x_ch_left, y_ch_top), ch_width, ch_height,
        linewidth=1, edgecolor='blue', facecolor='blue', alpha=0.5,
        label=label,
    )
    ax.add_patch(rect)
    first_ch = False

# Maximum temperature marker
ax.plot(x_max, y_max, marker='*', markersize=12, color='lime',
        markeredgecolor='black', markeredgewidth=0.5,
        label=f'T_max = {T_max:.2f} °C')

ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_title('2D Temperature Distribution — Solar Panel Cross-Section (T_w = 20 °C)')
ax.legend(loc='upper right', fontsize=8)

fig.tight_layout()
fig.savefig(
    '/Users/sichendeng/JHU/26Spring/26Spring_codingHW/Heat_Transfer/Computing HW 2/temperature_distribution.png',
    dpi=150,
)
plt.show()
print("Figure saved: temperature_distribution.png")
