# Heat Transfer Computing HW 2 - parameters and system assembly for parts 2 and 3.

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

# Physical parameters
Nx, Ny = 226, 31
dx = 0.002

k_pv  = 1.0
k_al  = 205.0
h     = 500.0
q_rad = 1500.0

q_term = q_rad * dx / k_pv
Bi     = h * dx / k_al

# Geometry (row indices)
J_INTERFACE = 10          # PV/Al interface  (y = 0.02 m)
J_BOT       = Ny - 1      # bottom row       (y = 0.06 m)

CH_J_TOP = 15             # channel top wall    (y = 0.03 m)
CH_J_BOT = 25             # channel bottom wall (y = 0.05 m)

CH_X_STARTS = [0.10, 0.16, 0.22, 0.28, 0.34]
CH_WIDTH    = 0.01


def build_node_map():
    nm = np.zeros((Ny, Nx), dtype=int)

    for j in range(Ny):
        for i in range(Nx):
            # Top PV surface - radiation flux
            if j == 0:
                if   i == 0:      nm[j, i] = 1
                elif i == Nx - 1: nm[j, i] = 3
                else:             nm[j, i] = 2

            # PV interior rows
            elif 0 < j < J_INTERFACE:
                if   i == 0:      nm[j, i] = 4
                elif i == Nx - 1: nm[j, i] = 6
                else:             nm[j, i] = 5

            # PV / Al interface
            elif j == J_INTERFACE:
                if   i == 0:      nm[j, i] = 7
                elif i == Nx - 1: nm[j, i] = 9
                else:             nm[j, i] = 8

            # Al interior rows
            elif J_INTERFACE < j < J_BOT:
                if   i == 0:      nm[j, i] = 10
                elif i == Nx - 1: nm[j, i] = 12
                else:             nm[j, i] = 11

            # Bottom row - insulated
            elif j == J_BOT:
                if   i == 0:      nm[j, i] = 14
                elif i == Nx - 1: nm[j, i] = 15
                else:             nm[j, i] = 13

    # Overwrite channel nodes
    for x_start in CH_X_STARTS:
        i_ws = int(round(x_start / dx))
        i_we = int(round((x_start + CH_WIDTH) / dx))
        j_ws, j_we = CH_J_TOP, CH_J_BOT

        nm[j_ws + 1:j_we, i_ws + 1:i_we] = 99   # fluid interior

        nm[j_ws,          i_ws + 1:i_we] = 16   # top  wall
        nm[j_we,          i_ws + 1:i_we] = 17   # bot  wall
        nm[j_ws + 1:j_we, i_ws         ] = 18   # left wall
        nm[j_ws + 1:j_we, i_we         ] = 19   # right wall

        nm[j_ws, i_ws] = 20   # TL solid corner
        nm[j_ws, i_we] = 21   # TR solid corner
        nm[j_we, i_ws] = 22   # BL solid corner
        nm[j_we, i_we] = 23   # BR solid corner

    return nm


def build_system(Tw):
    NodeMap = build_node_map()
    N = Nx * Ny
    A = lil_matrix((N, N))
    b = np.zeros(N)

    def idx(i, j):
        return j * Nx + i

    for j in range(Ny):
        for i in range(Nx):
            p = idx(i, j)
            t = NodeMap[j, i]

            # Top PV surface (radiation flux)
            if t == 1:     # top-left corner
                A[p, p] = -2
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j + 1)] = 1
                b[p] = -q_term

            elif t == 2:   # top plane
                A[p, p] = -4
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j + 1)] = 2
                b[p] = -2 * q_term

            elif t == 3:   # top-right corner
                A[p, p] = -2
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i, j + 1)] = 1
                b[p] = -q_term

            # PV insulated sides / interior
            elif t == 4:
                A[p, p] = -4
                A[p, idx(i + 1, j)] = 2
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1

            elif t == 5:
                A[p, p] = -4
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1

            elif t == 6:
                A[p, p] = -4
                A[p, idx(i - 1, j)] = 2
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1

            # PV / Al interface
            elif t == 7:
                A[p, p] = -2 * (k_pv + k_al)
                A[p, idx(i + 1, j)] = k_pv + k_al
                A[p, idx(i, j - 1)] = k_pv
                A[p, idx(i, j + 1)] = k_al

            elif t == 8:
                A[p, p] = -2 * (k_pv + k_al)
                A[p, idx(i - 1, j)] = (k_pv + k_al) / 2
                A[p, idx(i + 1, j)] = (k_pv + k_al) / 2
                A[p, idx(i, j - 1)] = k_pv
                A[p, idx(i, j + 1)] = k_al

            elif t == 9:
                A[p, p] = -2 * (k_pv + k_al)
                A[p, idx(i - 1, j)] = k_pv + k_al
                A[p, idx(i, j - 1)] = k_pv
                A[p, idx(i, j + 1)] = k_al

            # Al insulated sides / interior
            elif t == 10:
                A[p, p] = -4
                A[p, idx(i + 1, j)] = 2
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1

            elif t == 11:
                A[p, p] = -4
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1

            elif t == 12:
                A[p, p] = -4
                A[p, idx(i - 1, j)] = 2
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1

            # Bottom insulated
            elif t == 13:
                A[p, p] = -4
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j - 1)] = 2

            elif t == 14:
                A[p, p] = -2
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j - 1)] = 1

            elif t == 15:
                A[p, p] = -2
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i, j - 1)] = 1

            # Channel flat walls - convection
            elif t == 16:  # top wall, fluid below
                A[p, p] = -2 * (Bi + 2)
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j - 1)] = 2
                b[p] = -2 * Bi * Tw

            elif t == 17:  # bottom wall, fluid above
                A[p, p] = -2 * (Bi + 2)
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j + 1)] = 2
                b[p] = -2 * Bi * Tw

            elif t == 18:  # left wall, fluid to the right
                A[p, p] = -2 * (Bi + 2)
                A[p, idx(i - 1, j)] = 2
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1
                b[p] = -2 * Bi * Tw

            elif t == 19:  # right wall, fluid to the left
                A[p, p] = -2 * (Bi + 2)
                A[p, idx(i + 1, j)] = 2
                A[p, idx(i, j - 1)] = 1
                A[p, idx(i, j + 1)] = 1
                b[p] = -2 * Bi * Tw

            # Channel solid internal corners - convection
            elif t == 20:  # TL, fluid at bottom-right
                A[p, p] = -2 * (Bi + 3)
                A[p, idx(i - 1, j)] = 2
                A[p, idx(i, j - 1)] = 2
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j + 1)] = 1
                b[p] = -2 * Bi * Tw

            elif t == 21:  # TR, fluid at bottom-left
                A[p, p] = -2 * (Bi + 3)
                A[p, idx(i + 1, j)] = 2
                A[p, idx(i, j - 1)] = 2
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i, j + 1)] = 1
                b[p] = -2 * Bi * Tw

            elif t == 22:  # BL, fluid at top-right
                A[p, p] = -2 * (Bi + 3)
                A[p, idx(i - 1, j)] = 2
                A[p, idx(i, j + 1)] = 2
                A[p, idx(i + 1, j)] = 1
                A[p, idx(i, j - 1)] = 1
                b[p] = -2 * Bi * Tw

            elif t == 23:  # BR, fluid at top-left
                A[p, p] = -2 * (Bi + 3)
                A[p, idx(i + 1, j)] = 2
                A[p, idx(i, j + 1)] = 2
                A[p, idx(i - 1, j)] = 1
                A[p, idx(i, j - 1)] = 1
                b[p] = -2 * Bi * Tw

            # Fluid dummy
            elif t == 99:
                A[p, p] = 1
                b[p] = Tw

    return A.tocsr(), b


def solve_system(Tw):
    A, b = build_system(Tw)
    return spsolve(A, b).reshape((Ny, Nx))
