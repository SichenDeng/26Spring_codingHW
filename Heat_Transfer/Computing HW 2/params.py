"""
params.py — Shared parameters and system assembly for 2D steady-state heat
conduction in a solar panel cross-section (finite difference, energy balance).
"""

import numpy as np
from scipy.sparse import lil_matrix

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
k_P   = 1.0     # W/(m·K)  — PV panel thermal conductivity
k_A   = 205.0   # W/(m·K)  — aluminum thermal conductivity
h     = 500.0   # W/(m²·K) — convection coefficient inside water channels
q_pp  = 1500.0  # W/m²     — net radiation heat flux absorbed at top surface

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
L_x  = 0.45    # m — domain width
L_y  = 0.06    # m — total domain height  (H_PV=0.02 m + H_Al=0.04 m)
H_PV = 0.02    # m — PV layer thickness
H_Al = 0.04    # m — aluminum layer thickness

# ---------------------------------------------------------------------------
# Grid  (uniform spacing Δ = 0.002 m)
# ---------------------------------------------------------------------------
NX  = 226          # nodes in x-direction (i = 0 … 225)
NY  = 31           # nodes in y-direction (j = 0 … 30)
D   = 0.002        # m — Δx = Δy = Δ

# Derived convenience
x_nodes = np.arange(NX) * D   # x-coordinates of nodes
y_nodes = np.arange(NY) * D   # y-coordinates (j=0 is top surface)

# Layer boundary row index
J_INT = 10   # j index of PV/Al interface  (y = 0.02 m)

# ---------------------------------------------------------------------------
# Channel geometry — five rectangular water channels in the aluminum layer
#   j_top = 15  (y = 0.030 m),  j_bot = 25  (y = 0.050 m)
#   void nodes: strictly interior, i.e. iL < i < iR  AND  j_top < j < j_bot
# ---------------------------------------------------------------------------
J_CH_TOP = 15
J_CH_BOT = 25

# Left and right wall column indices for each of the 5 channels
CH_IL = [50,  80,  110, 140, 170]
CH_IR = [55,  85,  115, 145, 175]

# ---------------------------------------------------------------------------
# Helper: classify a node as void (inside a channel cavity)
# ---------------------------------------------------------------------------
def is_void(i, j):
    """Return True if node (i,j) is strictly inside a channel void (no material)."""
    if J_CH_TOP < j < J_CH_BOT:
        for iL, iR in zip(CH_IL, CH_IR):
            if iL < i < iR:
                return True
    return False


# ---------------------------------------------------------------------------
# Helper: classify channel-wall nodes
#   Returns a string tag or None if the node is not on a channel wall.
#   Priority: corner > top/bottom wall > left/right wall
# ---------------------------------------------------------------------------
def get_wall_type(i, j):
    """
    Classify node (i,j) with respect to channel walls.
    Returns one of:
        'corner_TL', 'corner_TR', 'corner_BL', 'corner_BR'
        'wall_top', 'wall_bot', 'wall_left', 'wall_right'
        None  — not a channel wall node
    """
    for iL, iR in zip(CH_IL, CH_IR):
        # corner nodes (Δ/2 × Δ/2 control volume, two convection faces)
        if j == J_CH_TOP and i == iL:
            return 'corner_TL'
        if j == J_CH_TOP and i == iR:
            return 'corner_TR'
        if j == J_CH_BOT and i == iL:
            return 'corner_BL'
        if j == J_CH_BOT and i == iR:
            return 'corner_BR'
        # flat-wall nodes (Δ/2 × Δ control volume, one convection face)
        if j == J_CH_TOP and iL < i < iR:
            return 'wall_top'
        if j == J_CH_BOT and iL < i < iR:
            return 'wall_bot'
        if i == iL and J_CH_TOP < j < J_CH_BOT:
            return 'wall_left'
        if i == iR and J_CH_TOP < j < J_CH_BOT:
            return 'wall_right'
    return None


# ---------------------------------------------------------------------------
# Build node map: assign a row index in the linear system to every live node.
# Void nodes get index -1.
# ---------------------------------------------------------------------------
def build_node_map():
    """
    Returns
    -------
    NODE_MAP  : (NY, NX) int array  — row index for each (j, i), -1 if void
    NODE_LIST : list of (i, j) tuples in row-index order
    """
    nm   = np.full((NY, NX), -1, dtype=int)
    lst  = []
    idx  = 0
    for j in range(NY):
        for i in range(NX):
            if not is_void(i, j):
                nm[j, i] = idx
                lst.append((i, j))
                idx += 1
    return nm, lst


NODE_MAP, NODE_LIST = build_node_map()
N_NODES = len(NODE_LIST)


# ---------------------------------------------------------------------------
# Assemble the sparse linear system  A · T_vec = b
# ---------------------------------------------------------------------------
def build_system(T_w):
    """
    Assemble the N_NODES × N_NODES sparse system for the given coolant
    temperature T_w (°C).

    Returns
    -------
    A : scipy CSR matrix
    b : 1-D numpy array
    """
    A  = lil_matrix((N_NODES, N_NODES))
    b  = np.zeros(N_NODES)
    nm = NODE_MAP

    hDkA  = h * D / k_A          # dimensionless convection parameter
    k_avg = 0.5 * (k_P + k_A)    # arithmetic mean conductivity at interface

    def row(i, j):
        return nm[j, i]

    def _fill_wall(A, b, r, i, j, wtype, hDkA, T_w):
        """
        Energy-balance FD equations for the eight channel-wall node types.
        Each derives from integrating the heat equation over the appropriate
        half or quarter control volume with one or two convective faces.
        """
        if wtype == 'wall_top':
            # Top wall: Δ×Δ/2 CV, convection on bottom face into channel
            # 2·T[i,j-1] + T[i-1,j] + T[i+1,j] − (4 + 2hΔ/k_A)·T[i,j] = −2(hΔ/k_A)·T_w
            A[r, row(i,   j-1)] += 2.0
            A[r, row(i-1, j  )] += 1.0
            A[r, row(i+1, j  )] += 1.0
            A[r, r             ] += -(4.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

        elif wtype == 'wall_bot':
            # Bottom wall: Δ×Δ/2 CV, convection on top face into channel
            # 2·T[i,j+1] + T[i-1,j] + T[i+1,j] − (4 + 2hΔ/k_A)·T[i,j] = −2(hΔ/k_A)·T_w
            A[r, row(i,   j+1)] += 2.0
            A[r, row(i-1, j  )] += 1.0
            A[r, row(i+1, j  )] += 1.0
            A[r, r             ] += -(4.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

        elif wtype == 'wall_left':
            # Left wall: Δ/2×Δ CV, convection on right face into channel
            # 2·T[i-1,j] + T[i,j-1] + T[i,j+1] − (4 + 2hΔ/k_A)·T[i,j] = −2(hΔ/k_A)·T_w
            A[r, row(i-1, j  )] += 2.0
            A[r, row(i,   j-1)] += 1.0
            A[r, row(i,   j+1)] += 1.0
            A[r, r             ] += -(4.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

        elif wtype == 'wall_right':
            # Right wall: Δ/2×Δ CV, convection on left face into channel
            # 2·T[i+1,j] + T[i,j-1] + T[i,j+1] − (4 + 2hΔ/k_A)·T[i,j] = −2(hΔ/k_A)·T_w
            A[r, row(i+1, j  )] += 2.0
            A[r, row(i,   j-1)] += 1.0
            A[r, row(i,   j+1)] += 1.0
            A[r, r             ] += -(4.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

        elif wtype == 'corner_TL':
            # Top-left corner: Δ/2×Δ/2 CV, convection on bottom AND right faces
            # T[iL,j-1] + T[iL-1,j] − (2 + 2hΔ/k_A)·T[iL,j] = −2(hΔ/k_A)·T_w
            A[r, row(i,   j-1)] += 1.0
            A[r, row(i-1, j  )] += 1.0
            A[r, r             ] += -(2.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

        elif wtype == 'corner_TR':
            # Top-right corner: Δ/2×Δ/2 CV, convection on bottom AND left faces
            # T[iR,j-1] + T[iR+1,j] − (2 + 2hΔ/k_A)·T[iR,j] = −2(hΔ/k_A)·T_w
            A[r, row(i,   j-1)] += 1.0
            A[r, row(i+1, j  )] += 1.0
            A[r, r             ] += -(2.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

        elif wtype == 'corner_BL':
            # Bottom-left corner: Δ/2×Δ/2 CV, convection on top AND right faces
            # T[iL,j+1] + T[iL-1,j] − (2 + 2hΔ/k_A)·T[iL,j] = −2(hΔ/k_A)·T_w
            A[r, row(i,   j+1)] += 1.0
            A[r, row(i-1, j  )] += 1.0
            A[r, r             ] += -(2.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

        elif wtype == 'corner_BR':
            # Bottom-right corner: Δ/2×Δ/2 CV, convection on top AND left faces
            # T[iR,j+1] + T[iR+1,j] − (2 + 2hΔ/k_A)·T[iR,j] = −2(hΔ/k_A)·T_w
            A[r, row(i,   j+1)] += 1.0
            A[r, row(i+1, j  )] += 1.0
            A[r, r             ] += -(2.0 + 2.0*hDkA)
            b[r]                 += -2.0*hDkA*T_w

    for idx, (i, j) in enumerate(NODE_LIST):
        r = idx   # row index equals position in NODE_LIST by construction

        # ------------------------------------------------------------------
        # Priority 1: channel wall nodes (override all other classifications)
        # ------------------------------------------------------------------
        wtype = get_wall_type(i, j)
        if wtype is not None:
            _fill_wall(A, b, r, i, j, wtype, hDkA, T_w)
            continue

        # ------------------------------------------------------------------
        # Priority 2: domain corners (quarter-node, two insulated/flux faces)
        # ------------------------------------------------------------------
        if i == 0 and j == 0:
            # Top-left corner: flux q_pp on top, insulated on left
            # T[1,0] + T[0,1] − 2·T[0,0] = −q_pp·Δ/k_P
            A[r, row(1, 0)] += 1.0
            A[r, row(0, 1)] += 1.0
            A[r, r         ] += -2.0
            b[r]              = -q_pp * D / k_P
            continue

        if i == NX-1 and j == 0:
            # Top-right corner: flux q_pp on top, insulated on right
            A[r, row(NX-2, 0)] += 1.0
            A[r, row(NX-1, 1)] += 1.0
            A[r, r             ] += -2.0
            b[r]                 = -q_pp * D / k_P
            continue

        if i == 0 and j == NY-1:
            # Bottom-left corner: insulated on bottom and left
            A[r, row(1,    NY-1)] += 1.0
            A[r, row(0,    NY-2)] += 1.0
            A[r, r               ] += -2.0
            b[r]                   = 0.0
            continue

        if i == NX-1 and j == NY-1:
            # Bottom-right corner: insulated on bottom and right
            A[r, row(NX-2, NY-1)] += 1.0
            A[r, row(NX-1, NY-2)] += 1.0
            A[r, r                ] += -2.0
            b[r]                    = 0.0
            continue

        # ------------------------------------------------------------------
        # Priority 3: top surface interior (j=0, absorbs solar flux)
        # ------------------------------------------------------------------
        if j == 0:
            # Half-CV: flux q_pp on top face, conduction on other three faces
            # T[i-1,0] + T[i+1,0] + 2·T[i,1] − 4·T[i,0] = −2·q_pp·Δ/k_P
            A[r, row(i-1, 0)] += 1.0
            A[r, row(i+1, 0)] += 1.0
            A[r, row(i,   1)] += 2.0
            A[r, r           ] += -4.0
            b[r]               = -2.0 * q_pp * D / k_P
            continue

        # ------------------------------------------------------------------
        # Priority 4: bottom surface interior (j=30, insulated)
        # ------------------------------------------------------------------
        if j == NY-1:
            # Half-CV: insulated bottom face
            # T[i-1,NY-1] + T[i+1,NY-1] + 2·T[i,NY-2] − 4·T[i,NY-1] = 0
            A[r, row(i-1, NY-1)] += 1.0
            A[r, row(i+1, NY-1)] += 1.0
            A[r, row(i,   NY-2)] += 2.0
            A[r, r               ] += -4.0
            b[r]                   = 0.0
            continue

        # ------------------------------------------------------------------
        # Priority 5: left edge (i=0), non-corner
        # ------------------------------------------------------------------
        if i == 0:
            if j == J_INT:
                # Interface left edge: half-CV in x, k-weighted normal fluxes
                # k_P·T[0,9] + k_A·T[0,11] + (k_P+k_A)·T[1,10] − 2(k_P+k_A)·T[0,10] = 0
                A[r, row(0, J_INT-1)] += k_P
                A[r, row(0, J_INT+1)] += k_A
                A[r, row(1, J_INT  )] += (k_P + k_A)
                A[r, r               ] += -2.0 * (k_P + k_A)
                b[r]                   = 0.0
            else:
                # Insulated left face: mirror image gives factor of 2
                # 2·T[1,j] + T[0,j-1] + T[0,j+1] − 4·T[0,j] = 0
                A[r, row(1, j  )] += 2.0
                A[r, row(0, j-1)] += 1.0
                A[r, row(0, j+1)] += 1.0
                A[r, r           ] += -4.0
                b[r]               = 0.0
            continue

        # ------------------------------------------------------------------
        # Priority 6: right edge (i=NX-1), non-corner
        # ------------------------------------------------------------------
        if i == NX-1:
            if j == J_INT:
                # Interface right edge: symmetric to left interface edge
                A[r, row(NX-1, J_INT-1)] += k_P
                A[r, row(NX-1, J_INT+1)] += k_A
                A[r, row(NX-2, J_INT  )] += (k_P + k_A)
                A[r, r                  ] += -2.0 * (k_P + k_A)
                b[r]                      = 0.0
            else:
                # Insulated right face
                A[r, row(NX-2, j  )] += 2.0
                A[r, row(NX-1, j-1)] += 1.0
                A[r, row(NX-1, j+1)] += 1.0
                A[r, r              ] += -4.0
                b[r]                  = 0.0
            continue

        # ------------------------------------------------------------------
        # Priority 7: PV/Al interface interior (j=10, i=1..224)
        # ------------------------------------------------------------------
        if j == J_INT:
            # Full-size CV straddling two materials; normal fluxes use actual k,
            # lateral fluxes use arithmetic mean k_avg for parallel conduction
            # k_P·T[i,9] + k_A·T[i,11] + k_avg·(T[i-1,10]+T[i+1,10]) − 2(k_P+k_A)·T[i,10] = 0
            A[r, row(i,   J_INT-1)] += k_P
            A[r, row(i,   J_INT+1)] += k_A
            A[r, row(i-1, J_INT  )] += k_avg
            A[r, row(i+1, J_INT  )] += k_avg
            A[r, r                 ] += -2.0 * (k_P + k_A)
            b[r]                     = 0.0
            continue

        # ------------------------------------------------------------------
        # Priority 8: fully interior node (PV or aluminum bulk)
        # ------------------------------------------------------------------
        # Standard 5-point Laplacian: sum of four neighbours − 4·T[i,j] = 0
        A[r, row(i-1, j  )] += 1.0
        A[r, row(i+1, j  )] += 1.0
        A[r, row(i,   j-1)] += 1.0
        A[r, row(i,   j+1)] += 1.0
        A[r, r              ] += -4.0
        b[r]                  = 0.0

    return A.tocsr(), b




# ---------------------------------------------------------------------------
# Map the solution vector back to the 2-D (NY × NX) field
# ---------------------------------------------------------------------------
def reconstruct_field(T_vec):
    """
    Parameters
    ----------
    T_vec : 1-D array of length N_NODES

    Returns
    -------
    T_field : (NY, NX) float array — NaN where void
    """
    T_field = np.full((NY, NX), np.nan)
    valid = NODE_MAP >= 0
    T_field[valid] = T_vec[NODE_MAP[valid]]
    return T_field
