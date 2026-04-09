import numpy as np
from scipy.optimize import brentq

# Geometry [m]
R0 = 0.30       # inner cavity radius
T1 = 0.002      # Layer 1
T3 = 0.003      # Layer 3
T2_MIN = 0.001  # aerogel thickness range
T2_MAX = 0.040

def get_radii(t2):
    r1 = R0 + T1
    r2 = r1 + t2
    r3 = r2 + T3
    return R0, r1, r2, r3

# Temperature dependent conductivities
def k1(T):
    return 160.0 * (1.0 - 2.0e-4 *(T - 300.0))

def k2(T):
    dT = T - 300.0
    return 0.035 + 2.0e-4 *dT + 1.0e-7 * dT**2

def k3(T):
    return 6.0 * (1.0 + 1.0e-4 * (T - 300.0))

# Contact resistances [m^2*K/W]
RC12_PP = 2.0e-4  # Al / Aerogel at r1
RC23_PP = 1.0e-4  # Aerogel / CFRP at r2

Q_INT_DEFAULT = 50.0  # electronics waste heat [W]

# Inner boundary: given gas propertise
KG = 0.030
NU_GAS= 1.7e-5
ALPHA_GAS = 2.4e-5
BETA = 1.0 / 300.0
LC = R0
G_EFF = 0.02

def h_inner(Tg, Tw):
    Ra = G_EFF * BETA * abs(Tg - Tw) * LC**3 / (NU_GAS * ALPHA_GAS)
    Nu = 2.0 + 0.55 * Ra**0.25
    return Nu * KG / LC

# Outer boundary: radiation
SIGMA = 5.670e-8
T_SPACE = 3.0
G_SOL = 9100.0
G_SOL_FAC  = 0.25
G_MERC= 1200.0
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

def solve(t2, Q_int=Q_INT_DEFAULT, N=20, tol=1e-6, max_iter=500, omega=1.0):
    r0, r1, r2, r3 = get_radii(t2)
    Q = Q_int

    # Outer surface temperature from radiation balance
    Ts = brentq(lambda T: Q_rad_net(T, r3) - Q, 100.0, 1500.0)

    # FDM node arrays for each layer
    r_L1 = np.linspace(r0, r1, N)
    r_L2 = np.linspace(r1, r2, N)
    r_L3 = np.linspace(r2, r3, N)
    dr1 = r_L1[1] - r_L1[0]
    dr2 = r_L2[1] - r_L2[0]
    dr3 = r_L3[1] - r_L3[0]

    # Unknowns: [Tg, T1[0..N-1], T2[0..N-1], T3[0..N-1]]
    n_total = 1 + 3 * N
    i1 = 1
    i2 = N + 1
    i3 = 2 * N + 1

    T_vec = np.zeros(n_total)
    T_vec[0] = Ts + 40
    T_vec[i1:i1+N] = np.linspace(Ts + 30, Ts + 20, N)
    T_vec[i2:i2+N] = np.linspace(Ts + 20, Ts + 5,  N)
    T_vec[i3:i3+N] = np.linspace(Ts + 5,  Ts,      N)

    converged = False
    for iteration in range(max_iter):
        T_old = T_vec.copy()

        Tg_c = T_old[0]
        T1c = T_old[i1:i1+N]
        T2c = T_old[i2:i2+N]
        T3c = T_old[i3:i3+N]

        A = np.zeros((n_total, n_total))
        b = np.zeros(n_total)

        # cavity gas energy balance
        hi = h_inner(Tg_c, T1c[0])
        a0 = hi * r0**2
        A[0, 0]  =  a0
        A[0, i1] = -a0
        b[0] = Q_int / (4.0 * np.pi)

        # layer 1
        r_ph = 0.5 * (r_L1[0] + r_L1[1])
        k_ph = k1(0.5 * (T1c[0] + T1c[1]))
        cp = r_ph**2 * k_ph / dr1
        A[i1, 0] = a0
        A[i1, i1] = -(a0 + cp)
        A[i1, i1+1] = cp

        # Interior nodes
        for i in range(1, N - 1):
            r_ph = 0.5 * (r_L1[i] + r_L1[i+1])
            r_mh = 0.5 * (r_L1[i] + r_L1[i-1])
            kp = k1(0.5 * (T1c[i] + T1c[i+1]))
            km = k1(0.5 * (T1c[i] + T1c[i-1]))
            cp = r_ph**2 * kp / dr1
            cm = r_mh**2 * km / dr1
            row = i1 + i
            A[row, row - 1] = cm
            A[row, row] = -(cm + cp)
            A[row, row + 1] = cp

        # Node N-1 at r1
        r_mh = 0.5 * (r_L1[N-2] + r_L1[N-1])
        km = k1(0.5 * (T1c[N-2] + T1c[N-1]))
        cm = r_mh**2 * km / dr1
        cr = r1**2 / RC12_PP
        row = i1 + N - 1
        A[row, row - 1] =cm
        A[row, row] = -(cm + cr)
        A[row, i2] =cr

        # Layer 2
        r_ph = 0.5 * (r_L2[0] + r_L2[1])
        kp = k2(0.5 * (T2c[0] + T2c[1]))
        cp = r_ph**2 * kp / dr2
        cr = r1**2/RC12_PP
        A[i2, i1 + N - 1] =  cr
        A[i2, i2]  = -(cr + cp)
        A[i2, i2 + 1] = cp

        # Interior nodes
        for i in range(1, N - 1):
            r_ph = 0.5 * (r_L2[i] + r_L2[i+1])
            r_mh = 0.5 * (r_L2[i] + r_L2[i-1])
            kp = k2(0.5* (T2c[i] + T2c[i+1]))
            km = k2(0.5 * (T2c[i] + T2c[i-1]))
            cp = r_ph**2 * kp / dr2
            cm = r_mh**2 * km / dr2
            row = i2 + i
            A[row, row -1] =  cm
            A[row, row] = -(cm + cp)
            A[row, row + 1] =  cp

        # Node N-1 at r2
        r_mh = 0.5 * (r_L2[N-2] + r_L2[N-1])
        km = k2(0.5* (T2c[N-2] + T2c[N-1]))
        cm = r_mh**2 * km / dr2
        cr = r2**2 / RC23_PP
        row = i2 + N - 1
        A[row, row -1] =  cm
        A[row, row] = -(cm + cr)
        A[row, i3] = cr

        # Layer 3 (CFRP)
        r_ph = 0.5 * (r_L3[0] + r_L3[1])
        kp = k3(0.5 * (T3c[0] + T3c[1]))
        cp = r_ph**2 * kp / dr3
        cr = r2**2/RC23_PP
        A[i3, i2 + N-1] = cr
        A[i3, i3] = -(cr + cp)
        A[i3, i3 + 1] = cp

        # Interior nodes
        for i in range(1, N - 1):
            r_ph = 0.5 * (r_L3[i] + r_L3[i+1])
            r_mh = 0.5 * (r_L3[i] + r_L3[i-1])
            kp = k3(0.5 * (T3c[i] + T3c[i+1]))
            km = k3(0.5 * (T3c[i] + T3c[i-1]))
            cp = r_ph**2 * kp / dr3
            cm = r_mh**2 * km / dr3
            row = i3 + i
            A[row, row -1] = cm
            A[row, row] = -(cm+ cp)
            A[row, row + 1] = cp

        # Node N-1 at r3; Ts
        row =i3 + N - 1
        A[row, row] = 1.0
        b[row] = Ts

        T_new =np.linalg.solve(A, b)
        T_vec = omega * T_new + (1.0 - omega) * T_old

        if np.max(np.abs(T_vec -T_old)) < tol:
            converged = True
            break

    if not converged:
        print(f"warning: Picard did not converge in {max_iter} iterations ")

    Tg = T_vec[0]
    T1_arr = T_vec[i1:i1+N]
    T2_arr = T_vec[i2:i2+N]
    T3_arr = T_vec[i3:i3+N]

    Tw  = T1_arr[0]
    T1m = T1_arr[-1]
    T1p = T2_arr[0]
    T2m = T2_arr[-1]
    T2p = T3_arr[0]

    dTc12 = T1m - T1p
    dTc23 = T2m- T2p

    r_profile = np.concatenate([r_L1, r_L2, r_L3])
    T_profile = np.concatenate([T1_arr, T2_arr, T3_arr])

    return {
        'Tg': Tg, 'Tw': Tw,
        'T1m': T1m, 'T1p': T1p,
        'T2m': T2m, 'T2p': T2p,
        'Ts': T3_arr[-1], 'Q': Q,
        'dTc12': dTc12, 'dTc23': dTc23,
        'r0': r0, 'r1': r1, 'r2': r2, 'r3': r3,
        't2': t2, 'Q_int': Q_int,
        'r_profile': r_profile, 'T_profile': T_profile, 'N': N,
        'iterations': iteration + 1,
    }

def temperature_profile(sol):
    return sol['r_profile'], sol['T_profile']

if __name__ == '__main__':
    sol = solve(t2=0.010)
    
