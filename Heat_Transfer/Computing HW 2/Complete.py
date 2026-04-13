import numpy as np
import matplotlib.pyplot as plt

# ======================
# 1. PARAMETERS
# ======================
Nx, Ny = 226, 31
Lx, Ly = 0.45, 0.06

dx = Lx / (Nx - 1)
dy = dx

N = Nx * Ny
Tw = 20.0

k_pv = 1.0
k_al = 205.0
q_rad = 1500.0
h = 500.0

q_term = q_rad * dx / k_pv
Bi = h * dx / k_al

# ======================
# INDEX
# ======================
def idx(i, j):
    return j * Nx + i

# ======================
# NODE MAP
# ======================
NodeMap = np.zeros((Nx, Ny), dtype=int)

for i in range(Nx):
    for j in range(Ny):

        if j == 0:
            if i == 0: NodeMap[i,j] = 1
            elif i == Nx-1: NodeMap[i,j] = 3
            else: NodeMap[i,j] = 2

        elif j < 5:
            if i == 0: NodeMap[i,j] = 4
            elif i == Nx-1: NodeMap[i,j] = 6
            else: NodeMap[i,j] = 5

        elif j == 5:
            if i == 0: NodeMap[i,j] = 7
            elif i == Nx-1: NodeMap[i,j] = 9
            else: NodeMap[i,j] = 8

        elif j < Ny-1:
            if i == 0: NodeMap[i,j] = 10
            elif i == Nx-1: NodeMap[i,j] = 12
            else: NodeMap[i,j] = 11

        elif j == Ny-1:
            if i == 0: NodeMap[i,j] = 14
            elif i == Nx-1: NodeMap[i,j] = 15
            else: NodeMap[i,j] = 13

# ======================
# CHANNELS (FIXED)
# ======================
for x in np.arange(0.10, 0.35, 0.06):
    i_ws = int(round(x/dx))
    i_we = int(round((x+0.01)/dx))

    j_ws = int(round(0.02/dy))
    j_we = int(round(0.04/dy))

    # FIXED HERE
    NodeMap[i_ws+1:i_we-1, j_ws+1:j_we-1] = 99

    NodeMap[i_ws+1:i_we-1, j_ws] = 16
    NodeMap[i_ws+1:i_we-1, j_we] = 17
    NodeMap[i_ws, j_ws+1:j_we-1] = 18
    NodeMap[i_we, j_ws+1:j_we-1] = 19

    NodeMap[i_ws, j_ws] = 20
    NodeMap[i_we, j_ws] = 21
    NodeMap[i_ws, j_we] = 22
    NodeMap[i_we, j_we] = 23

# ======================
# MATRIX
# ======================
A = np.zeros((N, N))
b = np.zeros(N)

for i in range(Nx):
    for j in range(Ny):

        p = idx(i,j)
        t = NodeMap[i,j]

        def I(ii, jj):
            return idx(ii, jj)

        if t == 1:
            A[p,p]=-2; A[p,I(i+1,j)]=1; A[p,I(i,j+1)]=1; b[p]=-q_term

        elif t == 2:
            A[p,p]=-4; A[p,I(i-1,j)]=1; A[p,I(i+1,j)]=1
            A[p,I(i,j+1)]=2; b[p]=-2*q_term

        elif t == 3:
            A[p,p]=-2; A[p,I(i-1,j)]=1; A[p,I(i,j+1)]=1; b[p]=-q_term

        elif t == 4:
            A[p,p]=-4; A[p,I(i+1,j)]=2; A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1

        elif t == 5:
            A[p,p]=-4; A[p,I(i-1,j)]=1; A[p,I(i+1,j)]=1
            A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1

        elif t == 6:
            A[p,p]=-4; A[p,I(i-1,j)]=2; A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1

        # ===== INTERFACE FIX =====
        elif t in [7,8,9]:
            A[p,p] = -2*(k_pv + k_al)
            A[p,I(i-1,j)] = (k_pv + k_al)/2
            A[p,I(i+1,j)] = (k_pv + k_al)/2
            A[p,I(i,j-1)] = k_pv
            A[p,I(i,j+1)] = k_al

        elif t == 10:
            A[p,p]=-4; A[p,I(i+1,j)]=2; A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1

        elif t == 11:
            A[p,p]=-4; A[p,I(i-1,j)]=1; A[p,I(i+1,j)]=1
            A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1

        elif t == 12:
            A[p,p]=-4; A[p,I(i-1,j)]=2; A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1

        elif t == 13:
            A[p,p]=-4; A[p,I(i-1,j)]=1; A[p,I(i+1,j)]=1; A[p,I(i,j-1)]=2

        elif t == 14:
            A[p,p]=-2; A[p,I(i+1,j)]=1; A[p,I(i,j-1)]=1

        elif t == 15:
            A[p,p]=-2; A[p,I(i-1,j)]=1; A[p,I(i,j-1)]=1

        elif t == 16:
            A[p,p]=-2*(Bi+2)
            A[p,I(i-1,j)]=1; A[p,I(i+1,j)]=1
            A[p,I(i,j-1)]=2
            b[p]=-2*Bi*Tw

        elif t == 17:
            A[p,p]=-2*(Bi+2)
            A[p,I(i-1,j)]=1; A[p,I(i+1,j)]=1
            A[p,I(i,j+1)]=2
            b[p]=-2*Bi*Tw

        elif t == 18:
            A[p,p]=-2*(Bi+2)
            A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1
            A[p,I(i-1,j)]=2
            b[p]=-2*Bi*Tw

        elif t == 19:
            A[p,p]=-2*(Bi+2)
            A[p,I(i,j-1)]=1; A[p,I(i,j+1)]=1
            A[p,I(i+1,j)]=2
            b[p]=-2*Bi*Tw

        # ===== CORNER FIX =====
        elif t in [20,21,22,23]:
            A[p,p]=-2*(Bi+3)
            A[p,I(i-1,j)]=2; A[p,I(i+1,j)]=1
            A[p,I(i,j-1)]=2; A[p,I(i,j+1)]=1
            b[p]=-2*Bi*Tw

        elif t == 99:
            A[p,p]=1; b[p]=Tw

# ======================
# SOLVE
# ======================
T = np.linalg.solve(A, b)
T = T.reshape((Ny, Nx))

# ======================
# OUTPUT (你要的格式)
# ======================
Tmax = np.max(T)
idx_max = np.unravel_index(np.argmax(T), T.shape)

x_max = idx_max[1] * dx
y_max = idx_max[0] * dy

print('--- Problem 2 ---')
print(f'Max Temp at 20C water: {Tmax:.2f} C')
print(f'Location: x = {x_max:.3f} m, y = {y_max:.3f} m\n')

delta_T = Tmax - Tw
Tw_opt = 35 - delta_T

print('--- Problem 3 ---')
print(f'To keep panel below 35C, max water temp is: {Tw_opt:.2f} C')

# ======================
# PLOT
# ======================
plt.figure(figsize=(10,3))
plt.contourf(np.linspace(0,0.45,Nx),
             np.linspace(0,0.06,Ny),
             T, 50)
plt.colorbar()
plt.title("Temperature Distribution")
plt.show()