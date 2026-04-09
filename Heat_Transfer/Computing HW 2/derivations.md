# CP2 Finite Difference Derivations

## Grid Parameters

- Domain: 0.45 m (x) × 0.06 m (y)
- Nodes: 226 (x, i=0…225) × 31 (y, j=0…30)
- **Δx = Δy = Δ = 0.002 m** (uniform grid)
- i increases left → right, j increases top → bottom

### Material Regions

| Region | j range | y range | k |
|--------|---------|---------|---|
| PV panel | j = 0…10 | 0 – 0.02 m | k_P = 1 W/m·K |
| Aluminum | j = 10…30 | 0.02 – 0.06 m | k_A = 205 W/m·K |
| Interface | j = 10 | y = 0.02 m | both |

### Water Channel Positions

Each channel: width 0.01 m, height 0.02 m, spanning **j = 15…25** (y = 0.03–0.05 m)

| Channel | i_L | i_R | x range |
|---------|-----|-----|---------|
| 1 | 50 | 55 | 0.10–0.11 m |
| 2 | 80 | 85 | 0.16–0.17 m |
| 3 | 110 | 115 | 0.22–0.23 m |
| 4 | 140 | 145 | 0.28–0.29 m |
| 5 | 170 | 175 | 0.34–0.35 m |

**No nodes are placed inside the channel void.** Only wall nodes at j=15, j=25, i=i_L, i=i_R.

---

## Parameters (all from problem statement)

- k_P = 1 W/m·K (PV panel)
- k_A = 205 W/m·K (aluminum)
- h = 500 W/m²·K (convection in channels)
- q'' = 1500 W/m² (radiation heat flux at top surface)
- T_w = water temperature (20°C for Part 2)

---

## Finite Difference Equations (Energy Balance Method)

### Type 1 — Interior node (PV or aluminum, uniform material)

Control volume: Δ × Δ

```
q_top    = k · Δ · (T[i,j-1] − T[i,j]) / Δ
q_bottom = k · Δ · (T[i,j+1] − T[i,j]) / Δ
q_left   = k · Δ · (T[i-1,j] − T[i,j]) / Δ
q_right  = k · Δ · (T[i+1,j] − T[i,j]) / Δ
```

Sum = 0, divide by k (k cancels — same form for both materials):

**T[i-1,j] + T[i+1,j] + T[i,j-1] + T[i,j+1] − 4·T[i,j] = 0**

---

### Type 2 — Top surface interior node (j=0, i=1…224)

Control volume: Δ × Δ/2 (half CV, radiation input at top face)

```
q_rad    = q'' · Δ                                    (heat in)
q_bottom = k_P · Δ · (T[i,1] − T[i,0]) / Δ
q_left   = k_P · (Δ/2) · (T[i-1,0] − T[i,0]) / Δ
q_right  = k_P · (Δ/2) · (T[i+1,0] − T[i,0]) / Δ
```

Sum = 0, divide by k_P/2:

**T[i-1,0] + T[i+1,0] + 2·T[i,1] − 4·T[i,0] = −2q''Δ/k_P**

---

### Type 3 — Top-left corner (j=0, i=0)

Control volume: Δ/2 × Δ/2 (radiation at top, insulated at left)

```
q_rad    = q'' · (Δ/2)
q_bottom = k_P · (Δ/2) · (T[0,1] − T[0,0]) / Δ
q_right  = k_P · (Δ/2) · (T[1,0] − T[0,0]) / Δ
```

Sum = 0, divide by k_P/2:

**T[1,0] + T[0,1] − 2·T[0,0] = −q''Δ/k_P**

Top-right corner (j=0, i=225) by symmetry:

**T[224,0] + T[225,1] − 2·T[225,0] = −q''Δ/k_P**

---

### Type 4 — Left side interior node (i=0, j in PV or Al interior, not interface)

Control volume: Δ/2 × Δ (insulated at left face)

```
q_right  = k · Δ · (T[1,j] − T[0,j]) / Δ
q_top    = k · (Δ/2) · (T[0,j-1] − T[0,j]) / Δ
q_bottom = k · (Δ/2) · (T[0,j+1] − T[0,j]) / Δ
```

Sum = 0, divide by k/2:

**2·T[1,j] + T[0,j-1] + T[0,j+1] − 4·T[0,j] = 0**

Right side (i=225) by symmetry:

**2·T[224,j] + T[225,j-1] + T[225,j+1] − 4·T[225,j] = 0**

---

### Type 5 — Bottom interior node (j=30, i=1…224)

Control volume: Δ × Δ/2 (insulated at bottom face)

```
q_top   = k_A · Δ · (T[i,29] − T[i,30]) / Δ
q_left  = k_A · (Δ/2) · (T[i-1,30] − T[i,30]) / Δ
q_right = k_A · (Δ/2) · (T[i+1,30] − T[i,30]) / Δ
```

Sum = 0, divide by k_A/2:

**T[i-1,30] + T[i+1,30] + 2·T[i,29] − 4·T[i,30] = 0**

Bottom-left corner (j=30, i=0):

**T[1,30] + T[0,29] − 2·T[0,30] = 0**

Bottom-right corner (j=30, i=225):

**T[224,30] + T[225,29] − 2·T[225,30] = 0**

---

### Type 6 — PV/Al interface interior node (j=10, i=1…224) ★

Control volume: Δ × Δ (upper Δ/2 is PV, lower Δ/2 is aluminum)

**y-direction:** materials are distinct on each side — use k_P and k_A directly (no averaging):
```
q_top    = k_P · Δ · (T[i,9]  − T[i,10]) / Δ
q_bottom = k_A · Δ · (T[i,11] − T[i,10]) / Δ
```

**x-direction:** lateral face spans Δ/2 of PV + Δ/2 of Al in parallel → arithmetic mean:
```
q_left  = (k_P+k_A)/2 · Δ · (T[i-1,10] − T[i,10]) / Δ
q_right = (k_P+k_A)/2 · Δ · (T[i+1,10] − T[i,10]) / Δ
```

Sum = 0, divide by Δ:

**k_P·T[i,9] + k_A·T[i,11] + (k_P+k_A)/2·(T[i-1,10] + T[i+1,10]) − 2(k_P+k_A)·T[i,10] = 0**

Interface left edge (j=10, i=0), left side insulated:

**(k_P+k_A)·T[1,10] + k_P·T[0,9] + k_A·T[0,11] − 2(k_P+k_A)·T[0,10] = 0**

Interface right edge (j=10, i=225) by symmetry:

**(k_P+k_A)·T[224,10] + k_P·T[225,9] + k_A·T[225,11] − 2(k_P+k_A)·T[225,10] = 0**

Note: harmonic mean is NOT needed here because the interface node sits exactly at the material boundary (j=10 = y=0.02m). The y-direction fluxes use k_P and k_A directly without any averaging.

---

### Type 7 — Channel top wall interior node (j=15, i=i_L+1…i_R−1) ★

Control volume: Δ × Δ/2 (only extends upward; bottom face exposed to water)

```
q_top    = k_A · Δ · (T[i,14] − T[i,15]) / Δ
q_conv   = h · Δ · (T_w − T[i,15])               (convection, bottom face)
q_left   = k_A · (Δ/2) · (T[i-1,15] − T[i,15]) / Δ
q_right  = k_A · (Δ/2) · (T[i+1,15] − T[i,15]) / Δ
```

Sum = 0, divide by k_A/2:

**T[i-1,15] + T[i+1,15] + 2·T[i,14] − (4 + 2hΔ/k_A)·T[i,15] = −2(hΔ/k_A)·T_w**

---

### Type 8 — Channel bottom wall interior node (j=25, i=i_L+1…i_R−1) ★

Control volume: Δ × Δ/2 (only extends downward; top face exposed to water)

```
q_bottom = k_A · Δ · (T[i,26] − T[i,25]) / Δ
q_conv   = h · Δ · (T_w − T[i,25])               (convection, top face)
q_left   = k_A · (Δ/2) · (T[i-1,25] − T[i,25]) / Δ
q_right  = k_A · (Δ/2) · (T[i+1,25] − T[i,25]) / Δ
```

Sum = 0, divide by k_A/2:

**T[i-1,25] + T[i+1,25] + 2·T[i,26] − (4 + 2hΔ/k_A)·T[i,25] = −2(hΔ/k_A)·T_w**

---

### Type 9 — Channel left wall interior node (i=i_L, j=16…24) ★

Control volume: Δ/2 × Δ (only extends left; right face exposed to water)

```
q_left   = k_A · Δ · (T[i_L-1,j] − T[i_L,j]) / Δ
q_conv   = h · Δ · (T_w − T[i_L,j])              (convection, right face)
q_top    = k_A · (Δ/2) · (T[i_L,j-1] − T[i_L,j]) / Δ
q_bottom = k_A · (Δ/2) · (T[i_L,j+1] − T[i_L,j]) / Δ
```

Sum = 0, divide by k_A/2:

**2·T[i_L-1,j] + T[i_L,j-1] + T[i_L,j+1] − (4 + 2hΔ/k_A)·T[i_L,j] = −2(hΔ/k_A)·T_w**

Channel right wall interior (i=i_R, j=16…24), right face exposed to water:

**2·T[i_R+1,j] + T[i_R,j-1] + T[i_R,j+1] − (4 + 2hΔ/k_A)·T[i_R,j] = −2(hΔ/k_A)·T_w**

---

### Type 10 — Channel corner nodes ★

Control volume: Δ/2 × Δ/2 (two faces exposed to water convection)

General form for all four corners:

**N1 + N2 − (2 + 2hΔ/k_A)·T_corner = −2(hΔ/k_A)·T_w**

where N1 and N2 are the two solid neighbors:

| Corner | N1 (vertical neighbor) | N2 (horizontal neighbor) |
|--------|------------------------|--------------------------|
| Top-left  (j=15, i=i_L) | T[i_L, 14]  | T[i_L-1, 15] |
| Top-right (j=15, i=i_R) | T[i_R, 14]  | T[i_R+1, 15] |
| Bot-left  (j=25, i=i_L) | T[i_L, 26]  | T[i_L-1, 25] |
| Bot-right (j=25, i=i_R) | T[i_R, 26]  | T[i_R+1, 25] |

Derivation (top-left corner as example):
```
q_top    = k_A · (Δ/2) · (T[i_L,14]   − T[i_L,15]) / Δ
q_left   = k_A · (Δ/2) · (T[i_L-1,15] − T[i_L,15]) / Δ
q_conv_1 = h · (Δ/2) · (T_w − T[i_L,15])    (bottom face)
q_conv_2 = h · (Δ/2) · (T_w − T[i_L,15])    (right face)
```
Sum = 0, divide by k_A/2.

---

## Summary Table

| Type | Location | RHS |
|------|----------|-----|
| 1 | Interior (PV or Al) | 0 |
| 2 | Top surface interior | −2q''Δ/k_P |
| 3 | Top corners | −q''Δ/k_P |
| 4 | Left/right side | 0 |
| 5 | Bottom surface + corners | 0 |
| 6 | PV/Al interface | 0 |
| 7 | Channel top wall | −2(hΔ/k_A)·T_w |
| 8 | Channel bottom wall | −2(hΔ/k_A)·T_w |
| 9 | Channel left/right wall | −2(hΔ/k_A)·T_w |
| 10 | Channel corners | −2(hΔ/k_A)·T_w |
