# Heat Transfer CP1 — Satellite Thermal Design (Mercury Orbit)

**Course:** 530.334 Heat Transfer — Spring 2026
**Due:** Gradescope, 3/6/26 10am
**Topic:** Steady-state, 1D radial, nonlinear conduction + radiation

---

## Problem Overview

A satellite orbits Mercury modeled as a **concentric 3-layer spherical shell**.

Heat path:
1. Electronics generate waste heat → cavity gas (natural convection)
2. Gas → inner wall (r₀) → radial conduction through 3 layers
3. Outer surface (r₃) → radiation to deep space (also absorbs solar + Mercury IR)

---

## Geometry

| Symbol | Value | Description |
|--------|-------|-------------|
| r₀ | 0.30 m | Inner cavity radius |
| t₁ | 2 mm | Layer 1: Al alloy liner (r₀ → r₁) |
| t₂ | **1–40 mm** (design var) | Layer 2: Aerogel insulation (r₁ → r₂) |
| t₃ | 3 mm | Layer 3: CFRP outer shell (r₂ → r₃) |

- r₁ = r₀ + t₁
- r₂ = r₁ + t₂
- r₃ = r₂ + t₃

---

## Parameters

### Internal Heat Source
- **Q_int = 50 W** (electronics waste heat)
  - Note: figure label in PDF shows 220 W — use text value 50 W
- Cavity gas temperature T_g (unknown, uniform)

### Natural Convection (inner boundary, nonlinear h_i)

$$Ra = \frac{g_\text{eff} \beta |T_g - T_w| L_c^3}{\nu \alpha}$$

$$Nu = 2 + 0.55 \, Ra^{1/4}$$

$$h_i(T_g, T_w) = \frac{Nu \cdot k_g}{L_c}$$

| Property | Value |
|----------|-------|
| k_g | 0.030 W/(m·K) |
| ν | 1.7 × 10⁻⁵ m²/s |
| α | 2.4 × 10⁻⁵ m²/s |
| β | 1/300 K⁻¹ |
| L_c | r₀ = 0.30 m |
| g_eff | 0.02 m/s² |

### Contact Resistances (area-specific)

| Interface | Location | R''_c (m²K/W) |
|-----------|----------|----------------|
| Layer 1 / Layer 2 | r = r₁ | 2.0 × 10⁻⁴ |
| Layer 2 / Layer 3 | r = r₂ | 1.0 × 10⁻⁴ |

Temperature jump at interface:
$$\Delta T_c = Q \cdot \frac{R''_c}{4\pi r^2}$$

### Temperature-Dependent Conductivities

$$k_1(T) = 160[1 - 2.0\times10^{-4}(T-300)] \quad \text{(Al, W/m·K)}$$

$$k_2(T) = 0.035 + 2.0\times10^{-4}(T-300) + 1.0\times10^{-7}(T-300)^2 \quad \text{(Aerogel)}$$

$$k_3(T) = 6.0[1 + 1.0\times10^{-4}(T-300)] \quad \text{(CFRP)}$$

### Outer Surface Radiation Boundary

| Parameter | Value |
|-----------|-------|
| T_space | 3 K |
| G_☉ (solar irradiance) | 9100 W/m² |
| g_☉ (geometric factor) | 0.25 |
| G_M (Mercury IR) | 1200 W/m² |
| g_m (geometric factor) | 0.25 |

Temperature-dependent optical properties of outer surface:

$$\alpha(T) = 0.10 + 1.5\times10^{-4}(T-300)$$

$$\varepsilon(T) = 0.92 - 1.0\times10^{-4}(T-300)$$

Outer surface energy balance (T_s = T(r₃)):

$$Q = 4\pi r_3^2 \left[ \alpha(T_s)(g_\odot G_\odot + g_m G_M) - \varepsilon(T_s)\sigma(T_s^4 - T_\text{space}^4) \right]$$

(σ = 5.67 × 10⁻⁸ W/m²K⁴)

---

## Design Requirements

| Constraint | Limit |
|-----------|-------|
| Cavity gas temp T_g | ≤ 360 K |
| Outer surface temp T_s | ≤ 450 K |

---

## Tasks Summary

### (A) Model Formulation
- Governing PDE: 1D spherical steady-state with k(T)
- Interface BCs: heat flux continuity + temperature jump from contact resistance
- Inner BC: convection coupling with h_i(T_g, T_w) + cavity energy balance
- Outer BC: nonlinear radiation balance
- **Submit:** equation set, unknowns list, numerical approach description, dominant resistance discussion

### (B) Baseline Case — t₂ = 10 mm
- Solve for: T_g, T_w = T(r₀), T_s = T(r₃), total Q
- Interface temperature jumps: ΔT_c,12 and ΔT_c,23
- **Submit:** numerical values with units + pass/fail for both design requirements

### (C) Maximum Insulation Thickness t₂
- Sweep t₂ ∈ [1, 40] mm
- Find maximum t₂ such that T_g ≤ 360 K AND T_s ≤ 450 K
- **Submit:** Plot T_g(t₂) with 360 K line, Plot T_s(t₂) with 450 K line, t₂_max value

### (D) Radial Temperature Distribution — use t₂_max from (C)
- Plot T(r) from r₀ to r₃, labeled by layer
- Show temperature discontinuities at interfaces
- **Submit:** T(r) plot + two interface jump magnitudes (K)

### (E) Sensitivity Study — fix t₂ = t₂_max
- Sweep Q_int = 20–200 W
- Plot T_g(Q_int), find Q_int,max such that T_g = 360 K
- **Submit:** plot + Q_int,max value

---

## Numerical Approach (recommended)

- **Discretization:** Finite difference or finite volume, 10–30 nodes per layer
- **Nonlinearities:** T⁴, k(T), α(T), ε(T), h_i(T_g, T_w) → use Picard iteration or `fsolve`
- **Convergence criterion:** max ΔT < 10⁻⁴ K or energy imbalance < 0.1%

### Unknowns
- T(r) at all nodes
- T_g (cavity gas temperature)
- Q (total heat flow, must be consistent throughout)

---

## Key Physics Notes

- Aerogel has very low k (~0.035 W/m·K) → expected to dominate thermal resistance
- Al liner has high k (~160 W/m·K) → negligible resistance
- CFRP has moderate k (~6 W/m·K) → minor resistance
- Contact resistances act as additional jumps at interfaces

---

## Submission Checklist

- [ ] (A) Equations, unknowns, method description, dominant resistance paragraph
- [ ] (B) T_g, T_w, T_s, Q values + pass/fail
- [ ] (C) Two plots + t₂_max
- [ ] (D) T(r) plot + jump magnitudes
- [ ] (E) T_g(Q_int) plot + Q_int,max
- [ ] Commented code for all parts
- [ ] One styled report on Gradescope
