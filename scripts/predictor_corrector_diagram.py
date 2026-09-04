"""
Visualises the predictor-corrector (PC) integration of the Bateman equations.

Model: single-nuclide depletion  dN/dt = -k·N  →  N(t) = N₀ e^{-kt}
  k  =  effective one-group reaction rate  (σ_eff · φ, normalised to 1)

Steps illustrated
-----------------
Predictor  – forward-Euler using the reaction rate computed at t = 0:
               N_pred = N₀ + (-k N₀) Δt

Corrector  – forward-Euler using the rate re-evaluated at the half-step
             composition predicted by the predictor (midpoint method):
               N_mid  = N₀ + (-k N₀) Δt/2         # half-step predictor
               N_corr = N₀ + (-k N_mid) Δt          # corrector
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Parameters ─────────────────────────────────────────────────────────────────
k  = 0.10   # effective reaction rate [1/EFPD]
N0 = 1.0    # normalised initial number density
dt = 8.0    # depletion step [EFPD]  — gives Δt/2 = 4, both on integer ticks

# ── True (reference) solution ──────────────────────────────────────────────────
t_fine     = np.linspace(0, dt, 1000)
N_true     = N0 * np.exp(-k * t_fine)
N_end_true = N0 * np.exp(-k * dt)          # ≈ 0.4966

# ── Predictor  (forward-Euler, initial rates) ──────────────────────────────────
slope_0 = -k * N0                           # dN/dt|t=0 = -k N₀ = -0.10
N_pred  = N0 + slope_0 * dt                 # 1 - 0.70 = 0.300

# ── Half-step midpoint estimate ────────────────────────────────────────────────
t_mid     = dt / 2.0                        # = 3.5 EFPD
N_mid     = N0 + slope_0 * t_mid            # half-step Euler  →  0.650
slope_mid = -k * N_mid                      # reaction rate at midpoint  →  -0.065

# ── Corrector  (midpoint rate for full step) ───────────────────────────────────
N_corr = N0 + slope_mid * dt               # 1 - 0.455 = 0.545

# ── Console summary ────────────────────────────────────────────────────────────
print(f"True N(Δt)      = {N_end_true:.4f}")
print(f"Predictor N(Δt) = {N_pred:.4f}  (error = {N_pred - N_end_true:+.4f})")
print(f"Corrector N(Δt) = {N_corr:.4f}  (error = {N_corr - N_end_true:+.4f})")
print(f"Error reduction : {abs(N_pred-N_end_true)/abs(N_corr-N_end_true):.1f}×")

# ── Colour palette ─────────────────────────────────────────────────────────────
c_true  = '#1565C0'   # dark blue
c_pred  = '#C62828'   # dark red
c_mid   = '#E65100'   # orange
c_corr  = '#2E7D32'   # dark green
c_bound = '#6A1B9A'   # purple for step-boundary lines

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 6.5))

# True solution
ax.plot(t_fine, N_true, color=c_true, lw=2.8, zorder=5,
        label=r'True solution  $N(t) = N_0\,e^{-kt}$')

# Predictor line
t_step = np.array([0.0, dt])
ax.plot(t_step, N0 + slope_0 * t_step,
        color=c_pred, lw=2.1, ls='--', zorder=4,
        label=r'Predictor  – initial rate  $k N_0$')

# Corrector line
ax.plot(t_step, N0 + slope_mid * t_step,
        color=c_corr, lw=2.1, ls='-.', zorder=4,
        label=r'Corrector  – midpoint rate  $k\,N(\Delta t/2)$')

# ── Key points ─────────────────────────────────────────────────────────────────
# Starting point
ax.scatter([0], [N0], s=90, color='black', zorder=10)
ax.text(-0.18, N0, '$N_0$', ha='right', va='center', fontsize=13, fontweight='bold')

# Midpoint composition (half-step predictor)
ax.scatter([t_mid], [N_mid], s=100, color=c_mid, zorder=10, marker='D',
           label=r'$N(\Delta t/2)$ – midpoint composition')
ax.axvline(t_mid, color=c_bound, ls='--', lw=1.2, alpha=0.7)
ax.text(t_mid - 0.12, 1.085, r'$\Delta t/2$', ha='right', va='bottom',
        fontsize=10.5, color=c_bound, fontweight='bold')

# Step boundary
ax.axvline(dt, color=c_bound, ls='--', lw=1.2, alpha=0.7)
ax.text(dt - 0.12, 1.085, r'$\Delta t$', ha='right', va='bottom',
        fontsize=10.5, color=c_bound, fontweight='bold')

# True endpoint
ax.scatter([dt], [N_end_true], s=130, color=c_true, zorder=10, marker='*')
ax.text(dt + 0.12, N_end_true + 0.008,
        f'$N_{{\\mathrm{{true}}}}$ = {N_end_true:.3f}',
        color=c_true, fontsize=10.5, va='bottom')

# Predictor endpoint
ax.scatter([dt], [N_pred], s=90, color=c_pred, zorder=10, marker='s')
ax.text(dt + 0.12, N_pred - 0.012,
        f'$N_{{\\mathrm{{pred}}}}$ = {N_pred:.3f}',
        color=c_pred, fontsize=10.5, va='top')

# Corrector endpoint
ax.scatter([dt], [N_corr], s=90, color=c_corr, zorder=10, marker='^')
ax.text(dt + 0.12, N_corr + 0.008,
        f'$N_{{\\mathrm{{corr}}}}$ = {N_corr:.3f}',
        color=c_corr, fontsize=10.5, va='bottom')

# ── Error arrows ───────────────────────────────────────────────────────────────
x_e1 = dt + 1.55   # predictor-error arrow x
x_e2 = dt + 2.50   # corrector-error arrow x

ax.annotate('', xy=(x_e1, N_end_true), xytext=(x_e1, N_pred),
            arrowprops=dict(arrowstyle='<->', color=c_pred, lw=1.8))
ax.text(x_e1 + 0.12, (N_end_true + N_pred) / 2,
        f'Predictor\nerror\n={abs(N_pred - N_end_true):.3f}',
        color=c_pred, fontsize=9.5, va='center', ha='left')

ax.annotate('', xy=(x_e2, N_end_true), xytext=(x_e2, N_corr),
            arrowprops=dict(arrowstyle='<->', color=c_corr, lw=1.8))
ax.text(x_e2 + 0.12, (N_end_true + N_corr) / 2,
        f'Corrector\nerror\n={abs(N_corr - N_end_true):.3f}',
        color=c_corr, fontsize=9.5, va='center', ha='left')

# Horizontal connector lines to arrows
ax.plot([dt, x_e1], [N_end_true, N_end_true], color='gray', lw=0.8, ls=':')
ax.plot([dt, x_e1], [N_pred, N_pred],         color=c_pred,  lw=0.8, ls=':')
ax.plot([dt, x_e2], [N_corr, N_corr],         color=c_corr,  lw=0.8, ls=':')

# ── Slope annotations ──────────────────────────────────────────────────────────
# Initial slope label – no arrow, placed just below the predictor line
ax.text(0.4, N0 + slope_0 * 0.4 - 0.17,
        f'Initial rate slope\n$= -k N_0 = {slope_0:.2f}$',
        fontsize=9.5, color=c_pred, ha='left', va='top')

# Midpoint slope label – shifted left so the slope arrow doesn't cross it
ax.annotate(
    f'Midpoint rate slope\n$= -k\\,N(\\Delta t/2) = {slope_mid:.3f}$',
    xy=(t_mid + 0.15, N_mid + slope_mid * 0.15),
    xytext=(1.5, 0.50),
    fontsize=9.5, color=c_corr,
    arrowprops=dict(arrowstyle='->', color=c_corr, lw=1.0, relpos=(1, 0)),
)

# Small slope-indicator segments
seg_len = 1.2
for x0, y0, slope, color in [
    (0.6, N0 + slope_0 * 0.6, slope_0, c_pred),
    (t_mid, N_mid,             slope_mid, c_corr),
]:
    ax.annotate('', xy=(x0 + seg_len, y0 + slope * seg_len),
                xytext=(x0, y0),
                arrowprops=dict(arrowstyle='->', color=color, lw=2.0))

# ── Formatting ─────────────────────────────────────────────────────────────────
ax.set_xlabel('Time  [EFPD]', fontsize=12)
ax.set_ylabel('Nuclide Number Density  (normalised)', fontsize=12)
ax.set_title(
    'Predictor–Corrector Integration of the Bateman Equations\n'
    r'(illustrative: $dN/dt = -kN$,  $k = 0.10\ \mathrm{EFPD}^{-1}$,  '
    r'$\Delta t = 8\ \mathrm{EFPD}$)',
    fontsize=11.5,
)
ax.legend(loc='upper right', fontsize=10, framealpha=0.92)
ax.set_xlim(-0.5, dt + 4.8)
ax.set_ylim(0.13, 1.12)
ax.grid(True, alpha=0.25)
ax.tick_params(labelsize=11)

plt.tight_layout()

out_path = 'scripts/predictor_corrector_diagram.png'
plt.savefig(out_path, dpi=150, bbox_inches='tight')
print(f'\nSaved → {out_path}')
