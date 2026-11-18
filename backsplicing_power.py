#!/usr/bin/env python3
"""Kinetic attenuation and detectability calculations for
"Does circadian timing modulate back-splicing? A testable hypothesis for circRNA regulation".

Standard library only (Python >= 3.8). Run:
    python3 backsplicing_power.py            # print every value reported in the manuscript
    python3 backsplicing_power.py --check    # assert that the printed values match the manuscript
    python3 backsplicing_power.py --design --mu 20 --phi 0.1 --timepoints 12 --reps 3 --half-life 24
                                             # evaluate your own design
"""
import argparse
import math

PERIOD_H = 24.0
OMEGA = 2 * math.pi / PERIOD_H          # angular frequency, rad h^-1
ALPHA = 0.05
POWER = 0.80


# ---------------- Section 6: first-order kinetics ----------------
def loss_rate(half_life_h, doubling_time_h=None):
    """Effective first-order loss constant k (h^-1): decay plus optional dilution by division."""
    k = math.log(2) / half_life_h
    if doubling_time_h:
        k += math.log(2) / doubling_time_h
    return k


def attenuation(half_life_h, doubling_time_h=None):
    """Steady-state amplitude attenuation H = 1 / sqrt(1 + (omega/k)^2) for sinusoidal synthesis."""
    k = loss_rate(half_life_h, doubling_time_h)
    return 1 / math.sqrt(1 + (OMEGA / k) ** 2)


def phase_lag_h(half_life_h, doubling_time_h=None):
    """Lag of the abundance rhythm behind synthesis, arctan(omega/k)/omega, in hours (always < 6 h)."""
    k = loss_rate(half_life_h, doubling_time_h)
    return math.atan(OMEGA / k) / OMEGA


# ---------------- Section 9.2: detectability ----------------
def coef_se(mu, phi, n_timepoints, reps):
    """Planning SE of each relative sine or cosine coefficient for negative-binomial counts:
    SE ~ sqrt(2 (1/mu + phi) / (n_timepoints * reps)), evenly sampled harmonic design."""
    return math.sqrt(2 * (1 / mu + phi) / (n_timepoints * reps))


def _chi2_df2_sf(x):
    return math.exp(-x / 2)


def _central_chi2_even_cdf(x, df):
    """CDF of a central chi-square with even df = 2m."""
    m = df // 2
    term, s = 1.0, 0.0
    for i in range(m):
        if i:
            term *= (x / 2) / i
        s += term
    return 1 - math.exp(-x / 2) * s


def _noncentral_chi2_df2_sf(x, lam, tol=1e-12):
    """P(X > x) for a noncentral chi-square with 2 df and noncentrality lam (Poisson mixture)."""
    total, j, w = 0.0, 0, math.exp(-lam / 2)
    while True:
        total += w * (1 - _central_chi2_even_cdf(x, 2 + 2 * j))
        j += 1
        w *= (lam / 2) / j
        if w < tol and j > lam:
            return total


def noncentrality(alpha=ALPHA, power=POWER):
    """Noncentrality giving the target power for the joint (2-df) sine-cosine Wald test."""
    crit = -2 * math.log(alpha)          # chi-square 2-df critical value
    lo, hi = 0.0, 100.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if _noncentral_chi2_df2_sf(crit, mid) < power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def detectable_amplitude(mu, phi, n_timepoints, reps, alpha=ALPHA, power=POWER):
    """Smallest steady-state fractional amplitude detectable at a single locus."""
    return coef_se(mu, phi, n_timepoints, reps) * math.sqrt(noncentrality(alpha, power))


def required_production_amplitude(a_obs_min, half_life_h, doubling_time_h=None):
    """Production amplitude needed for detection, A_prod,min = A_obs,min / H.
    Values > 1 are impossible for a non-negative sinusoidal production rate."""
    return a_obs_min / attenuation(half_life_h, doubling_time_h)


def amplitude_from_peak_trough(ratio):
    """Fractional amplitude A for a peak-to-trough ratio (1 + A) / (1 - A)."""
    return (ratio - 1) / (ratio + 1)


# ---------------- manuscript values ----------------
def manuscript_values():
    v = {}
    for t in (4, 12, 24, 48):
        v[f'H_{t}h'] = attenuation(t)
        v[f'lag_{t}h'] = phase_lag_h(t)
    v['H_48h_div24h'] = attenuation(48, doubling_time_h=24)
    v['k_dilution_24h'] = math.log(2) / 24
    v['lambda'] = noncentrality()
    for mu in (5, 20, 100):
        v[f'SE_mu{mu}'] = coef_se(mu, 0.1, 12, 3)
        v[f'Amin_mu{mu}'] = detectable_amplitude(mu, 0.1, 12, 3)
    for t in (4, 12, 24, 48):
        v[f'Aprod_0.24_{t}h'] = required_production_amplitude(0.24, t)
        v[f'Aprod_0.20_{t}h'] = required_production_amplitude(0.20, t)
    v['A_peak_trough_1.5'] = amplitude_from_peak_trough(1.5)
    return v


# value reported in the manuscript -> (key, decimals used when rounding)
MANUSCRIPT = {
    'H_4h': (0.552, 3), 'H_12h': (0.215, 3), 'H_48h': (0.055, 3),
    'lag_4h': (3.8, 1), 'lag_12h': (5.2, 1), 'lag_48h': (5.8, 1),
    'H_48h_div24h': (0.16, 2),
    'SE_mu5': (0.13, 2), 'SE_mu20': (0.09, 2), 'SE_mu100': (0.08, 2),
    'Amin_mu5': (0.40, 2), 'Amin_mu20': (0.28, 2), 'Amin_mu100': (0.24, 2),
    'Aprod_0.24_4h': (0.43, 2), 'Aprod_0.24_12h': (1.11, 2), 'Aprod_0.24_24h': (2.19, 2), 'Aprod_0.24_48h': (4.36, 2),
    'A_peak_trough_1.5': (0.2, 1),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true', help='assert agreement with the manuscript values')
    ap.add_argument('--design', action='store_true', help='evaluate a user-specified design')
    ap.add_argument('--mu', type=float, default=20, help='mean BSJ count per sample')
    ap.add_argument('--phi', type=float, default=0.1, help='negative-binomial dispersion')
    ap.add_argument('--timepoints', type=int, default=12)
    ap.add_argument('--reps', type=int, default=3, help='biological replicates per time point')
    ap.add_argument('--half-life', type=float, default=24, help='RNA half-life (h)')
    ap.add_argument('--doubling-time', type=float, default=None, help='cell doubling time (h), optional')
    ap.add_argument('--alpha', type=float, default=ALPHA)
    ap.add_argument('--power', type=float, default=POWER)
    a = ap.parse_args()

    if a.design:
        se = coef_se(a.mu, a.phi, a.timepoints, a.reps)
        amin = detectable_amplitude(a.mu, a.phi, a.timepoints, a.reps, a.alpha, a.power)
        h = attenuation(a.half_life, a.doubling_time)
        aprod = amin / h
        print(f'coefficient SE            {se:.3f}')
        print(f'detectable amplitude      {amin:.3f}  (alpha={a.alpha}, power={a.power})')
        print(f'attenuation H             {h:.3f}  (half-life {a.half_life} h)')
        print(f'phase lag                 {phase_lag_h(a.half_life, a.doubling_time):.1f} h')
        print(f'required production amp. {aprod:.2f}' + ('  -> exceeds 1: steady-state test uninformative' if aprod > 1 else ''))
        return

    v = manuscript_values()
    for k, x in v.items():
        print(f'{k:20s} {x:.4f}')
    if a.check:
        bad = []
        for k, (rep, d) in MANUSCRIPT.items():
            if round(v[k], d) != round(rep, d):
                bad.append(f'{k}: computed {v[k]:.4f}, manuscript {rep}')
        if bad:
            raise SystemExit('MISMATCH\n' + '\n'.join(bad))
        print(f'\nAll {len(MANUSCRIPT)} manuscript values reproduced.')


if __name__ == '__main__':
    main()
