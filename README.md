# circadian-backsplicing-power

Calculations supporting the manuscript **"Does circadian timing modulate back-splicing? A testable hypothesis for circRNA regulation"** (Li, Režen, Gao, Smulders-Srinivasan and Uliassi).

The script reproduces every kinetic and detectability value in the manuscript, and readers can rerun it for their own designs. It needs only the Python standard library (Python ≥ 3.8). There are no packages to install.

## Usage

```bash
python3 backsplicing_power.py --check     # print all values and verify them against the manuscript
python3 backsplicing_power.py --design --mu 20 --phi 0.1 --timepoints 12 --reps 3 --half-life 24
```

`--design` prints the following for a single locus:
- the coefficient standard error;
- the smallest detectable steady-state amplitude;
- the kinetic attenuation and phase lag;
- the production amplitude needed for detection. Values above 1 are impossible, so a steady-state test is uninformative.

Use `--doubling-time` (hours) to add dilution by cell division.

## Equations

**First-order kinetics with sinusoidal synthesis**
- Angular frequency: ω = 2π/24 h⁻¹.
- Loss rate: k = ln2/t½, plus ln2/T_d (the doubling time) in dividing cells.
- Amplitude attenuation: H = 1/√(1 + (ω/k)²).
- Phase lag of abundance behind synthesis: arctan(ω/k)/ω. This is always less than 6 h, a quarter cycle.

**Planning approximation for negative-binomial junction counts**
- Inputs: mean count μ, dispersion φ, n_T time points, r replicates.
- Standard error of each relative sine or cosine coefficient: SE ≈ √(2(1/μ + φ)/(n_T·r)).
- Smallest detectable steady-state amplitude, for a joint 2-df sine–cosine test at α = 0.05 with 80% power: A_min = SE·√λ, where λ = 9.63.
- Production amplitude required for detection: A_prod,min = A_min/H.
- Peak-to-trough ratio R corresponds to fractional amplitude A = (R − 1)/(R + 1). For example, R = 1.5 gives A = 0.2.

These are planning heuristics for an ideal balanced design, not a full power calculation. For confirmatory work, use simulations tailored to the exact sampling schedule, waveform and rhythm-detection method.

## Values reproduced (manuscript Sections 6 and 9.2)

| Quantity | Value |
|---|---|
| H at half-lives of 4, 12, 24, 48 h | 0.552, 0.215, 0.110, 0.055 |
| Phase lag at 4, 12, 48 h | 3.8, 5.2, 5.8 h |
| H at a 48-h half-life with a 24-h doubling time | 0.16 |
| SE at μ = 5, 20, 100 (12 time points × 3 replicates, φ = 0.1) | 0.13, 0.09, 0.08 |
| Detectable amplitude at the same depths | 0.40, 0.28, 0.24 |
| Required production amplitude at A_min = 0.24, for 4, 12, 24, 48 h | 0.43, 1.11, 2.19, 4.36 |

## Citation

See `CITATION.cff`. Please cite the archived release (Zenodo DOI) and the manuscript.

## Licence

MIT (see `LICENSE`).
