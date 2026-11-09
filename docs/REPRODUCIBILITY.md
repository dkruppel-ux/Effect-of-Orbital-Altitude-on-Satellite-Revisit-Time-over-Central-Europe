# Reproduction and provenance

The four CSVs were exported by the user from GMAT and cover 30 days starting 2026-01-01 UTC at approximately 10 s intervals. `trajectory_manifest.json` records their SHA-256 hashes and byte counts. The original analysis metadata also records input hashes and checks.

The CSV headers identify Earth-fixed positions in kilometres, elapsed seconds and A1 epoch fields. Native GMAT `.script` files were not supplied. Their force model, integrator and complete initial-state configuration therefore cannot be reproduced from this repository alone. No GMAT configuration is reconstructed or invented here.

The report and root-level code contain the current figure layout. Files in `archives/` preserve earlier source packages; some generate the previous figure style. Use the root-level entry points for new work.

The 2026-11-09 layout revision retained every saved numerical result table byte-for-byte and retained the numerical analysis functions. Only the plotting implementation and matching captions were changed. All ten figures were regenerated, Python files compiled, and the report compiled and visually inspected.

Recalculating the analysis overwrites `results/` and the analysis figures. Redrawing with `paper_figures.py` only reads saved results and inputs and overwrites `figures/`. Commit a checkpoint before changing settings or rerunning.

Natural Earth country boundaries are cached for offline plotting; source: https://www.naturalearthdata.com/ (Admin 0, 1:50m). The grid is the full study rectangle and is not clipped to national borders or land.
