# Orbital altitude and revisit time over Central Europe

Python analysis and paper figures for four 30-day GMAT cases (400, 500, 600 and 700 km).

## Contents

- `analyze_revisit.py`: visibility, observation windows and revisit metrics.
- `paper_figures.py`, `paper_style.py`: all ten figures and shared layout.
- `target_region.py`, `analysis_workflow.py`, `plot_maximum_gaps.py`: individual entry points.
- `config.json`, `requirements.txt`: study settings and Python dependencies.
- `data/`: four 30-day GMAT trajectory CSVs, target grid and boundary data.
- `results/`: per-target statistics, windows, daily coverage, summaries and analysis metadata.
- `figures/`: current PDF and PNG exports.
- `report/`: current LaTeX source and compiled paper.
- `archives/`: earlier ZIP packages, preserved for reference. Use the root-level code for the current figure layout.
- `docs/`: provenance and reproduction notes.

## Run on macOS

In Terminal, open this extracted project directory, then run:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 paper_figures.py
```

This repository includes the four trajectory CSVs in `data/trajectories/`, saved result tables, the target grid and cached Natural Earth boundaries. All figures can be redrawn from these files. The ground-track figure reads the trajectories; other result plots use the saved tables.

To redraw only the maps or workflow without trajectories:

```bash
python3 paper_figures.py --group maximum
python3 paper_figures.py --group target
python3 paper_figures.py --group workflow
```

To recalculate the analysis from trajectories, run `python3 analyze_revisit.py`. This takes longer than redrawing. It writes result tables, checks and metadata before generating figures. `python3 target_region.py` recreates the target grid from `config.json`.

## Figure layout

`paper_style.py` contains shared typography, line widths and export settings. `paper_figures.py` draws all ten figures. The original entry points `analysis_workflow.py` and `plot_maximum_gaps.py` use that shared layout.

Exports are in `figures/`: PDF for Overleaf, PNG for previews. Labels use serif type; curves use black line patterns; map panels share a sequential blue scale; the workflow uses simple rectangles. Narrative explanations belong in the report captions.

The layout revision does not change the numerical analysis or any saved result table. `results/analysis_metadata.json` records the original configuration, input hashes and numerical checks. The plotting script rejects a configuration that differs from those saved results.

## Overleaf

Upload `report/Orbit_Altitude_Report_Updated.tex` and the `figures/` folder into your existing AIAA project. Keep the `new-aiaa.cls` and support files already in your project. Set the uploaded TEX file as the main document. The separate Overleaf ZIP arranges these files at its top level.

Affiliation: Space Systems Engineer, AiShed Consulting.

## Scientific scope

The region is a user-defined rectangle (5–25° E, 45–55° N), including land and water, sampled at 800 cell centres. The ideal nadir sensor has a 20° half-angle. Access uses WGS 84 surface targets, a sensor-cone test and an ellipsoid-horizon test. Results exclude cloud cover, illumination restrictions, scheduling and downlink constraints.

Ten-second trajectory samples and a 0.5° target grid give preliminary geometric-access estimates. The 20 s subsample is a coarsening diagnostic, not evidence of convergence. Finer sampling, grid and orientation sensitivity remain outstanding. Native GMAT scripts were not supplied, so force-model and integrator settings have not been audited.

Map boundary data: Natural Earth Admin 0, 1:50m (public domain), https://www.naturalearthdata.com/.

## Authorship and reuse

Prepared for Danillo Krüppel's engineering portfolio. Report affiliation: Space Systems Engineer, AiShed Consulting. An affiliation does not imply company review or endorsement. No software reuse licence has been selected for this repository; contact the author before redistribution. Natural Earth boundaries have their own public-domain status as noted above.
