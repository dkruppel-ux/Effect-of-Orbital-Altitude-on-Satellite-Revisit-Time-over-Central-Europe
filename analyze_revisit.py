"""Analyse GMAT EarthFixed trajectories for a nadir-pointing conical sensor.
Run: python analyze_revisit.py
See README.txt for input placement, assumptions, outputs and limitations.
"""
import os
# Bound small-matrix BLAS parallelism for predictable memory use.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

BASE = Path(__file__).resolve().parent
DATA, RESULTS, FIGURES = BASE / "data", BASE / "results", BASE / "figures"
A = 6378.137  # WGS 84 semi-major axis in km
F = 1 / 298.257223563
E2 = F * (2-F)


def targets_ecef(lon_deg, lat_deg):
    """Zero-height WGS 84 geodetic targets and outward unit surface normals."""
    lon, lat = np.deg2rad(lon_deg), np.deg2rad(lat_deg)
    cl, sl, cp, sp = np.cos(lon), np.sin(lon), np.cos(lat), np.sin(lat)
    n = A / np.sqrt(1-E2*sp*sp)
    ground = np.column_stack((n*cp*cl, n*cp*sl, n*(1-E2)*sp))
    normal = np.column_stack((cp*cl, cp*sl, sp))
    return ground, normal


def visibility(r, ground, normal, half_angle_deg):
    """Conical-FOV and convex ellipsoid surface-horizon test.

    r: spacecraft EarthFixed position (km), shape (time,3).
    Boresight is -r/|r|, pointing at the Earth centre (geocentric nadir).
    A surface target is unobstructed when (r-ground).normal > 0.
    """
    r2 = np.sum(r*r, axis=1)[:, None]
    dot = r @ ground.T
    d2 = np.maximum(r2 + np.sum(ground*ground,axis=1)[None,:] - 2*dot, 1e-18)
    cos_off_nadir = (r2-dot) / np.sqrt(r2*d2)
    inside = cos_off_nadir >= np.cos(np.deg2rad(half_angle_deg))
    horizon = r @ normal.T - np.sum(ground*normal,axis=1)[None,:]
    return inside & (horizon > 0)


def geometry_checks(half_angle_deg):
    """Independent direct-vector checks plus analytic equatorial examples."""
    g,n=targets_ecef(np.array([0.,2.,4.,180.]),np.zeros(4))
    r=np.array([[A+500,0.,0.]])
    fast=visibility(r,g,n,half_angle_deg)[0]
    direct=[]
    for p,q in zip(g,n):
        v=p-r[0]
        c=np.dot(v,-r[0])/(np.linalg.norm(v)*np.linalg.norm(r[0]))
        direct.append(c>=np.cos(np.deg2rad(half_angle_deg)) and np.dot(r[0]-p,q)>0)
    if not np.array_equal(fast,direct):
        raise ValueError("Vectorised and direct geometry disagree.")
    if not fast[0] or fast[-1]:
        raise ValueError("Nadir/antipodal visibility check failed.")
    # On the equator, the ellipsoid section is a circle of radius A.
    alpha=np.deg2rad(half_angle_deg)
    limb=np.arcsin(A/(A+500))
    if alpha<limb:
        psi=np.arcsin((A+500)/A*np.sin(alpha))-alpha
        gg,nn=targets_ecef(np.rad2deg(np.array([psi*.99,psi*1.01])),np.zeros(2))
        if not np.array_equal(visibility(r,gg,nn,half_angle_deg)[0],[True,False]):
            raise ValueError("Analytic footprint-boundary check failed.")
    return {"direct_vector_comparison":"passed", "nadir_and_antipodal":"passed",
            "equatorial_analytic_boundary":"passed"}


def read_trajectory(path, expected_altitude, expected_days):
    with path.open(newline="") as f:
        reader=csv.reader(f); header=next(reader); first=next(reader)
    elapsed=next(k for k in header if k.endswith(".ElapsedSecs"))
    prefix=elapsed.rsplit(".",1)[0]
    wanted=[prefix+".A1ModJulian",elapsed]+[prefix+".EarthFixed."+c for c in "XYZ"]
    missing=[k for k in wanted if k not in header]
    if missing:
        raise ValueError(f"{path.name}: missing {missing}; export spacecraft EarthFixed coordinates.")
    arr=np.loadtxt(path,delimiter=",",skiprows=1,usecols=[header.index(k) for k in wanted],ndmin=2)
    t,r=arr[:,1],arr[:,2:]
    if len(t)<2 or not np.isfinite(arr).all() or np.any(np.diff(t)<=0):
        raise ValueError(f"{path.name}: nonfinite data, duplicate/reversed times, or too few rows.")
    duration=expected_days*86400
    if abs(t[0])>1e-3 or abs(t[-1]-duration)>1e-3:
        raise ValueError(f"{path.name}: ends at {t[-1]:g} s; expected {duration:g} s.")
    initial_radius=np.linalg.norm(r[0])
    if abs(initial_radius-(A+expected_altitude))>0.1:
        raise ValueError(f"{path.name}: initial radius does not match labelled altitude.")
    # Numerical consistency between elapsed seconds and absolute A1 epochs.
    if np.max(np.abs((arr[:,0]-arr[0,0])*86400-t))>0.02:
        raise ValueError(f"{path.name}: elapsed and absolute timestamps disagree.")
    check={"rows":len(t),"duration_days":float(t[-1]/86400),
           "initial_a1_modjulian":float(arr[0,0]),"initial_radius_km":float(initial_radius),
           "median_step_s":float(np.median(np.diff(t))),"max_step_s":float(np.diff(t).max()),
           "first_reported_date":first[header.index(prefix+".A1Gregorian")] if prefix+".A1Gregorian" in header else None}
    return t,r,arr[:,0],check


def observation_windows(t, visible):
    """Midpoint boundary estimate; never claims sub-sample event accuracy.

    A detection at sample k starts halfway between k-1 and k and ends halfway
    between the last visible and next invisible sample. Boundary windows
    are clipped to the available time interval, with flags in the CSV.
    """
    ids=np.flatnonzero(visible)
    if not len(ids):
        return np.empty(0),np.empty(0),np.empty(0,dtype=int),np.empty(0,dtype=int)
    breaks=np.flatnonzero(np.diff(ids)>1)
    first=ids[np.r_[0,breaks+1]]
    last=ids[np.r_[breaks,len(ids)-1]]
    starts=t[first].copy();ends=t[last].copy()
    interior=first>0;starts[interior]=(t[first[interior]-1]+t[first[interior]])/2
    interior=last<len(t)-1;ends[interior]=(t[last[interior]]+t[last[interior]+1])/2
    return starts,ends,first,last


def weighted_mean_available(values,weights):
    valid=np.isfinite(values)
    return float(np.sum(values[valid]*weights[valid])/weights[valid].sum()) if valid.any() else np.nan


def save_csv(path, fieldnames, rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames);w.writeheader();w.writerows(rows)


def analyse_case(t,mask,targets,days,altitude,write_windows=True):
    weights=targets["area_weight"]
    daily=np.zeros((days,len(targets)),dtype=bool)
    per=[];all_windows=[]
    for j,target in enumerate(targets):
        starts,ends,first,last=observation_windows(t,mask[:,j])
        gaps=(starts[1:]-ends[:-1])/3600
        mean=float(gaps.mean()) if len(gaps) else np.nan
        maximum=float(gaps.max()) if len(gaps) else np.nan
        per.append({"target_id":int(target["target_id"]),"longitude_deg":float(target["longitude_deg"]),
                    "latitude_deg":float(target["latitude_deg"]),"area_weight":float(weights[j]),
                    "window_count":len(starts),"complete_gap_count":len(gaps),"mean_gap_h":mean,
                    "max_gap_h":maximum,
                    "leading_unobserved_h":float((starts[0]-t[0])/3600) if len(starts) else float((t[-1]-t[0])/3600),
                    "trailing_unobserved_h":float((t[-1]-ends[-1])/3600) if len(starts) else float((t[-1]-t[0])/3600),
                    "no_observations":not bool(len(starts))})
        for k,(s,e) in enumerate(zip(starts,ends)):
            if e>s:
                day0=max(0,int(np.floor(s/86400)))
                day1=min(days-1,int(np.ceil(e/86400))-1)
                if day1>=day0:daily[day0:day1+1,j]=True
            if write_windows:
                all_windows.append({"target_id":int(target["target_id"]),"window_index":k+1,
                                    "start_elapsed_s":float(s),"end_elapsed_s":float(e),
                                    "duration_s":float(e-s),"start_clipped":bool(first[k]==0),
                                    "end_clipped":bool(last[k]==len(t)-1)})
    means=np.array([p["mean_gap_h"] for p in per]);maxes=np.array([p["max_gap_h"] for p in per])
    coverage=daily @ weights * 100
    summary={"altitude_km":altitude,"area_weighted_mean_target_gap_h":weighted_mean_available(means,weights),
             "worst_observed_complete_gap_h":float(np.nanmax(maxes)) if np.isfinite(maxes).any() else np.nan,
             "mean_daily_coverage_pct":float(coverage.mean()),"min_daily_coverage_pct":float(coverage.min()),
             "targets_without_observations":sum(p["no_observations"] for p in per),
             "targets_without_complete_gaps":int((~np.isfinite(means)).sum()),
             "area_without_complete_gaps_pct":float(weights[~np.isfinite(means)].sum()*100),
             "observation_window_count":sum(p["window_count"] for p in per)}
    if write_windows:
        save_csv(RESULTS/f"targets_{altitude}km.csv",list(per[0]),per)
        winfields=["target_id","window_index","start_elapsed_s","end_elapsed_s","duration_s","start_clipped","end_clipped"]
        save_csv(RESULTS/f"windows_{altitude}km.csv",winfields,all_windows)
        save_csv(RESULTS/f"daily_coverage_{altitude}km.csv",["elapsed_day","coverage_pct"],
                 [{"elapsed_day":d+1,"coverage_pct":float(c)} for d,c in enumerate(coverage)])
    return summary,means,maxes,coverage,per


def save_figure(fig,name):
    fig.savefig(FIGURES/(name+".pdf"),bbox_inches="tight")
    fig.savefig(FIGURES/(name+".png"),dpi=180,bbox_inches="tight")
    plt.close(fig)


def make_plots(config,targets,cases,timelines,probe,groundtracks):
    # Presentation only: reuse saved tables and the shared paper figure layout.
    from paper_figures import setup, load_targets, analysis_plots, maps
    plot_config, nearby = setup()
    records, lon, lat = load_targets(plot_config)
    analysis_plots(plot_config, nearby, records, lon, lat)
    maps(plot_config, nearby, records, lon, lat, "max_gap_h",
         "maximum_gap_maps", "Maximum complete gap [h]", mark_max=True)


def main():
    begun=time.time()
    config=json.loads((BASE/"config.json").read_text())
    half=config["sensor_half_angle_deg"]
    if not 0<half<90:raise ValueError("Sensor half-angle must lie between 0 and 90 degrees.")
    for p in [RESULTS,FIGURES]:p.mkdir(exist_ok=True)
    paths={alt:DATA/"trajectories"/f"trajectory_{alt}km_30days.csv" for alt in config["altitudes_km"]}
    absent=[p.name for p in paths.values() if not p.exists()]
    if absent:raise FileNotFoundError("Copy these Drive CSVs into data/trajectories: "+", ".join(absent))
    targets=np.genfromtxt(DATA/"target_grid.csv",delimiter=",",names=True)
    if not np.isclose(targets["area_weight"].sum(),1):raise ValueError("Grid weights must sum to one.")
    expected=np.array([(x,y) for y in np.unique(targets["latitude_deg"]) for x in np.unique(targets["longitude_deg"])])
    if expected.shape!=(len(targets),2) or not np.allclose(expected,np.column_stack((targets["longitude_deg"],targets["latitude_deg"]))):
        raise ValueError("Map plotting expects the supplied latitude-major rectangular grid.")
    ground,normal=targets_ecef(targets["longitude_deg"],targets["latitude_deg"])
    checks=geometry_checks(half)
    cases={};probe={};validation={};timelines={};groundtracks={};reference_t=None;reference_epoch=None
    example=int(np.argmin((targets["longitude_deg"]-15)**2+(targets["latitude_deg"]-50)**2))
    for alt,path in paths.items():
        print(f"{alt} km: reading trajectory...",flush=True)
        t,r,epoch,validation[alt]=read_trajectory(path,alt,config["expected_duration_days"])
        if reference_t is None:reference_t=t.copy();reference_epoch=float(epoch[0])
        elif len(t)!=len(reference_t) or not np.allclose(t,reference_t,rtol=0,atol=1e-4) or abs(float(epoch[0])-reference_epoch)>1e-9:
            raise ValueError("The altitude cases must have matching epochs and timestamps.")
        mask=np.empty((len(t),len(targets)),dtype=bool)
        chunk=config["chunk_rows"]
        for start in range(0,len(t),chunk):
            end=min(len(t),start+chunk)
            mask[start:end]=visibility(r[start:end],ground,normal,half)
        cases[alt]=analyse_case(t,mask,targets,config["expected_duration_days"],alt)
        probe[alt]=analyse_case(t[::2],mask[::2],targets,config["expected_duration_days"],alt,False)[0]
        s,e,_,_=observation_windows(t,mask[:,example]);timelines[alt]=(s,e)
        first_day=t<=86400+1e-3;rr=r[first_day];rad=np.linalg.norm(rr,axis=1)
        groundtracks[alt]=(np.rad2deg(np.arctan2(rr[:,1],rr[:,0])),np.rad2deg(np.arcsin(np.clip(rr[:,2]/rad,-1,1))))
        print(f"{alt} km: mean gap {cases[alt][0]['area_weighted_mean_target_gap_h']:.3f} h; "
              f"daily coverage {cases[alt][0]['mean_daily_coverage_pct']:.1f}%",flush=True)
        del mask,r,epoch
    rows=[cases[a][0] for a in cases]
    save_csv(RESULTS/"comparison_summary.csv",list(rows[0]),rows)
    probe_rows=[probe[a] for a in probe]
    save_csv(RESULTS/"coarse_20s_summary.csv",list(probe_rows[0]),probe_rows)
    report={"configuration":config,"input_validation":validation,"geometry_checks":checks,
            "target_count":len(targets),"example_target":{k:float(targets[k][example]) for k in ["longitude_deg","latitude_deg"]},
            "geometry":"WGS84 zero-height geodetic targets, EarthFixed positions, geocentric nadir cone, ellipsoid horizon",
            "window_boundaries":"midpoints of visible/invisible sample brackets; endpoint windows clipped",
            "gap_definition":"end of one window to start of next; boundary unobserved intervals reported separately",
            "sampling_diagnostic":"20 s subsampling is coarsening only; finer GMAT output required for convergence",
            "scope":"geometric access only; no clouds, daylight requirement, scheduling, downlink or image-quality model",
            "input_sha256":{str(a):hashlib.sha256(p.read_bytes()).hexdigest() for a,p in paths.items()},
            "python_version":sys.version,"numpy_version":np.__version__,"matplotlib_version":matplotlib.__version__,
            "platform":platform.platform()}
    (RESULTS/"analysis_metadata.json").write_text(json.dumps(report,indent=2)+"\n")
    print("Generating figures...",flush=True)
    make_plots(config,targets,cases,timelines,probe,groundtracks)
    print(f"Finished in {time.time()-begun:.1f} s. Results: {RESULTS}\nFigures: {FIGURES}",flush=True)


if __name__=="__main__":
    main()
