"""Redraw study figures from existing result tables without recalculating access.

Run: python paper_figures.py
Individual groups: --group target, workflow, maximum, or analysis.
No numerical results or raw trajectory files are modified.
"""
from pathlib import Path
import argparse
import csv
import json
import string
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch
import matplotlib.patheffects as pe
from paper_style import apply_style, label_panel, light_grid, save, LINE_STYLES

BASE = Path(__file__).resolve().parent
OUT, DATA, RESULTS = BASE / "figures", BASE / "data", BASE / "results"


def setup():
    apply_style()
    config = json.loads((BASE / "config.json").read_text())
    west, south, east, north = config["region_bounds_lon_lat"]
    countries = gpd.read_file(DATA / "ne_50m_admin_0_countries.zip").to_crs("EPSG:4326")
    nearby = countries.cx[west-2:east+2, south-2:north+2]
    return config, nearby


def load_targets(config):
    metadata = json.loads((RESULTS / "analysis_metadata.json").read_text())
    if metadata["configuration"] != config:
        raise ValueError("Configuration differs from saved results. Run analyze_revisit.py first.")
    alts = config["altitudes_km"]
    records = {a: np.genfromtxt(RESULTS / f"targets_{a}km.csv", delimiter=",", names=True) for a in alts}
    first = records[alts[0]]
    lon, lat = np.unique(first["longitude_deg"]), np.unique(first["latitude_deg"])
    expected = np.array([(x,y) for y in lat for x in lon])
    if not np.array_equal(expected, np.column_stack((first["longitude_deg"], first["latitude_deg"]))):
        raise ValueError("Expected a latitude-major rectangular grid.")
    for a in alts:
        for key in ["target_id", "longitude_deg", "latitude_deg"]:
            if not np.array_equal(records[a][key], first[key]):
                raise ValueError("Altitude cases use different target grids.")
    return records, lon, lat


def map_axes(ax, config, geocentric=False):
    w,s,e,n = config["region_bounds_lon_lat"]
    ax.set(xlim=(w,e), ylim=(s,n), xticks=np.linspace(w,e,5), yticks=np.linspace(s,n,3),
           xlabel="Longitude [°E]", ylabel=("Geocentric latitude [°N]" if geocentric else "Latitude [°N]"))
    ax.set_aspect(1/np.cos(np.deg2rad((s+n)/2)))


def target_map(config, nearby):
    grid = np.genfromtxt(DATA/"target_grid.csv", delimiter=",", names=True)
    w,s,e,n = config["region_bounds_lon_lat"]
    fig,ax = plt.subplots(figsize=(6.8,4.4),layout="constrained")
    nearby.plot(ax=ax,facecolor="0.97",edgecolor="0.55",linewidth=0.45)
    ax.scatter(grid["longitude_deg"],grid["latitude_deg"],s=1.8,c="0.50",linewidths=0,zorder=3)
    ax.add_patch(Rectangle((w,s),e-w,n-s,fill=False,edgecolor="0.05",linestyle="--",linewidth=.8,zorder=4))
    labels={"Germany":(10.5,51),"Poland":(19,52),"Czechia":(15.4,49.8),"Austria":(13,47.5),
            "Switzerland":(8.1,46.7),"Slovakia":(19.3,48.7),"Hungary":(19.4,47.1),"Slovenia":(14.7,46.1)}
    for name,(x,y) in labels.items():
        if w<x<e and s<y<n:
            ax.text(x,y,name,ha="center",va="center",fontsize=8,zorder=5,
                    path_effects=[pe.withStroke(linewidth=2,foreground="white")])
    ax.set(xlim=(w-1,e+1),ylim=(s-.7,n+.7),xlabel="Longitude [°E]",ylabel="Latitude [°N]",
           xticks=np.linspace(w,e,5),yticks=np.linspace(s,n,3))
    ax.set_aspect(1/np.cos(np.deg2rad((s+n)/2)))
    ax.legend(handles=[Line2D([],[],color="0.1",linestyle="--",linewidth=.8,label="Study boundary"),
                       Line2D([],[],color="0.5",marker=".",linestyle="none",markersize=3,label="Target centres")],
              loc="upper center",bbox_to_anchor=(.5,-.12),ncol=2,handlelength=2)
    save(fig,OUT,"target_region")


def workflow(config):
    fig,ax = plt.subplots(figsize=(6.8,5.1))
    ax.set_xlim(0,10);ax.set_ylim(0,8.2);ax.axis("off")
    def box(x,y,w,h,text):
        ax.add_patch(Rectangle((x-w/2,y-h/2),w,h,facecolor="white",edgecolor="0.15",linewidth=.65))
        ax.text(x,y,text,ha="center",va="center",fontsize=9,linespacing=1.4)
    def arrow(a,b,vertices=None):
        if vertices:
            xx,yy=zip(*vertices);ax.plot(xx,yy,color="0.15",linewidth=.65);a=vertices[-1]
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle="-|>",mutation_scale=8,color="0.15",linewidth=.65))
    alt_text=", ".join(str(a) for a in config["altitudes_km"])
    box(2.5,7.5,4.4,.85,"Fixed inputs\nTarget grid, sensor geometry, epoch\nand orbit orientation")
    box(7.5,7.5,4.4,.85,f"Altitude cases\n{alt_text} km")
    arrow((2.5,7.04),(4,6.53));arrow((7.5,7.04),(6,6.53))
    box(5,6.05,6.5,.85,"GMAT trajectory output\nEarth-fixed spacecraft position versus time")
    arrow((5,5.59),(5,5.15))
    box(5,4.68,6.5,.85,"Target visibility\nSensor-cone test and Earth-occultation test")
    arrow((5,4.22),(5,3.78))
    box(5,3.25,6.5,.95,"Observation windows\n"+r"$\tau_{j,k}=t^{\mathrm{start}}_{j,k+1}-t^{\mathrm{end}}_{j,k}$")
    ax.plot([5,5],[2.74,2.28],color="0.15",linewidth=.65)
    ax.plot([1.8,8.2],[2.28,2.28],color="0.15",linewidth=.65)
    for x in [1.8,5,8.2]:arrow((x,2.28),(x,1.83))
    box(1.8,1.32,2.8,.95,"Mean revisit gap\n"+r"$\overline{\tau}_j$")
    box(5,1.32,2.8,.95,"Maximum revisit gap\n"+r"$\tau_j^{\max}$")
    box(8.2,1.32,2.8,.95,"Daily coverage\n"+r"$C_d$")
    fig.subplots_adjust(left=.015,right=.985,top=.99,bottom=.01)
    save(fig,OUT,"analysis_workflow")


def maps(config,nearby,records,lon,lat,key,name,colorlabel,mark_max=False):
    alts=config["altitudes_km"]
    vals=np.concatenate([records[a][key] for a in alts]);finite=vals[np.isfinite(vals)]
    if not len(finite):raise ValueError("No finite data for map.")
    vmax=max(1,float(np.ceil(finite.max()/10)*10))
    norm=Normalize(0,vmax)
    cmap=plt.get_cmap("Blues").copy();cmap.set_bad("0.75")
    step=config["grid_step_deg"]
    xx=np.r_[lon-step/2,lon[-1]+step/2];yy=np.r_[lat-step/2,lat[-1]+step/2]
    fig,axes=plt.subplots(2,2,figsize=(7.0,5.2),layout="constrained")
    for i,(ax,a) in enumerate(zip(axes.flat,alts)):
        d=records[a];z=d[key].reshape(len(lat),len(lon))
        mesh=ax.pcolormesh(xx,yy,np.ma.masked_invalid(z),cmap=cmap,norm=norm,shading="flat",rasterized=True)
        nearby.boundary.plot(ax=ax,color="0.4",linewidth=.35,alpha=.9)
        if mark_max and np.isfinite(d[key]).any():
            j=int(np.nanargmax(d[key]));p=ax.scatter(d["longitude_deg"][j],d["latitude_deg"][j],
                    s=22,facecolors="none",edgecolors="white",linewidths=.7,zorder=5)
            p.set_path_effects([pe.Stroke(linewidth=1.5,foreground="0.15"),pe.Normal()])
        map_axes(ax,config);label_panel(ax,f"({string.ascii_lowercase[i]}) {a} km")
    cb=fig.colorbar(mesh,ax=list(axes.flat),shrink=.84,pad=.02,fraction=.035)
    cb.set_label(colorlabel,fontsize=9);cb.ax.tick_params(labelsize=8,width=.6,length=3)
    save(fig,OUT,name)


def read_summary(name):
    with (RESULTS/name).open(newline="") as f:
        return {int(r["altitude_km"]):{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)}


def analysis_plots(config,nearby,records,lon,lat):
    alts=config["altitudes_km"];summary=read_summary("comparison_summary.csv");coarse=read_summary("coarse_20s_summary.csv")
    fig,axes=plt.subplots(1,2,figsize=(6.9,2.6),layout="constrained")
    for ax,key,label,panel in [(axes[0],"area_weighted_mean_target_gap_h","Regional mean gap [h]","(a)"),
                               (axes[1],"worst_observed_complete_gap_h","Largest complete gap [h]","(b)")]:
        ax.plot(alts,[summary[a][key] for a in alts],color="0.15",marker="o",markersize=3)
        ax.set(xlabel="Initial altitude [km]",ylabel=label,xticks=alts);label_panel(ax,panel);light_grid(ax)
    save(fig,OUT,"altitude_revisit")

    fig,ax=plt.subplots(figsize=(6.9,2.9),layout="constrained")
    for i,a in enumerate(alts):
        d=np.genfromtxt(RESULTS/f"daily_coverage_{a}km.csv",delimiter=",",names=True)
        ax.plot(d["elapsed_day"],d["coverage_pct"],color="0.15",linestyle=LINE_STYLES[i],label=f"{a} km")
    ax.set(xlabel="Elapsed day",ylabel="Daily region coverage [%]",xlim=(1,30),ylim=(0,100),xticks=[1,5,10,15,20,25,30])
    light_grid(ax);ax.legend(loc="lower left",ncol=4,handlelength=3)
    save(fig,OUT,"daily_coverage")

    maps(config,nearby,records,lon,lat,"mean_gap_h","revisit_maps","Mean complete gap [h]")
    fig,ax=plt.subplots(figsize=(6.9,2.9),layout="constrained")
    for i,a in enumerate(alts):
        d=records[a];valid=np.isfinite(d["mean_gap_h"]);idx=np.argsort(d["mean_gap_h"][valid])
        ax.step(d["mean_gap_h"][valid][idx],np.cumsum(d["area_weight"][valid][idx])*100,where="post",
                color="0.15",linestyle=LINE_STYLES[i],label=f"{a} km")
    ax.set(xlabel="Per-target mean complete gap [h]",ylabel="Cumulative region fraction [%]",ylim=(0,100))
    light_grid(ax);ax.legend(loc="lower right",ncol=2,handlelength=3)
    save(fig,OUT,"revisit_distribution")

    metadata=json.loads((RESULTS/"analysis_metadata.json").read_text())
    ex=metadata["example_target"];first=records[alts[0]]
    j=int(np.argmin((first["longitude_deg"]-ex["longitude_deg"])**2+(first["latitude_deg"]-ex["latitude_deg"])**2))
    target_id=int(first["target_id"][j])
    fig,ax=plt.subplots(figsize=(6.9,2.5),layout="constrained")
    for i,a in enumerate(alts):
        windows=np.genfromtxt(RESULTS/f"windows_{a}km.csv",delimiter=",",names=True)
        keep=(windows["target_id"]==target_id)&(windows["start_elapsed_s"]<72*3600)&(windows["end_elapsed_s"]>0)
        ss,ee=windows["start_elapsed_s"][keep],windows["end_elapsed_s"][keep]
        ax.hlines(i,0,72,color="0.85",linewidth=.4)
        ax.broken_barh([(s/3600,(min(e,72*3600)-s)/3600) for s,e in zip(ss,ee)],(i-.18,.36),facecolors="0.15")
        ax.scatter((ss+ee)/7200,np.full(len(ss),i),marker="|",s=50,color="0.1",linewidths=.9)
    ax.set(xlabel="Elapsed time [h]",yticks=range(4),yticklabels=[f"{a} km" for a in alts],xlim=(0,72),xticks=np.arange(0,73,12),ylim=(-.5,3.5))
    save(fig,OUT,"observation_timeline")

    fig,axes=plt.subplots(1,2,figsize=(6.9,2.7),layout="constrained")
    for ax,key,label,panel in [(axes[0],"area_weighted_mean_target_gap_h","Regional mean gap [h]","(a)"),
                               (axes[1],"worst_observed_complete_gap_h","Largest complete gap [h]","(b)")]:
        ax.plot(alts,[summary[a][key] for a in alts],color="0.15",marker="o",markersize=3,label="10 s")
        ax.plot(alts,[coarse[a][key] for a in alts],color="0.15",linestyle="--",marker="s",markersize=3,markerfacecolor="white",label="20 s")
        ax.set(xlabel="Initial altitude [km]",ylabel=label,xticks=alts);label_panel(ax,panel);light_grid(ax)
        ax.legend(loc="upper right",handlelength=2.5)
    save(fig,OUT,"sampling_diagnostic")

    # Original first-day plotted samples are reused (same geocentric latitude).
    fig,axes=plt.subplots(2,2,figsize=(7.0,5.2),layout="constrained")
    for i,(ax,a) in enumerate(zip(axes.flat,alts)):
        p=DATA/"trajectories"/f"trajectory_{a}km_30days.csv"
        with p.open() as f:header=next(csv.reader(f))
        tk=next(k for k in header if k.endswith(".ElapsedSecs"));prefix=tk.rsplit(".",1)[0]
        cols=[header.index(k) for k in [tk]+[prefix+".EarthFixed."+c for c in "XYZ"]]
        arr=np.loadtxt(p,delimiter=",",skiprows=1,usecols=cols)
        r=arr[arr[:,0]<=86400+1e-3,1:];radius=np.linalg.norm(r,axis=1)
        glon=np.rad2deg(np.arctan2(r[:,1],r[:,0]));glat=np.rad2deg(np.arcsin(np.clip(r[:,2]/radius,-1,1)))
        jump=np.r_[False,np.abs(np.diff(glon))>180];glon[jump]=np.nan;glat[jump]=np.nan
        nearby.plot(ax=ax,facecolor="0.97",edgecolor="0.65",linewidth=.35)
        ax.plot(glon,glat,color="0.1",linewidth=.65)
        map_axes(ax,config,geocentric=True);label_panel(ax,f"({string.ascii_lowercase[i]}) {a} km")
    save(fig,OUT,"ground_tracks")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--group",choices=["all","target","workflow","maximum","analysis"],default="all")
    args=p.parse_args();config,nearby=setup()
    if args.group in ["all","target"]:target_map(config,nearby)
    if args.group in ["all","workflow"]:workflow(config)
    if args.group in ["all","maximum","analysis"]:
        records,lon,lat=load_targets(config)
        if args.group in ["all","analysis"]:analysis_plots(config,nearby,records,lon,lat)
        if args.group in ["all","maximum"]:
            maps(config,nearby,records,lon,lat,"max_gap_h","maximum_gap_maps","Maximum complete gap [h]",mark_max=True)


if __name__=="__main__":
    main()
