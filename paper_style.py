"""Shared publication figure settings. Change appearance here, not the data."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

LINE_STYLES = ["-", "--", "-.", ":"]


def apply_style():
    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.linewidth": 0.6, "lines.linewidth": 0.9,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "axes.edgecolor": "0.15", "axes.labelcolor": "0.1",
        "text.color": "0.1", "xtick.color": "0.15", "ytick.color": "0.15",
        "legend.frameon": False, "figure.facecolor": "white",
        "axes.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def label_panel(ax, label):
    ax.set_title(label, loc="left", pad=5, fontsize=10, fontweight="normal")


def light_grid(ax, axis="y"):
    ax.grid(axis=axis, color="0.85", linewidth=0.35, zorder=0)
    ax.set_axisbelow(True)


def save(fig, output_dir, name):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "png"]:
        fig.savefig(output_dir / f"{name}.{ext}", dpi=240,
                    bbox_inches="tight", pad_inches=0.045)
    plt.close(fig)
    print(f"Saved {name}.pdf and .png", flush=True)
