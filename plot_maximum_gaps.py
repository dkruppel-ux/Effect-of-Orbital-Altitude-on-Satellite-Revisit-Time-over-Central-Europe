"""Generate maximum-gap maps from saved results. Run: python plot_maximum_gaps.py"""
from paper_figures import setup, load_targets, maps

if __name__ == "__main__":
    config, nearby = setup()
    records, lon, lat = load_targets(config)
    maps(config, nearby, records, lon, lat, "max_gap_h", "maximum_gap_maps",
         "Maximum complete gap [h]", mark_max=True)
