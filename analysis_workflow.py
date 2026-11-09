"""Generate the paper workflow. Run: python analysis_workflow.py"""
from paper_figures import setup, workflow

if __name__ == "__main__":
    config, _ = setup()
    workflow(config)
