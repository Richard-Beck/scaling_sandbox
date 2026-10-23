"""Rebuild this layout from unchanged committed source data."""
from pathlib import Path
import runpy
import sys
root=next(p for p in Path(__file__).resolve().parents if (p/'manuscript_outline/outline_v11.json').exists())
sys.argv=[str(root/'manuscript_outline/build_round6_layout_figures_20261023.py'),'--figure','2']
runpy.run_path(sys.argv[0],run_name='__main__')
