#!/home/4473331/.conda/envs/scaling_sandbox/bin/python
"""Single entry point for the completed Figure 4 presentation rebuild.

The executed numerical analysis is archived unchanged in
supporting/metadata/final_analysis_source_snapshot/rebuild_figure4.py.
This post-review entry point reads frozen numerical exports and runs no analysis
or native models. See supporting/metadata/postreview_presentation_audit.json.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent/'scripts'))
from presentation_postreview import main
if __name__=='__main__':main()
