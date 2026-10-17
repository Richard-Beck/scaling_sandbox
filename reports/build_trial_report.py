"""Build the standalone report from completed analyses; never launch simulations."""
from pathlib import Path
import sys
ASSETS=Path(__file__).resolve().parent/'trial_report_assets'
sys.path.insert(0,str(ASSETS))
if __name__=='__main__':
    from build_figures_v3 import figure2,figure3,supplement1
    figure2();figure3();supplement1()
    from zmurchok2020.summarize import main as spring_summary
    spring_summary()
    from render_v3 import main
    main()
