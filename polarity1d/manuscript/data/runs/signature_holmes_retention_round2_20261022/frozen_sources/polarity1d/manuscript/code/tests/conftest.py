from pathlib import Path
import sys
CORE=Path(__file__).resolve().parents[1]/'core'
sys.path.insert(0,str(CORE))
sys.path.insert(0,str(CORE/'src'))
