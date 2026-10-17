#!/usr/bin/env python3
"""Build only retained WP1 Figure 1, Figure S1, and Figure S2."""

from __future__ import annotations

from pathlib import Path

import build_wp1_manuscript_figures as figures


def main() -> None:
    output = figures.ROOT / "results" / "manuscript_wp1_v3"
    source = output / "source_data"
    output.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    figures._style()
    figures.build_figure1(output, source)
    figures.build_supplement_otsuji(output, source)
    figures.build_supplement_tracking(output, source)
    print(output)


if __name__ == "__main__":
    main()
