"""Build the revised outline with explicitly sized, readable continued sheets.

The first-edition renderer remains unchanged. Optional asset page_size and
page_title fields control presentation only; source tracking and validation
still run through the original builder.
"""
from pathlib import Path

import build_manuscript_outline as base


class RevisionBuilder(base.OutlineBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.track(Path(__file__), 'revision_renderer')
        self._asset = None
        self._printed_width = None

    def new_page(self, title, kind='prose', landscape=False, continued=False):
        custom = self._asset.get('page_size') if self._asset and kind == 'figure' else None
        if custom is None:
            return super().new_page(title, kind, landscape, continued)
        if len(custom) != 2 or not all(6 <= value <= 24 for value in custom):
            raise ValueError(f'Invalid figure page size: {custom}')
        original_page, original_landscape = base.PAGE, base.LANDSCAPE
        try:
            base.PAGE = base.LANDSCAPE = tuple(custom)
            return super().new_page(title, kind, landscape, continued)
        finally:
            base.PAGE, base.LANDSCAPE = original_page, original_landscape

    def figure(self, section, asset, index, count):
        self._asset = asset
        display = dict(section, title=asset.get('page_title', section['title']))
        try:
            super().figure(display, asset, index, count)
            self._printed_width = self.fig.axes[-1].get_position().width * self.w
            minimum = asset.get('minimum_printed_width_inches', 0)
            if self._printed_width + 1e-6 < minimum:
                raise ValueError(f"Figure {asset['path']} prints at {self._printed_width:.3f}in; minimum {minimum}in")
        finally:
            self._asset = None

    def finish_page(self):
        before = len(self.pages)
        width = self._printed_width
        super().finish_page()
        if len(self.pages) > before and width is not None:
            self.pages[-1]['printed_figure_width_inches'] = width
        self._printed_width = None


if __name__ == '__main__':
    base.OutlineBuilder = RevisionBuilder
    base.main()
