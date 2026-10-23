"""Use the exact native Figure 2 PDF page in otherwise unchanged outline layout.

The normal renderer supplies notes, links, bookmarks and input checks. Only an
asset explicitly declaring native_pdf_page replaces its temporary raster sheet,
before the final PDF and build manifest are generated. This keeps the 7.2 × 10
inch figure at its designed font sizes without adding a title to its artwork.
"""
from pathlib import Path
import sys
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import build_manuscript_outline_v2 as revised


class NativeFigureBuilder(revised.RevisionBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.track(Path(__file__), 'native_figure_page_renderer')
        self.native_pages = {}

    def figure(self, section, asset, index, count):
        source = asset.get('native_pdf_page')
        if not source:
            return super().figure(section, asset, index, count)
        path = self.resolve(source)
        self.track(path, 'native_vector_figure')
        native = PdfReader(path)
        assert len(native.pages) == 1
        page = native.pages[0]
        dimensions = [float(page.mediabox.width)/72, float(page.mediabox.height)/72]
        assert all(abs(a-b) < 1e-9 for a, b in zip(dimensions, asset['page_size']))
        assert dimensions[0] + 1e-9 >= asset['minimum_printed_width_inches']
        placeholder = dict(asset, minimum_printed_width_inches=0)
        super().figure(section, placeholder, index, count)
        # figure() finishes the previous page before opening this figure sheet.
        self.native_pages[len(self.pages)] = dict(path=path, dimensions=dimensions)

    def build(self, destination):
        super().build(destination)
        rendered = PdfReader(destination)
        writer = PdfWriter()
        for index, page in enumerate(rendered.pages):
            native = self.native_pages.get(index)
            if native:
                original = PdfReader(native['path']).pages[0]
                writer.add_page(original)
                self.pages[index].update(size_in=native['dimensions'],
                    printed_figure_width_inches=native['dimensions'][0], text_blocks=0,
                    native_vector_pdf=str(native['path'].relative_to(ROOT)),
                    native_page_content_preserved=True)
                self.links[index] = []
            else:
                writer.add_page(page)
        writer.add_metadata({key: str(value) for key, value in rendered.metadata.items()})
        temporary = destination.with_suffix('.native.pdf')
        writer.write(temporary)
        temporary.replace(destination)


if __name__ == '__main__':
    revised.base.OutlineBuilder = NativeFigureBuilder
    revised.base.main()
