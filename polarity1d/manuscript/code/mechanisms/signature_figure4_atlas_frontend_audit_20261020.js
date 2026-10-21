// Deterministic DOM harness; this is not a browser screenshot/render audit.
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const crypto = require('crypto');
const root = path.resolve(process.argv[2]);
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
class Element {
  constructor(tag) { this.tag = tag; this.children = []; this.value = ''; this.listeners = {}; }
  append(...children) { this.children.push(...children); }
  replaceChildren(...children) { this.children = children; if (this.tag === 'select') this.value = children.length ? String(children[0].value) : ''; }
  addEventListener(name, callback) { this.listeners[name] = callback; }
  removeAttribute(name) { delete this[name]; }
}
const elements = Object.fromEntries(['search', 'page', 'title', 'image', 'files', 'archives'].map(id => [id, new Element(id === 'page' ? 'select' : 'div')]));
const document = { getElementById: id => elements[id], createElement: tag => new Element(tag) };
const context = vm.createContext({ document });
vm.runInContext(script, context);
const requireTrue = (condition, message) => { if (!condition) throw new Error(message); };
const pages = vm.runInContext('pages', context);
const pageIndex = JSON.parse(fs.readFileSync(path.join(root, 'page_data', 'page_data_index.json'), 'utf8'));
const expectedPages = pageIndex.summary.pages;
requireTrue(pages.length === expectedPages && expectedPages > 0, 'Expected complete frozen-context page census');
let portableLinks = 0;
for (let i = 0; i < pages.length; i++) {
  elements.page.value = String(i); vm.runInContext('show()', context);
  requireTrue(elements.title.textContent === pages[i].title, 'Selection/title mismatch');
  requireTrue(elements.image.src === pages[i].image, 'Selection/image mismatch');
  requireTrue(fs.existsSync(path.join(root, elements.image.src)), 'Missing portable image');
  for (const item of elements.files.children[0].children) {
    const anchor = item.children[0];
    requireTrue(fs.existsSync(path.join(root, anchor.href.split('#')[0])), 'Missing portable CSV/PDF'); portableLinks++;
  }
  for (const item of elements.archives.children[0].children) {
    requireTrue(item.children[0].tag === 'code' && !item.children[0].href, 'HPC source path must be plain provenance, not an absent GitHub href');
  }
}
for (const query of ['Spring', 'Holmes 3', 'Wave-pinning', 'Goryachev', '0.0001', 'mean', 'zero']) {
  elements.search.value = query; vm.runInContext('filter()', context);
  const expected = pages.filter(p => p.title.toLowerCase().includes(query.toLowerCase())).length;
  requireTrue(elements.page.children.length === expected, 'Filter coverage mismatch');
}
elements.search.value = 'no matching condition at all'; vm.runInContext('filter()', context);
requireTrue(elements.page.children.length === 0 && !elements.image.src && elements.files.children.length === 0 && elements.title.textContent === 'No matching condition', 'Empty filters must clear stale data');
elements.search.value = ''; vm.runInContext('filter()', context); requireTrue(elements.page.children.length === expectedPages, 'Reset must restore every page');
const result = { passed: true, scope: 'Node DOM harness: generated JavaScript syntax, all page/image/source selections, physical-asset links, filtering and no-results reset; not a browser screenshot/render audit.', pages: pages.length, portableLinks, noAbsentHPCArchiveHrefs: true, generated_utc: new Date().toISOString(), reader_sha256: crypto.createHash('sha256').update(html).digest('hex'), audit_source_sha256: crypto.createHash('sha256').update(fs.readFileSync(__filename)).digest('hex') };
fs.copyFileSync(__filename, path.join(root, 'metadata', 'frozen_frontend_audit.js'));
fs.writeFileSync(path.join(root, 'metadata', 'frontend_audit.json'), JSON.stringify(result, null, 2) + '\n');
process.stdout.write(JSON.stringify(result) + '\n');
