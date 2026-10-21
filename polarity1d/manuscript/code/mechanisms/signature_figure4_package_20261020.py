"""Assemble a sourced Figure 4 candidate, then qualify the exact reviewed export.

This module does not choose a framing or recompute a scientific estimand. Each
S0--S12 page comes from a completed scientific renderer and retains an exact CSV.
Candidate generation precedes independent root review. Finalization preserves
the reviewed pixels and writes final integration metadata only after checking
the review's explicit candidate bindings.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from datetime import datetime, timezone
import hashlib
import html
import json
import math
from pathlib import Path
import re
import shutil
import textwrap

csv.field_size_limit(64*1024*1024)

from pypdf import PdfReader, PdfWriter


REQUIRED_EVIDENCE = {"E0", "E1", "E2", "E3", "E4", *{f"A{i}" for i in range(1, 8)},
                     "common_horizon", "long_horizon", "windows_and_coverage", "framing",
                     "independent_final_decision_audit"}
SECTIONS = tuple(f"S{i}" for i in range(13))
CONTENT_FIELDS = ("title", "role", "legend", "narrative", "measurements", "limitations")
PORTABLE_LIMIT=90*1024**2


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def verified_file(record, base):
    path = Path(record["path"])
    if not path.is_absolute():
        path = Path(base) / path
    path = path.resolve()
    if not path.is_file() or sha(path) != record["sha256"]:
        raise ValueError(f"Missing or changed source: {path}")
    return path


def get_field(value, dotted):
    for field in dotted.split("."):
        if not isinstance(value, dict) or field not in value:
            raise ValueError(f"Required receipt field is absent: {dotted}")
        value = value[field]
    return value


def verify_receipt(record, base):
    path = verified_file(record, base)
    payload = json.loads(path.read_text())
    assertions = record.get("assertions", {})
    if not assertions or not any(value is True for value in assertions.values()):
        raise ValueError("Evidence requires an explicit affirmative completion/audit predicate")
    for field, expected in assertions.items():
        actual = get_field(payload, field)
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError(f"Receipt predicate failed: {path}: {field}")
    if not record.get("scope"):
        raise ValueError("Evidence must declare its scientific scope")
    return dict(path=str(path), sha256=sha(path), assertions=assertions, scope=record["scope"])


def repo_path(path, repo):
    try:
        return str(Path(path).resolve().relative_to(Path(repo).resolve()))
    except ValueError as error:
        raise ValueError("Integration assets must be inside the repository") from error


def validate_atlas_coverage(spec,base):
    gate=verify_receipt(dict(spec["atlas_coverage"],assertions=dict(passed=True,all_required_contexts_complete=True),
        scope="All required native-atlas contexts from the frozen canonical campaign census"),base)
    master=json.loads(Path(gate["path"]).read_text());expected={};indices=[]
    for record in master["page_data_indices"]:
        path=verified_file(record,base);payload=json.loads(path.read_text())
        pages=payload["pages"]
        if len(pages)!=record["pages"] or payload["summary"].get("passed") is not True:
            raise ValueError("An atlas context omits declared native pages")
        if record["sha256"] in expected:raise ValueError("The same atlas context was declared twice")
        expected[record["sha256"]]={}
        for page in pages:
            key=(page["pdf_volume"],page["pdf_page"])
            if key in expected[record["sha256"]]:raise ValueError("An atlas source page is duplicated")
            expected[record["sha256"]][key]=page
        indices.append(dict(path=str(path),sha256=record["sha256"],pages=len(pages)))
    if not expected:raise ValueError("The complete atlas context census is required")
    observed={key:set() for key in expected}
    for page in spec["sections"]["S11"]:
        index=page.get("atlas_index_sha256")
        key=(Path(page["pdf"]["path"]).name,page["page"])
        if index not in expected or key not in expected[index] or key in observed[index]:
            raise ValueError("A native atlas page is omitted, duplicated or outside its canonical context")
        source=expected[index][key]
        actual_csvs=page.get("csv_shards") or [page["csv"]]
        if sorted(row["sha256"] for row in actual_csvs)!=sorted(row["sha256"] for row in source["native_field_CSV_shards"]):
            raise ValueError("Every exact plotted native CSV shard must be retained")
        if "condition" in source and (page.get("native_condition")!=source["condition"] or
                                      page.get("selected_cases")!=source["selected_cases"]):
            raise ValueError("Atlas case identities cannot change during assembly")
        metadata={row["sha256"] for row in page.get("page_metadata",[])}
        required={source["selection_CSV"]["sha256"],*[row["sha256"] for row in source["companions"].values()]}
        if not required<=metadata:raise ValueError("Atlas selection/snapshot/scalar/native-grid companions are incomplete")
        observed[index].add(key)
    if any(observed[key]!=set(expected[key]) for key in expected):
        raise ValueError("An entire native atlas volume or context was omitted")
    return gate,indices


def validate_no_clean_mapping_memo(spec,base,producer,independent):
    """A no-clean-mapping result requires the brief's concrete reviewed memo."""
    if type(spec['main'].get('no_clean_mapping')) is not bool:
        raise ValueError('Main framing must explicitly declare no_clean_mapping true or false')
    if spec['main']['no_clean_mapping'] is False:return {}
    if 'reviewed_interpretation_memo' not in spec:
        raise ValueError('No-clean-mapping final package requires the mandatory reviewed S0 memo')
    record=spec['reviewed_interpretation_memo']
    source=verify_receipt(dict(record['source'],assertions=dict(approved_for_candidate_framing=True),
        scope='Root-reviewed no-clean-mapping interpretation, numeric shortfalls and alternatives'),base)
    memo=json.loads(Path(source['path']).read_text())
    if (memo.get('no_clean_mapping') is not True or memo.get('framing')!=spec['main']['framing'] or
            memo.get('producer_complete_sha256')!=producer['sha256'] or
            memo.get('independent_final_decision_audit_sha256')!=independent['sha256'] or
            not memo.get('overall_interpretation') or not memo.get('retained_qualifications')):
        raise ValueError('The mandatory S0 memo must bind the exact framing, completed criteria and independent audit')
    normalize=lambda value:'spring_ring' if value=='spring' else value
    summaries=memo.get('model_summaries',[])
    if len(summaries)!=4 or {normalize(row['model']) for row in summaries}!={'wave_pinning','goryachev','holmes_model3','spring_ring'}:
        raise ValueError('The mandatory S0 memo must cover all four actual model patterns')
    complete=json.loads(Path(producer['path']).read_text());criteria_path=Path(producer['path']).parent/'criteria.csv'
    verified_file(dict(path=str(criteria_path),sha256=complete['outputs']['criteria.csv']),base)
    truth=lambda value:value is True or str(value).lower()=='true'
    with criteria_path.open(newline='') as stream:
        criteria={(row['arm'],row['view'],truth(row['forced_inclusive']),normalize(row['model']),row['criterion']):row
                  for row in csv.DictReader(stream)}
    for row in summaries:
        if not all(row.get(key) for key in ('pattern','numeric_shortfalls','qualifications','criterion_refs')):
            raise ValueError('Every model memo must retain numerical shortfalls, qualifications and actual criterion references')
        for ref in row['criterion_refs']:
            key=(ref['arm'],ref['view'],truth(ref['forced_inclusive']),normalize(ref['model']),ref['criterion'])
            if key not in criteria or ref.get('status')!=criteria[key]['status']:
                raise ValueError('The mandatory S0 memo substituted a different criterion status')
    alternatives=memo.get('alternative_framings',[])
    if not 2<=len(alternatives)<=3 or any(not all(row.get(key) for key in
            ('framing','supporting_panels','weaknesses','further_runs')) for row in alternatives):
        raise ValueError('No-clean-mapping S0 memo requires two or three supported alternatives with weaknesses and resolving runs')
    display=verify_receipt(dict(record['display_audit'],assertions=dict(passed=True,complete=True),
        scope='Every reviewed S0 memo page and exact per-page CSV'),base)
    actual=json.loads(Path(display['path']).read_text())
    if actual.get('source_files',{}).get('reviewed_interpretation_memo',{}).get('sha256')!=source['sha256']:
        raise ValueError('The S0 memo display does not bind the exact reviewed memo bytes')
    expected={row['pdf_page']:row for row in actual['pages']}
    supplied=[page for page in spec['sections']['S0'] if page['display_audit']['sha256']==display['sha256']]
    if (not expected or len(expected)!=actual['pdf'].get('pages') or len(expected)!=len(actual['pages']) or
            len(supplied)!=len(expected) or {page['page'] for page in supplied}!=set(expected)):
        raise ValueError('Every mandatory S0 memo page must appear exactly once in the package')
    for page in supplied:
        if (page['pdf']['sha256']!=actual['pdf']['sha256'] or page.get('csv',{}).get('sha256')!=expected[page['page']]['csv_sha256']):
            raise ValueError('A mandatory S0 memo PDF page or numeric CSV differs from its reviewed display')
    return dict(reviewed_interpretation_memo=source,reviewed_interpretation_memo_display=display)


def validate_spec(spec, base):
    if spec.get("schema") != "figure4_reviewable_package_v1":
        raise ValueError("Unknown Figure 4 package schema")
    evidence = spec.get("evidence", {})
    if set(evidence) != REQUIRED_EVIDENCE:
        raise ValueError(f"Scientific gates must cover exactly {sorted(REQUIRED_EVIDENCE)}")
    sources = {key: verify_receipt(value, base) for key, value in evidence.items()}
    statistics=spec["final_statistics"]
    producer=verify_receipt(dict(statistics["producer"],assertions=dict(passed=True,complete=True),
        scope="Complete frozen common/long H1/H2/H3 criterion producer; independent decision replay required"),base)
    independent=sources["independent_final_decision_audit"]
    audited=json.loads(Path(independent["path"]).read_text())
    if audited.get("producer_complete_sha256")!=producer["sha256"]:
        raise ValueError("Independent final decision audit does not bind the exact criterion producer")
    sources["final_statistics_producer"]=producer
    sections = spec.get("sections", {})
    if set(sections) != set(SECTIONS) or any(not sections[key] for key in SECTIONS):
        raise ValueError("Every S0--S12 section needs at least one audited page")
    for section in SECTIONS:
        for page in sections[section]:
            verify_receipt(page["display_audit"], base)
            if section == "S11" and (page.get("csv_scope") not in ("plotted_values", "native_plotted_values") or
                                      not page.get("page_metadata")):
                raise ValueError("S11 requires actual native plotted-value CSV plus separate selection/snapshot metadata")
    atlas_gate,atlas_indices=validate_atlas_coverage(spec,base)
    sources["S11_full_atlas_coverage"]=atlas_gate
    main = spec.get("main", {})
    minimum_width=main.get("asset_metadata",{}).get("minimum_printed_width_inches",12.)
    if isinstance(minimum_width,bool) or not isinstance(minimum_width,(int,float)) or not math.isfinite(minimum_width) or minimum_width<12.:
        raise ValueError("Figure4 main needs an explicit readable minimum printed width of at least12 inches")
    panels = main.get("panels", [])
    if not 1 <= len(panels) <= 5 or len({p["id"] for p in panels}) != len(panels):
        raise ValueError("The main figure needs one to five distinct actual panels")
    if panels[0].get("id") != "A" or panels[0].get("kind") != "native_small_large_kymographs":
        raise ValueError("Main panel A must show small/large native example kymographs")
    if any(not panel.get("measurement") or not panel.get("source") for panel in panels):
        raise ValueError("Every main panel needs its measured observable and exact source")
    if any(not panel.get("csv") and not panel.get("csv_shards") for panel in panels):
        raise ValueError("Every main panel needs its exact plotted numeric CSV sources")
    if panels[0].get("csv_scope") != "native_plotted_values":
        raise ValueError("Main panel A needs actual native plotted-value sources")
    sources.update(validate_no_clean_mapping_memo(spec,base,producer,independent))
    sources["main_display_audit"]=verify_receipt(main["display_audit"],base)
    for name in ("png", "pdf"):
        verified_file(main[name], base)
    content = spec.get("content", {})
    if any(not content.get(key) for key in CONTENT_FIELDS):
        raise ValueError("Integration prose must include all scientific content fields")
    if any(not isinstance(content[key], list) for key in CONTENT_FIELDS[2:]):
        raise ValueError("Legend, narrative, measurements and limitations must be lists")
    historical = verify_receipt(spec["historical_numerics"], base)
    qualification = json.loads(Path(historical["path"]).read_text())
    if (qualification.get("original_full_gate_passed") is not False or
            qualification.get("stable_scope_accepted") is not True or
            qualification.get("gate_counts", {}).get("all_q90_excluded") != 18 or
            qualification.get("gate_counts", {}).get("stable_scope_passed") != 180):
        raise ValueError("The material historical Spring cadence qualification was lost")
    numerical_rows = [item for item in qualification.get("outputs", [])
                      if Path(item["path"]).name == "all_198_gate_rows.csv"]
    for record in qualification.get("outputs", []):
        verified_file(record, Path(historical["path"]).parent)
    if numerical_rows and not any(item["csv"]["sha256"] == numerical_rows[0]["sha256"]
                                  for item in sections["S12"]):
        raise ValueError("S12 must display and retain all 198 historical numerical gate rows")
    immutable = spec["delivered_interim"]
    verified_file(immutable, base)
    return sources, historical,atlas_indices


def copy_source(record, target, base):
    source = verified_file(record, base)
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if source != target.resolve():
        shutil.copyfile(source, target)
    return dict(path=str(target.resolve()), sha256=sha(target))


def page_csvs(item, section, sequence, data, base):
    """Keep portable tables local and verify full native archives in place."""
    records=item.get("csv_shards") or [dict(item["csv"],rows=item.get("rows",item["csv"].get("rows")))]
    if not records:raise ValueError("A page cannot omit its plotted numeric values")
    portable,archives,sources=[],[],[]
    for index,record in enumerate(records,1):
        source=verified_file(record,base)
        compressed=source.name.endswith(".csv.gz")
        opener=gzip.open if compressed else open
        with opener(source,"rt",newline="") as stream:
            csv_reader=csv.DictReader(stream)
            if not csv_reader.fieldnames or len(csv_reader.fieldnames)!=len(set(csv_reader.fieldnames)) or next(csv_reader,None) is None:
                raise ValueError("Every page requires nonempty complete numeric CSV sources")
            if compressed:
                row_count=record.get("rows")
                if type(row_count) is not int or row_count<1:
                    raise ValueError("Compressed native tables need their independently audited row count")
            else:row_count=1+sum(1 for _ in csv_reader)
        entry=dict(record,path=str(source),sha256=record["sha256"],rows=row_count)
        sources.append(entry)
        distribution=record.get("distribution","portable")
        if distribution=="HPC_native_archive":
            if section not in ("S11","MAIN_A"):raise ValueError("Statistical per-page CSVs must remain portable")
            archives.append(entry)
        elif distribution=="portable":
            part=f"_part{index:03d}" if len(records)>1 else ""
            recompress=not compressed and source.stat().st_size>=PORTABLE_LIMIT
            target=data/(f"figure4_complete_{sequence:04d}_{section}{part}.csv"+(".gz" if compressed or recompress else ""))
            if recompress:
                with source.open("rb") as original,gzip.open(target,"wb",compresslevel=6) as packed:
                    shutil.copyfileobj(original,packed,1024*1024)
                copied=dict(path=str(target),sha256=sha(target),uncompressed_source_sha256=record["sha256"])
            else:copied=copy_source(record,target,base)
            if target.stat().st_size>=PORTABLE_LIMIT:
                raise ValueError("A portable numeric table needs additional lossless shards below90MiB")
            portable.append(dict(copied,rows=row_count,distribution="portable"))
        else:raise ValueError("Unknown numeric CSV distribution")
    inventory=None
    if archives:
        path=data/f"figure4_complete_{sequence:04d}_{section}_HPC_archive_inventory.csv"
        fields=list(dict.fromkeys(key for row in archives for key in row))
        with path.open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(archives)
        inventory=dict(path=str(path),sha256=sha(path),rows=len(archives),distribution="portable_archive_inventory")
    return dict(csv=portable[0] if portable else inventory,csv_files=portable,native_archive_csvs=archives,
                archive_inventory=inventory,source_csvs=sources,rows=sum(row["rows"] for row in sources))


def atlas_index(pages,volumes,metadata,data,core_writer):
    """A compact core index links every exact, separately retained S11 page."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    records=[]
    for page in pages:
        if page["section"]!="S11":continue
        records.append(dict(registry_page=page["page"],title=page["title"],
            appendix=Path(page["export_pdf"]["path"]).name,appendix_sha256=page["export_pdf"]["sha256"],
            appendix_page=page["pdf_page"],numeric_rows=page["rows"],
            numeric_csv_inventory=Path(page["archive_inventory"]["path"]).name if page["archive_inventory"] else "",
            portable_csvs=";".join(Path(row["path"]).name for row in page["csv_files"]),
            selection_snapshot_metadata=";".join(Path(row["path"]).name for row in page["page_metadata"])))
    if not records:raise ValueError("S11 cannot be replaced by an empty appendix index")
    registry=data/"figure4_S11_complete_appendix_page_registry.csv"
    with registry.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    csv_record=dict(path=str(registry),sha256=sha(registry),rows=len(records))
    index_pages=[]
    for start in range(0,len(volumes),24):
        batch=volumes[start:start+24]
        fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis("off")
        fig.suptitle("S11 · Native activity atlas appendices",x=.06,ha="left",y=.97,fontsize=14)
        explanation=(f"All {len(records)} native atlas pages are retained in {len(volumes)} numbered PDF volumes. "
            "This core PDF indexes those appendices. Full plotted native fields remain in verified HPC wideCSV.gz archives; "
            "portable per-page inventories give exact paths, hashes, native grids and row counts. Statistical CSVs and selection/snapshot metadata are portable. "
            "Native profiles retain their saved grids and colour scales; no interpolation or quantization is introduced.")
        ax.text(0,.94,"\n".join(textwrap.wrap(explanation,145)),va="top",fontsize=9)
        cells=[[str(start+i+1),Path(record["path"]).name,str(record["pages"]),f"{Path(record['path']).stat().st_size/1024**2:.2f}"]
               for i,record in enumerate(batch)]
        table=ax.table(cellText=cells,colLabels=["Volume","PDF filename","Native pages","MiB"],
            colWidths=[.10,.63,.15,.12],cellLoc="left",bbox=[0,.18,1,.58])
        table.auto_set_font_size(False);table.set_fontsize(8)
        for i,record in enumerate(batch,1):table[(i,1)].get_text().set_url(Path(record["path"]).name)
        ax.text(0,.10,"Complete per-page registry: "+registry.name+"\nEvery appendix and native-field CSV is hash-verified. The reader links compact archive inventories; full-field files are available on HPC.",fontsize=8,va="top")
        target=metadata/f"figure4_S11_appendix_index_{len(index_pages)+1:02d}.pdf"
        fig.savefig(target);plt.close(fig)
        core_writer.add_page(PdfReader(target).pages[0])
        index_pages.append(dict(section="S11",title="Native atlas appendix index",pdf_page=len(core_writer.pages),
            source_pdf=dict(path=str(target),sha256=sha(target)),csv=csv_record,rows=len(records)))
    return index_pages


def atlas_section(index_path,pdf_root,display_audit):
    """Adapt the atlas owner's complete native-value index without reselection."""
    index_path,pdf_root=Path(index_path).resolve(),Path(pdf_root).resolve()
    envelope=json.loads((index_path.parent/"complete.json").read_text())
    if envelope.get("page_data_index_sha256")!=sha(index_path):
        raise ValueError("The atlas per-page numeric index changed")
    payload=json.loads(index_path.read_text());summary=payload["summary"];pages=payload["pages"]
    if (envelope.get("passed") is not True or summary.get("passed") is not True or
            envelope.get("float64_roundtrip_equal_every_field_time_actualC") is not True or
            envelope.get("sourcefield_resampling") is not False or envelope.get("sourcefield_quantization") is not False or
            len(pages)!=envelope["pages"] or len({page["slot"] for page in pages})!=len(pages)):
        raise ValueError("Every atlas page needs a complete exact native numeric export")
    verify_receipt(display_audit,index_path.parent)
    pdf_hashes={};result=[]
    for page in pages:
        if (page.get("passed") is not True or page.get("decimal_significant_digits")!=17 or
                page.get("float64_roundtrip_equal") is not True or page.get("sourcefield_resampling") is not False or
                page.get("sourcefield_quantization") is not False):
            raise ValueError("Atlas page values must retain native float64 roundtrip and grids")
        pdf=pdf_root/page["pdf_volume"]
        if pdf not in pdf_hashes:pdf_hashes[pdf]=sha(pdf)
        native=[]
        for record in page["native_field_CSV_shards"]:
            if record.get("roundtrip_float64_equal") is not True:raise ValueError("A native CSV shard failed roundtrip")
            native.append(dict(record,distribution="HPC_native_archive"))
        if not native:raise ValueError("An atlas page omits its actual plotted native field values")
        companions=[dict(page["selection_CSV"],label="selection")]
        companions.extend(dict(record,label=label) for label,record in page["companions"].items())
        condition=page["condition"]
        title=(f"{condition['experiment']} · {condition['model']} · C={condition['C_label']:g}, "
            f"f={condition['f']:g}, {condition['kind']}, σ={condition['sigma']:g}; "
            f"horizon {condition['native_horizon']:g}, profile cadence {condition['profile_cadence']:g}")
        result.append(dict(title=title,pdf=dict(path=str(pdf),sha256=pdf_hashes[pdf]),page=page["pdf_page"],
            csv_shards=native,csv_scope="native_plotted_values",page_metadata=companions,display_audit=display_audit,
            native_condition=condition,selected_cases=page["selected_cases"],atlas_index_sha256=sha(index_path)))
    return result


def assemble(spec_path, output, repo):
    spec_path, output, repo = Path(spec_path).resolve(), Path(output).resolve(), Path(repo).resolve()
    spec = json.loads(spec_path.read_text())
    sources, historical,atlas_indices = validate_spec(spec, spec_path.parent)
    output.mkdir(parents=True, exist_ok=True)
    metadata = output / "supporting/metadata"
    data = output / "supporting/data"
    vectors = output / "supporting/vectors"
    for folder in (metadata, data, vectors):
        folder.mkdir(parents=True, exist_ok=True)
    if (metadata / "content.json").exists():
        raise ValueError("A final export already exists; use a fresh candidate destination")
    portable_scientific_sources={}
    for key,record in sources.items():
        target=metadata/f"scientific_evidence_{key}.json"
        portable_scientific_sources[key]=copy_source(record,target,spec_path.parent)
    portable_historical=copy_source(historical,metadata/"historical_Spring_numerical_qualification.json",spec_path.parent)
    immutable = verified_file(spec["delivered_interim"], spec_path.parent)
    supplement = vectors / "figure4_revision_complete_supplement.pdf"
    main_pdf = vectors / "figure4_revision_main.pdf"
    main_png = output / "figure4_revision_main.png"
    if immutable in {supplement.resolve(), main_pdf.resolve(), main_png.resolve()}:
        raise ValueError("The delivered interim cannot be replaced")
    main_records = [copy_source(spec["main"]["png"], main_png, spec_path.parent),
                    copy_source(spec["main"]["pdf"], main_pdf, spec_path.parent)]
    main_source_assets=[dict(spec["main"][key],path=str(verified_file(spec["main"][key],spec_path.parent))) for key in ("png","pdf")]
    main_panel_data=[]
    for sequence,panel in enumerate(spec["main"]["panels"],1):
        numeric=page_csvs(panel,"MAIN_"+panel["id"],sequence,data,spec_path.parent)
        main_panel_data.append(dict(id=panel["id"],measurement=panel["measurement"],source=panel["source"],
                                    csv_scope=panel.get("csv_scope","plotted_values"),**numeric))
    writer, pages, source_pdfs,checked_pdfs,checked_display_audits = PdfWriter(), [], {},{},{}
    appendix_by_source={};appendix_page_use={};index_pages=[]
    for section in SECTIONS:
        writer.add_outline_item(section,len(writer.pages))
        for item in spec["sections"][section]:
            pdf_key=(item["pdf"]["path"],item["pdf"]["sha256"])
            if pdf_key not in checked_pdfs:
                checked_pdfs[pdf_key]=verified_file(item["pdf"],spec_path.parent)
            pdf=checked_pdfs[pdf_key]
            page_number = item["page"]
            if type(page_number) is not int or page_number < 1:
                raise ValueError("Supplement pages use positive one-based PDF page numbers")
            if pdf not in source_pdfs:
                source_pdfs[pdf] = PdfReader(pdf)
            pdf_reader = source_pdfs[pdf]
            if page_number > len(pdf_reader.pages):
                raise ValueError("Requested source PDF page does not exist")
            numeric=page_csvs(item,section,len(pages)+1,data,spec_path.parent)
            page_metadata=[]
            for index,record in enumerate(item.get("page_metadata",[]),1):
                source=verified_file(record,spec_path.parent)
                label=re.sub(r"[^A-Za-z0-9_-]","_",record.get("label",f"source{index:02d}"))
                target=metadata/f"figure4_complete_{len(pages)+1:04d}_{section}_{index:02d}_{label}{''.join(source.suffixes)}"
                page_metadata.append(dict(copy_source(record,target,spec_path.parent),label=record.get("label",source.name),
                    source_path=str(source),source_sha256=record["sha256"]))
            if section=="S11":
                if pdf not in appendix_by_source:
                    if pdf.stat().st_size>=PORTABLE_LIMIT:
                        raise ValueError("S11 source volumes must be lossless PDFs below90MiB")
                    target=vectors/f"figure4_S11_native_atlas_appendix_{len(appendix_by_source)+1:03d}.pdf"
                    appendix_by_source[pdf]=dict(copy_source(item["pdf"],target,spec_path.parent),pages=len(pdf_reader.pages),
                                                label=f"S11 native atlas appendix {len(appendix_by_source)+1:03d}")
                    appendix_page_use[pdf]=set()
                if page_number in appendix_page_use[pdf]:raise ValueError("An S11 native page is declared twice")
                appendix_page_use[pdf].add(page_number)
                export_pdf=appendix_by_source[pdf];export_page=page_number
            else:
                writer.add_page(pdf_reader.pages[page_number - 1])
                export_pdf=dict(path=str(supplement));export_page=len(writer.pages)
            audit_key=(item["display_audit"]["path"],item["display_audit"]["sha256"])
            if audit_key not in checked_display_audits:
                checked_display_audits[audit_key]=verify_receipt(item["display_audit"],spec_path.parent)
            pages.append(dict(section=section, title=item["title"], page=len(pages)+1,
                              source_pdf=dict(path=str(pdf), sha256=item["pdf"]["sha256"], page=page_number),
                              csv_scope=item.get("csv_scope", "plotted_values"),
                              page_metadata=page_metadata,
                              native_condition=item.get("native_condition"),selected_cases=item.get("selected_cases"),
                              atlas_index_sha256=item.get("atlas_index_sha256"),
                              export_pdf=export_pdf,pdf_page=export_page,
                              display_audit=checked_display_audits[audit_key],**numeric))
        if section=="S11":
            for pdf,used in appendix_page_use.items():
                if used!=set(range(1,len(source_pdfs[pdf].pages)+1)):
                    raise ValueError("Every S11 appendix source page must be declared; extra or missing pages are prohibited")
            index_pages=atlas_index(pages,list(appendix_by_source.values()),metadata,data,writer)
    writer.add_metadata({"/Title": spec["content"]["title"],
                         "/Subject": "Completed Figure 4 scientific analyses and source evidence"})
    with supplement.open("wb") as stream:
        writer.write(stream)
    if supplement.stat().st_size>=PORTABLE_LIMIT:
        raise ValueError("The statistical core needs additional lossless volumes below90MiB")
    for page in pages:
        if page["section"]!="S11":page["export_pdf"]["sha256"]=sha(supplement)
    reader_path = output / "figure4_complete_reader.html"
    main_numeric_html='<h2>Main source data</h2><ul>'+''.join(
        '<li>Panel '+html.escape(panel['id'])+': '+html.escape(panel['measurement'])+' · '+
        ' · '.join(f'<a href="supporting/data/{Path(record["path"]).name}">numeric CSV</a>' for record in panel['csv_files'])+
        (f'<a href="supporting/data/{Path(panel["archive_inventory"]["path"]).name}">Full native values: HPC archive inventory</a>' if panel['archive_inventory'] else '')+
        '</li>' for panel in main_panel_data)+'</ul>'
    sections_html = []
    for section in SECTIONS:
        selected = [page for page in pages if page["section"] == section]
        sections_html.append("<h2>" + section + "</h2><ul>" + "".join(
            f'<li><a href="supporting/vectors/{Path(p["export_pdf"]["path"]).name}#page={p["pdf_page"]}">{html.escape(p["title"])}</a> '
            '('+' · '.join(f'<a href="supporting/data/{Path(record["path"]).name}">CSV part{index}</a>' for index,record in enumerate(p["csv_files"],1))+
            ((' · ' if p["csv_files"] else '')+f'<a href="supporting/data/{Path(p["archive_inventory"]["path"]).name}">Full native values: HPC archive inventory</a>' if p["archive_inventory"] else '')+
            ''.join(f' · <a href="supporting/metadata/{Path(record["path"]).name}">{html.escape(record["label"])}</a>' for record in p["page_metadata"])+
            ')</li>' for p in selected) + "</ul>")
    generated = datetime.now(timezone.utc).isoformat()
    historical_text = ("For Spring C0=3, f=1, seed202610192320800, the original full numerical gate failed "
        "12 scalar-observation cadence comparisons. The accepted stable scope passes 180 rows and excludes all "
        "18 q90 multiregion-span rows, including six passing ODE-only comparisons. These are consecutive saved "
        "multiregion-mask durations, not tracked individual lifetimes. The full-window q90 changes from "
        "10 to 7.15 when scalar cadence changes from 0.5 to 0.05. Exact sigma-one replay preserves an existing "
        "path and forcing; it does not establish numerical or observation convergence. Numerical companions "
        "and dense samples are not independent scientific seeds.")
    content_prose = dict(spec["content"])
    if historical_text not in content_prose["limitations"]:
        content_prose["limitations"] = [*content_prose["limitations"], historical_text]
    prose_html = '<h2>Legend</h2>' + ''.join('<p>' + html.escape(text) + '</p>' for text in content_prose["legend"])
    for field in ("narrative", "measurements", "limitations"):
        prose_html += '<h2>' + field.title() + '</h2>' + ''.join('<p>' + html.escape(text) + '</p>' for text in content_prose[field])
    minimum_width=spec["main"].get("asset_metadata",{}).get("minimum_printed_width_inches",12.)
    display_width=math.ceil(96*minimum_width)
    main_asset_metadata=dict(spec["main"].get("asset_metadata",{}),minimum_printed_width_inches=minimum_width,
                             recommended_printed_width_inches=max(14.,minimum_width),minimum_display_width_px=display_width)
    reader_path.write_text('<!doctype html><meta charset="utf-8"><title>' + html.escape(spec["content"]["title"]) +
        '</title><style>body{max-width:1344px;margin:3rem auto;font:16px/1.5 sans-serif}.main-scroll{overflow-x:auto}.main-scroll img{width:100%;min-width:'+str(display_width)+'px;height:auto}</style>' +
        '<h1>' + html.escape(spec["content"]["title"]) + '</h1><p>Scientific analyses and source evidence. Export qualification is recorded in metadata. Generated ' +
        generated + '. Producer stopped.</p><p><a href="supporting/vectors/' + main_pdf.name + '">Main PDF</a> · '
        '<a href="supporting/vectors/' + supplement.name + '">S0–S12 statistical core and native-atlas index</a></p>' +
        '<p>All native atlas pages are in numbered appendix volumes: '+" · ".join(
            f'<a href="supporting/vectors/{Path(record["path"]).name}">{html.escape(record["label"])}</a>' for record in appendix_by_source.values())+
        '.</p><p>Read the main figure at a printed width of at least '+f'{minimum_width:g}'+
        ' inches (14 inches recommended). On smaller screens, scroll the figure horizontally or open the vector PDF.</p><div class="main-scroll"><img alt="Figure 4 main" src="' +
        main_png.name + '"></div>' + main_numeric_html + prose_html + ''.join(sections_html))
    content = dict(content_prose, id="figure4", status="candidate", assets=[dict(
        **main_asset_metadata, path=repo_path(main_png, repo), caption=spec["content"]["title"])],
        main_pdf=repo_path(main_pdf, repo), supplement_pdf=repo_path(supplement, repo), reader=repo_path(reader_path, repo),
        main_source_data=[dict(panel=panel["id"],measurement=panel["measurement"],rows=panel["rows"],
            portable_CSVs=[dict(record,path=repo_path(record["path"],repo)) for record in panel["csv_files"]],
            HPC_archive_inventory=(dict(panel["archive_inventory"],path=repo_path(panel["archive_inventory"]["path"],repo)) if panel["archive_inventory"] else None)) for panel in main_panel_data],
        supplement_scope="S0–S12 statistical core and complete native-atlas appendix index; native pages are in all linked appendix volumes",
        supplement_volumes=[dict(path=repo_path(record["path"],repo),label=record["label"],pages=record["pages"],sha256=record["sha256"])
                            for record in appendix_by_source.values()],
        sources=[*content_prose.get("sources", []), *[dict(label=key, path=repo_path(value["path"], repo), sha256=value["sha256"])
                 for key, value in portable_scientific_sources.items()] + [dict(label="Historical Spring observation-cadence qualification",
                 path=repo_path(portable_historical["path"], repo), sha256=portable_historical["sha256"])]] )
    candidate_content = metadata / "content_candidate.json"
    write_json(candidate_content, content)
    scientific_readme = metadata / "scientific_readme.md"
    scientific_readme.write_text("# " + content["title"] + "\n\nGenerated " + generated + ". Producer stopped.\n\n" +
        "The candidate and final export audits record review status. The delivered E0 interim is retained separately. "
        "Every S0–S12 page has actual plotted-value CSVs. Statistical tables are portable. Full native-atlas fields retain "
        "native grids and 17-significant-digit roundtrip precision in wideCSV.gz in the HPC native archive, with exact per-page "
        "paths/hashes/rows and separate selection/snapshot metadata in the compact export. The reader links these inventories, "
        "not absent Git files. Original NPZ archives remain local/HPC evidence.\n\n" +
        "The supplement PDF is the statistical core and S11 appendix index. All native pages are retained in separately numbered "
        "lossless appendix volumes below90MiB; the reader and content metadata link every volume.\n\n"+
        "\n\n".join("## " + field.title() + "\n\n" + "\n\n".join(content[field])
                       for field in ("legend", "narrative", "measurements", "limitations")))
    records = main_records + [dict(path=str(supplement), sha256=sha(supplement)),
                              dict(path=str(reader_path), sha256=sha(reader_path)),
                              dict(path=str(candidate_content), sha256=sha(candidate_content)),
                              dict(path=str(scientific_readme), sha256=sha(scientific_readme))] + [record for page in pages for record in page["csv_files"]] + [
                                  page["archive_inventory"] for page in pages if page["archive_inventory"]] + [
                                  record for page in pages for record in page["page_metadata"]]+list(appendix_by_source.values())+[
                                  record for page in index_pages for record in (page["csv"],page["source_pdf"])]
    records.extend([*portable_scientific_sources.values(),portable_historical])
    records.extend(record for panel in main_panel_data for record in panel["csv_files"])
    records.extend(panel["archive_inventory"] for panel in main_panel_data if panel["archive_inventory"])
    if any(Path(record["path"]).stat().st_size>=100*1024**2 for record in records):
        raise ValueError("Every compact publication artifact must remain below100MiB")
    report = dict(schema="figure4_candidate_export_v1", generated_utc=generated, passed=True, status="candidate",
        final_export_eligible=False, framing=spec["main"]["framing"], panels=spec["main"]["panels"],
        scientific_evidence=sources, historical_numerics=historical,main_panel_data=main_panel_data,main_source_assets=main_source_assets,
        portable_scientific_evidence=portable_scientific_sources,
        delivered_interim=dict(path=str(immutable), sha256=sha(immutable)), pages=pages, artifacts=records,
        native_archive_artifacts=[record for page in [*pages,*main_panel_data] for record in page["native_archive_csvs"]],
        atlas_source_indices=atlas_indices,
        supplement_core=dict(path=str(supplement),sha256=sha(supplement),pages=len(writer.pages)),
        supplement_volumes=list(appendix_by_source.values()),native_atlas_index_pages=index_pages,
        coverage=dict(registered_scientific_pages=len(pages),native_atlas_pages=sum(p["section"]=="S11" for p in pages),
            statistical_core_pages=len(writer.pages),atlas_index_pages=len(index_pages),
            appendix_pages=sum(record["pages"] for record in appendix_by_source.values()),
            native_pages_are_in_appendices=True,all_registered_pages_have_verified_numeric_csvs=True),
        candidate_content=dict(path=str(candidate_content), sha256=sha(candidate_content)),
        package_inputs=dict(path=str(spec_path), sha256=sha(spec_path)), repo_root=str(repo), source_sha256=sha(__file__))
    write_json(metadata / "candidate_export_audit.json", report)
    return report


def historical_numerics(qualification_path, output):
    """Display the retained full-gate failure without converting it to a pass."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    qualification_path, output = Path(qualification_path).resolve(), Path(output).resolve()
    qualification = json.loads(qualification_path.read_text())
    if (qualification.get("original_full_gate_passed") is not False or
            qualification.get("stable_scope_accepted") is not True):
        raise ValueError("Historical numerical scope must retain both outcomes explicitly")
    for record in qualification["sources"]:
        verified_file(dict(path=record["frozen"], sha256=record["sha256"]), qualification_path.parent)
    tables = {}
    for record in qualification["outputs"]:
        path = verified_file(record, qualification_path.parent)
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        if len(rows) != record["rows"]:
            raise ValueError("Historical numerical row coverage changed")
        tables[path.name] = rows
    rows = tables["all_198_gate_rows.csv"]
    q90 = [row for row in rows if "q90" in row["metric"]]
    stable = [row for row in rows if "q90" not in row["metric"]]
    if (len(rows) != 198 or len(q90) != 18 or len(stable) != 180 or
            sum(row["status"] == "fail" for row in rows) != 12 or
            any(row["status"] != "pass" for row in stable)):
        raise ValueError("Historical full and stable numerical scopes no longer match")
    vectors, data, metadata = (output / "supporting" / name for name in ("vectors", "data", "metadata"))
    for directory in (vectors, data, metadata): directory.mkdir(parents=True, exist_ok=True)
    csv_path = data / "figure4_S12_historical_numerics.csv"
    source_record = next(record for record in qualification["outputs"]
                         if Path(record["path"]).name == "all_198_gate_rows.csv")
    copy_source(source_record, csv_path, qualification_path.parent)
    fig = plt.figure(figsize=(10, 8.2))
    fig.suptitle("S12 · Retained Spring numerical and observation qualification", x=.08, ha="left", y=.97, fontsize=13)
    fig.text(.08, .90, "Spring C₀ = 3 · input fraction f = 1 · seed 202610192320800\n"
             "Original full gate: failed. Accepted stable scope: 180/180 rows; all 18 q90 rows excluded.", fontsize=10)
    ax = fig.add_axes([.10, .48, .82, .33])
    windows = ("full", "first_half", "last_half")
    selected = [next(row for row in rows if row["check"] == "observation_standard" and
                     row["metric"] == "multiregion_run_q90" and row["window"] == window) for window in windows]
    ax.bar([i-.18 for i in range(3)], [float(row["reference"]) for row in selected], width=.35,
           color="#4477aa", label="Scalar cadence 0.5")
    ax.bar([i+.18 for i in range(3)], [float(row["comparison"]) for row in selected], width=.35,
           color="#dd8844", label="Same path at cadence 0.05")
    ax.set_xticks(range(3), ("Full window", "First half", "Last half"))
    ax.set_ylabel("90th percentile multiregion-mask span\n(native time)")
    ax.set_ylim(0, 13); ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    for i, row in enumerate(selected):
        ax.text(i-.18, float(row["reference"])+.15, f"{float(row['reference']):g}", ha="center", fontsize=9)
        ax.text(i+.18, float(row["comparison"])+.15, f"{float(row['comparison']):g}", ha="center", fontsize=9)
    paragraphs = [qualification["definition"], qualification["failure"], qualification["qualification"],
                  qualification["replay_scope"], qualification["replication"]]
    y = .40
    for paragraph in paragraphs:
        wrapped = "\n".join(textwrap.wrap(paragraph, 112))
        fig.text(.08, y, wrapped, va="top", fontsize=8.1)
        y -= .019 * len(wrapped.splitlines()) + .012
    fig.text(.08, .025, "CSV retains all 198 rows, including all 12 failures and all 18 excluded q90 comparisons.\n" + csv_path.name,
             fontsize=7)
    pdf_path = vectors / "figure4_S12_historical_numerics.pdf"
    fig.savefig(pdf_path); plt.close(fig)
    page = dict(title="Retained Spring numerical and observation qualification", page=1,
                pdf=dict(path=str(pdf_path), sha256=sha(pdf_path)),
                csv=dict(path=str(csv_path), sha256=sha(csv_path)), rows=len(rows))
    report = dict(passed=True, original_full_gate_passed=False, stable_scope_accepted=True,
                  qualification=dict(path=str(qualification_path), sha256=sha(qualification_path)),
                  pages=[page], source_sha256=sha(__file__))
    write_json(metadata / "S12_historical_numerics_display_audit.json", report)
    return report


def finalize(candidate_path, review_record):
    candidate_path = Path(candidate_path).resolve()
    candidate = json.loads(candidate_path.read_text())
    if candidate.get("schema") != "figure4_candidate_export_v1" or candidate.get("passed") is not True:
        raise ValueError("Finalization requires a passed concrete candidate export")
    for record in candidate["artifacts"]:
        verified_file(record, candidate_path.parent)
    verified_file(candidate["package_inputs"], candidate_path.parent)
    for record in candidate["atlas_source_indices"]:verified_file(record,candidate_path.parent)
    for record in candidate["scientific_evidence"].values():
        verify_receipt(record, candidate_path.parent)
    checked_sources=set()
    for page in candidate["pages"]:
        verify_receipt(page["display_audit"], candidate_path.parent)
        for record in [page["source_pdf"],*page["source_csvs"]]:
            key=(record["path"],record["sha256"])
            if key not in checked_sources:
                verified_file(record,candidate_path.parent);checked_sources.add(key)
    for record in [*candidate["main_source_assets"],*[record for panel in candidate["main_panel_data"] for record in panel["source_csvs"]]]:
        key=(record["path"],record["sha256"])
        if key not in checked_sources:
            verified_file(record,candidate_path.parent);checked_sources.add(key)
    verify_receipt(candidate["historical_numerics"], candidate_path.parent)
    verified_file(candidate["delivered_interim"], candidate_path.parent)
    review_path = verified_file(review_record, candidate_path.parent)
    if review_record["repo_path"] != repo_path(review_path, candidate["repo_root"]):
        raise ValueError("Completion review must use its exact repository-relative integration path")
    review = json.loads(review_path.read_text())
    if (review.get("approved_for_final_export") is not True or
            review.get("candidate_export_audit_sha256") != sha(candidate_path) or
            review.get("candidate_content_sha256") != candidate["candidate_content"]["sha256"] or
            review.get("scientific_qualifications_retained") is not True or
            review.get("main_and_supplement_visual_review_passed") is not True):
        raise ValueError("Root review must approve the exact candidate and retain its scientific qualifications")
    content = json.loads(Path(candidate["candidate_content"]["path"]).read_text())
    content["status"] = "final"
    content["sources"].append(dict(label="Independent Figure 4 completion review", path=review_record["repo_path"],
                                   sha256=sha(review_path)))
    final_content = candidate_path.parent / "content.json"
    write_json(final_content, content)
    report = dict(status="final", passed=True, final_export_eligible=True,
        generated_utc=datetime.now(timezone.utc).isoformat(),
        candidate_export=dict(path=str(candidate_path), sha256=sha(candidate_path)),
        completion_review=dict(path=str(review_path), sha256=sha(review_path)),
        content=dict(path=str(final_content), sha256=sha(final_content)),
        reviewed_pixels_unchanged=True, historical_full_gate_failure_retained=True,
        artifacts=candidate["artifacts"], native_archive_artifacts=candidate["native_archive_artifacts"],source_sha256=sha(__file__))
    write_json(candidate_path.parent / "final_export_audit.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("candidate")
    build.add_argument("--inputs", required=True, type=Path)
    build.add_argument("--output", required=True, type=Path)
    build.add_argument("--repo", required=True, type=Path)
    final = commands.add_parser("finalize")
    final.add_argument("--candidate-audit", required=True, type=Path)
    final.add_argument("--review", required=True, type=Path)
    final.add_argument("--review-sha256", required=True)
    final.add_argument("--review-repo-path", required=True)
    numerical = commands.add_parser("historical-numerics")
    numerical.add_argument("--qualification", required=True, type=Path)
    numerical.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "candidate":
        result = assemble(args.inputs, args.output, args.repo)
    elif args.command == "finalize":
        result = finalize(args.candidate_audit, dict(path=str(args.review), sha256=args.review_sha256,
                                                   repo_path=args.review_repo_path))
    else:
        result = historical_numerics(args.qualification, args.output)
    print(json.dumps({key: result[key] for key in ("status", "passed", "final_export_eligible", "original_full_gate_passed")
                      if key in result}, indent=2))
