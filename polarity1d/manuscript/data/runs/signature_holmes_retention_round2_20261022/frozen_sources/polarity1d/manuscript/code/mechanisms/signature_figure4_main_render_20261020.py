"""Compose a reviewed Figure 4 main candidate from exact measured sources.

This adapter neither selects a framing nor chooses examples. It reads the
accepted S11 selections and round-trip native CSVs named by the root's framing
record. Statistical panels retain their separately audited vector PDFs and
numeric CSVs. A concrete main candidate still requires independent root review.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pypdf import PdfReader, PdfWriter, Transformation

import signature_figure4_atlas_20261020 as Atlas
import signature_figure4_interim_render_20261020 as P
import signature_figure4_package_20261020 as Q
import signature_figure4_stage_render_20261020 as S


def rows(record, base):
    path = Q.verified_file(record, base)
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("Native source CSV needs unique declared columns")
        yield from reader


def binary_sha(value):
    return hashlib.sha256(np.ascontiguousarray(value, dtype=np.float64).tobytes()).hexdigest()


def selected_native(index_record, slot, case_id, base):
    """Read one already selected history, preserving every saved float64 value."""
    path = Q.verified_file(index_record, base)
    payload = json.loads(path.read_text())
    if payload.get("summary", {}).get("passed") is not True:
        raise ValueError("A native page needs its accepted round-trip index")
    matches = [page for page in payload["pages"] if page["slot"] == slot]
    if len(matches) != 1:
        raise ValueError("The exact declared native atlas slot must exist once")
    page = matches[0]
    if (page.get("passed") is not True or page.get("float64_roundtrip_equal") is not True or
            page.get("sourcefield_resampling") is not False or page.get("sourcefield_quantization") is not False or
            page.get("decimal_significant_digits") != 17 or case_id not in page["selected_cases"]):
        raise ValueError("The native example must retain an existing exact selected case")
    selections = [row for row in rows(page["selection_CSV"], path.parent) if row["case_id"] == case_id]
    if len(selections) != 1:
        raise ValueError("The native example identity is absent or duplicated")
    selection = selections[0]; case_index = int(selection["case_index"])
    expected = page["expected_native_binary_hashes"][str(case_index)]
    time=np.empty(expected["rows"],dtype=np.float64)
    circumference=np.empty_like(time)
    field=np.empty((expected["rows"],expected["cells"]),dtype=np.float64)
    seen=np.zeros(expected["rows"],dtype=bool);source_records=[]
    for source in page["native_field_CSV_shards"]:
        if source.get("roundtrip_float64_equal") is not True:
            raise ValueError("Every original native CSV shard must pass float64 round trip")
        count = 0
        for row in rows(source, path.parent):
            count += 1
            if int(row["case_index"]) == case_index:
                # Python float parses the17-significant-digit text by correct
                # round-to-nearest. The independent native binary hashes below
                # are decisive; no pandas default-precision parser is used.
                fields = sorted(key for key in row if key.startswith("active_cell_"))
                if fields != [f"active_cell_{i:04d}" for i in range(expected["cells"])]:
                    raise ValueError("Native cells were omitted, renamed or resampled")
                frame=int(row["native_frame"])
                if frame<0 or frame>=expected["rows"] or seen[frame]:
                    raise ValueError("Native frames were duplicated or lie outside the declared history")
                time[frame]=float(row["time"]);circumference[frame]=float(row["actualC"])
                field[frame]=np.fromiter((float(row[key]) for key in fields),dtype=np.float64,count=expected["cells"])
                seen[frame]=True
        if count != source["rows"]:
            raise ValueError("Native CSV shard row coverage changed")
        source_records.append(dict(source, distribution="HPC_native_archive",
                                   selected_case_id=case_id, selected_case_index=case_index))
    if not seen.all():raise ValueError("Native frames were omitted")
    if (field.shape != (expected["rows"], expected["cells"]) or not np.isfinite(field).all() or
            not np.isfinite(time).all() or np.any(np.diff(time) <= 0) or
            not np.isfinite(circumference).all() or np.any(circumference <= 0)):
        raise ValueError("Native history must be complete, finite and increasing in time")
    actual_hashes = dict(field_sha256=binary_sha(field), time_sha256=binary_sha(time),
                         actualC_sha256=binary_sha(circumference))
    if any(actual_hashes[key] != expected[key] for key in actual_hashes):
        raise ValueError("Main native values differ from the independent native binary hashes")
    qrows = [row for row in rows(page["companions"]["q_grid"], path.parent)
             if int(row["case_index"]) == case_index]
    qrows.sort(key=lambda row:int(row["native_cell_index"]))
    if [int(row["native_cell_index"]) for row in qrows] != list(range(expected["cells"])):
        raise ValueError("The native material grid is incomplete")
    q = np.asarray([float(row["q"]) for row in qrows])
    np.testing.assert_allclose(q, (np.arange(len(q))+.5)/len(q), rtol=0, atol=1e-12)
    if (not math.isclose(time[-1], float(selection["native_horizon"]), abs_tol=1e-8) or
            not np.allclose(np.diff(time), float(selection["profile_cadence"]), rtol=0, atol=1e-8)):
        raise ValueError("The example native horizon or profile cadence changed")
    evidence = dict(index=dict(path=str(path), sha256=P.sha(path)), slot=slot, case_id=case_id,
        selected_case_index=case_index, condition=page["condition"], selection=selection,
        selected_values_binary_sha256=actual_hashes, native_field_CSV_shards=source_records,
        q_grid=page["companions"]["q_grid"], selection_CSV=page["selection_CSV"],
        native_frames=len(time), native_cells=field.shape[1], native_values=len(time)*field.shape[1],
        saved_field_and_geometry_array_bytes=field.nbytes+time.nbytes+circumference.nbytes+q.nbytes,
        expanded_physical_face_vertices=(2*len(time)-1)*(field.shape[1]+1),
        minimum_actualC=float(circumference.min()), maximum_actualC=float(circumference.max()))
    return dict(time=time, C=circumference, field=field, q=q, evidence=evidence)


def paired_limits(small, large):
    a, b = small["evidence"]["selection"], large["evidence"]["selection"]
    if (P.model(a) != P.model(b) or a["active_field"] != b["active_field"] or
            a["kind"] != "noise" or b["kind"] != "noise" or
            any(not math.isclose(float(a[key]), float(b[key])) for key in ("f", "sigma")) or
            not float(a["C_label"]) < float(b["C_label"])):
        raise ValueError("Each native pair needs the same model/field/noise/input and increasing size")
    if (not math.isclose(small['time'][0],large['time'][0],abs_tol=1e-8) or
            not math.isclose(small['time'][-1],large['time'][-1],abs_tol=1e-8)):
        raise ValueError('Compared main examples need the same full native observation horizon; no hidden clock rescaling')
    low = min(0., float(small["field"].min()), float(large["field"].min()))
    high = max(float(small["field"].max()), float(large["field"].max()))
    if high <= low:
        raise ValueError("A constant null native field cannot be a qualified main example")
    return dict(vmin=low, vmax=high, strategy="Full saved native range shared within this model/field pair",
                physical_arc_extent=[0.,max(float(small['C'].max()),float(large['C'].max()))],
                native_time_extent=[float(small['time'][0]),float(small['time'][-1])],
                native_binary_source_hashes=[movie["evidence"]["selected_values_binary_sha256"]
                                             for movie in (small, large)])


def native_panel(pairs, output, viewports=None):
    """Render the exact atlas face geometry; never interpolate concentration."""
    if not 1 <= len(pairs) <= 2:
        raise ValueError("Panel A permits one or two reviewed native pairs at the declared12-inch minimum main width")
    fig, axes = plt.subplots(len(pairs), 2, figsize=(12, 3.4*len(pairs)+.2), squeeze=False)
    records = []
    for i, (small, large) in enumerate(pairs):
        limits = paired_limits(small, large)
        viewport=(viewports or {}).get(P.model(small['evidence']['selection']),limits['native_time_extent'])
        if (len(viewport)!=2 or not limits['native_time_extent'][0]<=viewport[0]<viewport[1]<=limits['native_time_extent'][1]):
            raise ValueError('A reviewed viewport must lie inside both full native histories')
        limits['display_native_time_extent']=list(viewport)
        for j, movie in enumerate((small, large)):
            selection = movie["evidence"]["selection"]; ax = axes[i,j]
            x, t, nearest = Atlas.physical_faces(movie["time"], movie["C"], len(movie["q"]))
            image = ax.pcolormesh(x, t, movie["field"][nearest], shading="flat", cmap="viridis",
                                 vmin=limits["vmin"], vmax=limits["vmax"], rasterized=True)
            ax.axhline(.5*movie["time"][-1], color="white", ls="--", lw=.7)
            ax.set_xlim(*limits['physical_arc_extent'])
            ax.set_ylim(*reversed(limits['display_native_time_extent']))
            label = P.LABELS[P.model(selection)]
            role=str(selection.get('role','')).lower()
            display_role='medoid' if 'medoid' in role else 'distinct seed' if 'random' in role else 'selected history'
            ax.set_title(f"{label} · C₀={float(selection['C_label']):g} · f={float(selection['f']):g}\n"
                         f"Atlas slot {movie['evidence']['slot']} · {display_role} · Δt={float(selection['profile_cadence']):g}", fontsize=13)
            ax.set_xlabel("Physical arc from fixed material origin", fontsize=12.5)
            ax.set_ylabel("Native time" + (" (Spring)" if P.model(selection)=="spring_ring" else ""), fontsize=12.5)
            ax.tick_params(labelsize=12)
            records.append(dict(movie["evidence"], displayed_role="small" if j==0 else "large",
                shared_color_limits=limits, physical_faces_source=dict(path=str(Path(Atlas.__file__).resolve()),
                sha256=P.sha(Atlas.__file__)), concentration_interpolated=False,
                geometry="Exact native C extrema retained as temporal faces; geometry-only interpolation; nearest native concentration strips",
                paired_axis_scales='Shared physical arc and native time limits; blank outside the actual cell remains outside the domain',
                full_native_time_extent=limits['native_time_extent'],display_native_time_extent=limits['display_native_time_extent'],
                viewport_changes_activity_values_or_regime_analysis=False,
                discarded_first_half_boundary=.5*movie["time"][-1]))
        colorbar=fig.colorbar(image, ax=axes[i,:].tolist(), fraction=.022, pad=.025,
                             label=small["evidence"]["selection"]["active_field"])
        colorbar.ax.tick_params(labelsize=12)
        colorbar.ax.set_ylabel(small["evidence"]["selection"]["active_field"],fontsize=12.5)
    fig.subplots_adjust(left=.09, right=.85, top=.86, bottom=.17, hspace=.65, wspace=.35)
    bounds = S.presentation_bounds(fig)
    if not bounds["passed"]:
        plt.close(fig)
        raise ValueError("Native example axes exceed the panel bounds: "+str(
            [row for row in bounds["items"] if not row["within_page"]]))
    fig.savefig(output, dpi=220); plt.close(fig)
    return dict(path=str(output), sha256=P.sha(output), pages=1), records, bounds


def validated_inputs(spec, base):
    if spec.get("schema") != "figure4_main_candidate_v1":
        raise ValueError("Unknown main-candidate input schema")
    producer = Q.verify_receipt(dict(spec["final_statistics"]["producer"], assertions=dict(passed=True,complete=True),
        scope="Complete final criterion producer"), base)
    audit = Q.verify_receipt(dict(spec["final_statistics"]["independent_audit"], assertions=dict(passed=True),
        scope="Independent final decision replay"), base)
    if json.loads(Path(audit["path"]).read_text()).get("producer_complete_sha256") != producer["sha256"]:
        raise ValueError("Independent audit must bind the exact completed final statistics")
    framing = Q.verify_receipt(dict(spec["framing_decision"], assertions=dict(approved_for_candidate_framing=True),
        scope="Root-reviewed actual final criteria/window/horizon/framing decision"), base)
    decision = json.loads(Path(framing["path"]).read_text())
    parent_memo=None
    if 'parent_interpretation_memo' in decision:
        parent_memo=Q.verify_receipt(dict(decision['parent_interpretation_memo'],assertions=dict(approved_for_candidate_framing=True),
            scope='Immutable reviewed interpretation and ordered native identities, before measured-panel source review'),base)
        original=json.loads(Path(parent_memo['path']).read_text())
        if any(decision.get(key)!=value for key,value in original.items()):
            raise ValueError('The main source review must preserve every field of the immutable interpretation memo')
    if (decision.get("producer_complete_sha256") != producer["sha256"] or
            decision.get("independent_final_decision_audit_sha256") != audit["sha256"] or
            decision.get("framing") != spec["framing"] or type(spec.get('no_clean_mapping')) is not bool or
            decision.get('no_clean_mapping') is not spec['no_clean_mapping']):
        raise ValueError("The framing record must bind the exact scientific decision sources")
    declared_viewports=({row['model']:[row['start'],row['end']] for row in decision['native_display_windows']}
                        if 'native_display_windows' in decision else decision.get('native_viewports',{}))
    if spec.get('native_viewports',{})!=declared_viewports:
        raise ValueError('The native viewport must match the explicit root-reviewed display declaration')
    requested = [item[side]["case_id"] for item in spec["native_pairs"] for side in ("small","large")]
    if requested != decision.get("selected_native_case_ids") or len(set(requested)) != len(requested):
        raise ValueError("Main native identities must follow the reviewed ordered selection")
    for pair in spec["native_pairs"]:
        for side in ("small","large"):
            item=pair[side]
            gate=Q.verify_receipt(item["atlas_completion_gate"],base)
            gate_payload=json.loads(Path(gate["path"]).read_text())
            index_field=item.get("atlas_index_binding_field","native_CSV_page_index_sha256")
            if Q.get_field(gate_payload,index_field)!=item["atlas_index"]["sha256"]:
                raise ValueError("The native atlas completion gate must bind the selected exact page index")
    panels=spec.get("statistical_panels",[])
    if len(panels)>4 or [row["id"] for row in panels]!=list("BCDE")[:len(panels)]:
        raise ValueError("Use at most four sequential measured statistical panels after A")
    for panel in panels:
        if panel.get("small_large_comparison") is not True or not panel.get("measurement"):
            raise ValueError("Every measured main comparison must retain small versus large")
        Q.verify_receipt(panel["display_audit"],base);Q.verified_file(panel["pdf"],base)
        if len(PdfReader(Q.verified_file(panel["pdf"],base)).pages)!=1:
            raise ValueError("Each main measurement component must be a single vector page")
        for record in panel["csvs"]:Q.verified_file(record,base)
        if not panel["csvs"]:raise ValueError("Every main measurement needs actual numeric sources")
    panel_sources=[dict(id=panel["id"],pdf_sha256=panel["pdf"]["sha256"],
        CSV_sha256s=[record["sha256"] for record in panel["csvs"]],
        display_audit_sha256=panel["display_audit"]["sha256"]) for panel in panels]
    if decision.get("statistical_panel_sources")!=panel_sources:
        raise ValueError("Statistical panels must match the exact reviewed framing source identities")
    return dict(producer=producer, independent_audit=audit, framing=framing,parent_interpretation_memo=parent_memo)


def render(spec_path, output, name):
    spec_path,output=Path(spec_path).resolve(),Path(output).resolve();spec=json.loads(spec_path.read_text())
    gates=validated_inputs(spec,spec_path.parent)
    vectors=output/"supporting/vectors";metadata=output/"supporting/metadata";data=output/"supporting/data"
    for directory in (vectors,metadata,data):directory.mkdir(parents=True,exist_ok=True)
    pdf=vectors/f"{name}.pdf";png=output/f"{name}.png";panel_pdf=vectors/f"{name}_panelA.pdf"
    if any(path.exists() for path in (pdf,png,panel_pdf)):raise ValueError("Use a fresh named candidate closure")
    if spec.get('reuse_native_panel'):
        prior_path=Q.verified_file(spec['reuse_native_panel'],spec_path.parent);prior=json.loads(prior_path.read_text())
        if prior.get('passed') is not True or prior.get('complete') is not True:raise ValueError('Native panel reuse needs complete source/presentation proof')
        original=prior['native_panel'];identities=original['selected_histories'];bounds=original['presentation_bounds']
        requested=[pair[side]['case_id'] for pair in spec['native_pairs'] for side in ('small','large')]
        if requested!=[row['case_id'] for row in identities]:raise ValueError('Native panel reuse changed ordered selected identities')
        for identity in identities:
            declared=spec.get('native_viewports',{}).get(P.model(identity['selection']),identity['shared_color_limits']['native_time_extent'])
            if list(declared)!=identity.get('display_native_time_extent',identity['shared_color_limits']['native_time_extent']):
                raise ValueError('A changed reviewed viewport requires a fresh native panel')
            for record in [*identity['native_field_CSV_shards'],identity['q_grid'],identity['selection_CSV']]:Q.verified_file(record,spec_path.parent)
        native=Q.copy_source(original['pdf'],panel_pdf,spec_path.parent);native['pages']=1
    else:
        pairs=[tuple(selected_native(pair[side]["atlas_index"],pair[side]["slot"],pair[side]["case_id"],spec_path.parent)
                     for side in ("small","large")) for pair in spec["native_pairs"]]
        native,identities,bounds=native_panel(pairs,panel_pdf,spec.get('native_viewports'))
    selection_csv=data/f"{name}_panelA_native_selection.csv"
    flat=[];native_files={};companions={}
    for identity in identities:
        flat.append(dict(case_id=identity["case_id"],displayed_role=identity["displayed_role"],
            atlas_slot=identity["slot"],atlas_index_path=identity["index"]["path"],atlas_index_sha256=identity["index"]["sha256"],
            selection=json.dumps(identity["selection"],sort_keys=True),
            selected_values_binary_sha256=json.dumps(identity["selected_values_binary_sha256"],sort_keys=True),
            shared_color_limits=json.dumps(identity["shared_color_limits"],sort_keys=True),
            native_frames=identity["native_frames"],native_cells=identity["native_cells"],
            minimum_actualC=identity["minimum_actualC"],maximum_actualC=identity["maximum_actualC"],
            full_native_time_extent=json.dumps(identity.get('full_native_time_extent',identity['shared_color_limits']['native_time_extent'])),
            display_native_time_extent=json.dumps(identity.get('display_native_time_extent',identity['shared_color_limits']['native_time_extent'])),
            exact_native_CSV_shards=json.dumps(identity["native_field_CSV_shards"],sort_keys=True),
            concentration_interpolated=False,geometry=identity["geometry"]))
        for record in identity["native_field_CSV_shards"]:
            key=(record["path"],record["sha256"])
            native_files.setdefault(key,dict(record,selected_case_ids=[]))["selected_case_ids"].append(identity["case_id"])
        for key in ("q_grid","selection_CSV"):
            record=identity[key];companions[(record["path"],record["sha256"])]=record
    with selection_csv.open("w",newline="") as stream:
        csv_writer=csv.DictWriter(stream,fieldnames=list(flat[0]));csv_writer.writeheader();csv_writer.writerows(flat)
    selection_record=dict(path=str(selection_csv),sha256=P.sha(selection_csv),rows=len(flat))
    components=[dict(id="A",pdf=native)]+spec["statistical_panels"]
    width=1008.
    # The reserved label strip keeps panel letters outside the exact measured
    # component bounds. Vector statistical artwork is retained in the PDF.
    if len(components)==4:
        # Four-panel framing: the complete all-input overview spans the page,
        # while the two selected comparisons retain independent compact slots.
        def natural_height(component,slot_width):
            source=PdfReader(Q.verified_file(component['pdf'],spec_path.parent)).pages[0]
            return (slot_width-20.)*float(source.mediabox.height)/float(source.mediabox.width)
        full_width=width-50.;half_width=(width-75.)/2
        top_height=natural_height(components[0],full_width);middle_height=natural_height(components[1],full_width)
        bottom_height=max(natural_height(item,half_width) for item in components[2:])
        height=50.+top_height+middle_height+bottom_height+70.
        rectangles=[(25.,25.+bottom_height+35.+middle_height+35.,full_width,top_height),
                    (25.,25.+bottom_height+35.,full_width,middle_height),
                    (25.,25.,half_width,bottom_height),(50.+half_width,25.,half_width,bottom_height)]
    else:
        height=900.
        rectangles=[(25.,height*.51,width-50.,height*.45)]
        rectangles.extend([(25.+(i%2)*width*.50,25.+(1-i//2)*height*.23,width*.46,height*.215) for i in range(4)])
    writer=PdfWriter();page=writer.add_blank_page(width=width,height=height)
    for component,rect in zip(components,rectangles):
        source=PdfReader(Q.verified_file(component["pdf"],spec_path.parent)).pages[0]
        x,y,w,h=rect;sw,sh=float(source.mediabox.width),float(source.mediabox.height);scale=min((w-20.)/sw,h/sh)
        page.merge_transformed_page(source,Transformation().scale(scale).translate(x+20.,y+(h-sh*scale)/2))
    labels_path=vectors/f"{name}_labels.pdf";fig=plt.figure(figsize=(width/72,height/72))
    for component,(x,y,w,h) in zip(components,rectangles):fig.text(x/width,(y+h-12)/height,component["id"],weight="bold",fontsize=15)
    label_bounds=S.presentation_bounds(fig);fig.savefig(labels_path,transparent=True);plt.close(fig)
    page.merge_page(PdfReader(labels_path).pages[0])
    with pdf.open("wb") as stream:writer.write(stream)
    subprocess.run(["pdftoppm","-singlefile","-r","160","-png",str(pdf),str(png.with_suffix(""))],check=True,capture_output=True)
    result=dict(passed=True,complete=True,status="candidate",final_export_eligible=False,
        generated_utc=datetime.now(timezone.utc).isoformat(),framing=spec["framing"],no_clean_mapping=spec['no_clean_mapping'],scientific_gates=gates,
        main=dict(pdf=dict(path=str(pdf),sha256=P.sha(pdf),pages=1),png=dict(path=str(png),sha256=P.sha(png))),
        composition=dict(width_points=width,height_points=height,slots=[dict(id=item['id'],rectangle=list(rect)) for item,rect in zip(components,rectangles)],aspect_ratio_preserved=True),
        native_panel=dict(pdf=native,selected_histories=identities,presentation_bounds=bounds,
            selection_CSV=selection_record,native_field_CSV_shards=list(native_files.values()),
            portable_page_metadata=[selection_record,*companions.values()]),
        components=components,panel_label_bounds=label_bounds,
        minimum_printed_width_inches=12.,recommended_printed_width_inches=14.,
        required_legend_notes=list(dict.fromkeys(note for panel in spec['statistical_panels'] for note in panel.get('required_legend_notes',[]))),
        no_native_integrations=True,no_framing_or_native_selection_by_renderer=True,
        source_spec=dict(path=str(spec_path),sha256=P.sha(spec_path)),renderer_sha256=P.sha(__file__),
        shared_source_sha256={str(Path(module.__file__).resolve()):P.sha(module.__file__) for module in (Atlas,P,Q,S)})
    P.write_json(metadata/f"{name}_display_audit.json",result)
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--spec",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True);parser.add_argument("--name",required=True)
    args=parser.parse_args();print(json.dumps(render(args.spec,args.output,args.name)["main"]))
