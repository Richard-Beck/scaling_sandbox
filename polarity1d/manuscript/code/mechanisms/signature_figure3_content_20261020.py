"""Outline-ready scientific Figure3 content from the final audited exports."""
from pathlib import Path
from datetime import datetime,timezone
import json
import hashlib
import pandas as pd

REPO=Path(__file__).resolve().parents[4]
TOPIC=Path('polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_program_20261020/figure3_distance')
OUT=REPO/TOPIC
NAMES=dict(wave_pinning='Wave-pinning',goryachev='Goryachev',holmes_model3='Holmes3',spring='Spring')


def main():
    A1=pd.read_csv(OUT/'supporting/data/A1_complete_region_evidence.csv')
    B=pd.read_csv(OUT/'supporting/data/B_complete_region_evidence_and_geometry.csv')
    heading=pd.read_csv(OUT/'supporting/data/A1_direct_heading_check.csv')
    evidence=[]
    mainA=A1[A1.control.isna()]
    for (model,size),part in mainA[(mainA.classification=='external')&mainA.handover_time.notna()].groupby(['model','size_parameter']):
        part=part.sort_values('gap');first,last=part.iloc[0],part.iloc[-1]
        evidence.append(dict(model=model,size_label=float(size),first_valid_external_gap=float(first.gap),
            first_valid_handover_time=float(first.handover_time),antipode_gap=float(last.gap),antipode_handover_time=float(last.handover_time)))
    angle={model:float(part.handover_angle_fraction_of_requested.median()) for model,part in
        heading[heading.control.isna()&(heading.classification=='external')&(heading.coherence_at_handover>=.05)].groupby('model')}
    content=dict(id='figure3',title='Figure 3. Separation affects native competition differently across models',
        role='Test whether greater incumbent–challenger separation slows competition at fixed circumference, then examine its geometric translation into measured activity direction.',
        purpose='Hold each native basal and cell size fixed while rotating one native cue. Separate cue-driven takeover, cue-free two-region resolution and actual global-heading reorientation.',
        status='final',main_pages=2,supplement_pages=7,
        assets=[dict(path=str(TOPIC/'figure3_distance_main.png'),caption='Fixed-size cue placement, native histories and separation responses',page_size=[10.5,13.5],minimum_printed_width_inches=9.5,page_title='Figure 3. Separation and native competition'),
            dict(path=str(TOPIC/'figure3_distance_translation.png'),caption='Measured-width geometry remaps and direct native-heading qualifications',page_size=[10.5,10.5],minimum_printed_width_inches=9.5,page_title='Figure 3. Geometry and measured activity direction')],
        main_pdf=str(TOPIC/'supporting/vectors/figure3_distance_main.pdf'),
        supplement_pdf=str(TOPIC/'supporting/vectors/figure3_distance_supplement.pdf'),
        reader=str(TOPIC/'index.html'),
        legend=[
            'A, A native cue is rotated at fixed size. The primary initial gap g runs from the incumbent initial 10%-contrast edge to the cue support near edge along the shorter route. Chemical circumferences are30,37.8 and47.6; Spring resting-size labels are2.5,3.15 and3.97, with actual circumference recorded. All12 fresh single-region basals qualify without retuning. External placements alone enter the distance claim; borderline negative gaps use open symbols and internal references are supplementary.',
            'B, Native Goryachev and Holmes3 active-field histories at circumference37.8 and external gaps2 and8. The dashed white line marks independently replayed complete-region handover. Holmes3 uses active Rac. Chemical A1 histories complete5,000 or20,000 native units according to the retained continuation rule; Spring completes2,000. Every bound uses its actual completed horizon.',
            'C1, Continuous-cue handover is the first time challenger excess active amount exceeds incumbent excess while two resolved, spatially associated regions exist. C2, Cue-free resolution follows seeding to complete-region challenger/incumbent ratio1; all0.5/1/2 results are supplementary. Its axis uses the measured current global10%-contrast gap at withdrawal. Triangles retain completed finite-cap time lower bounds. Unseeded, overlapping-contour, merged and single-region outcomes remain explicit. Lines do not bridge unavailable, censored or topology-changing placements.',
            'D, Requested center angle combines the physical edge gap, half the measured incumbent width and half the cue support, divided by circumference and expressed in radians. The relation is typeset in the main panel. A fixed edge gap is not automatically a fixed center displacement when native incumbent width changes with size. Current B widths, actual circumference, short/long routes and measured center displacement are retained.',
            'E, Solid curves show measured handover time for requested cue angle; dashed curves evaluate the reference-size handover curve at larger-cell measured-width gaps. The geometry-only remap stops within observed reference support and does not interpolate through missing/censored neighbors. The inset directly measures actual global-heading angle at handover, typically about half the requested angle. Handover time is not time to complete the requested angular turn.',
            'F, The fixed-physical-scale view divides measured region-center angular step by measured cue-free resolution time, with censoring propagated. This conditional potential rate is distinct from the directly observed vector-heading change. Full-angle target-time brackets,50/90/99% sensitivities, final headings, coherence coverage and cancellation gaps are supplementary.',
            'Complete active amounts integrate native active concentration above a common instantaneous minimum over the two full valley-bounded basins, including signal outside the cue. Continuous boundary position uses a local three-point parabolic minimum; actual field amounts use exact piecewise-linear fractional-cell integration. Active flat/nonunique saddles are unscorable. Single moving regions and merged fields do not supply invented two-region handovers or competition times.',
            'The global vector is the first spatial harmonic divided by instantaneous mean active field; R = |M₁|/M₀ is coherence, unlike Figure4 basal-M1-normalized magnitude. Raw exact-direction availability uses R > 10⁻¹⁰, with displayed R ≥ 0.05 and0.1/0.2 sensitivities. Direction-quality/cancellation gaps break angle unwrapping. All1,150 histories include saved-interval principal angular-rate time courses; these are not instantaneous RHS rates or counts of unobserved complete turns.',
            'A2 uses0.5,0.75,1,1.5 and2 times each placement own corrected handover, then13,000 cue-free chemical units or1,300 Spring units. B withdraws at first valid ratios0.5,1 or2, then20,000 chemical units or2,000 Spring units. Native finite-field stationarity is reported independently of scalar/spatial fate. No finite cap substitutes for an unavailable own handover.',
            'All native positivity and conserved-pool audits pass. All24 B mirrors agree. A1 handover and B seed-time numerical differences are at most0.0449%;35/36 B resolution brackets overlap. Goryachev C37.8, gap4, ratio1 shifts from1040–1060 to1080–1100 units on the finer mesh, about4% in resolution timing. Native clocks differ; no absolute cross-model time ranking or cell-velocity law is imposed.'
        ],
        narrative=[
            'Under continuous forcing, Goryachev handover increases with separation, while Holmes3 is nearly flat over valid two-region placements. Wave-pinning and Spring show near-edge increases followed by broad plateaux. These are model-dependent responses under the declared cue, not a universal absence of distance effects in Wave-pinning.',
            'Cue-free competition has a different distance response: Wave-pinning, Holmes3 and Spring can coarsen very slowly, with many large-size resolution lower bounds. Goryachev includes a nonmonotone near-edge response at ratio1; the full curve and every ratio are retained.',
            'Actual global-heading measurements qualify the proposed geometry-to-persistence inference. Complete-region handover typically occurs at about half the requested cue angle. Cue-free fields can resolve to an intermediate heading, so a potential center step is not always the realized angular event. These data do not establish an unconditional geometry-driven persistence gain.'
        ],
        measurements=[
            'At fixed native size, vary only cue placement and measure complete-region handover, incumbent50% excess loss and challenger reaching50% of basal incumbent excess.',
            'Seed two complete regions to a declared excess ratio, withdraw the cue and measure their actual10%-contrast separation, first single-peak resolution time, winner association and exact cue-free RHS exchange budget.',
            'Measure the same-run global activity vector, its coherence, actual angular displacement, full-target time brackets and saved-interval angular rates; compare them with geometric remaps.'
        ],
        limitations=[
            'Basal width, center displacement and edge gap differ; fixed edge gap is not automatically fixed physical center step across sizes.',
            'A native activity-vector heading is not observed migration direction, a cell velocity law or a conditional OU persistence estimate.',
            'Finite caps produce lower bounds, not proof of permanent arrest. Scalar/spatial winner and full native-field stationarity remain separate.',
            'Two resolved peaks with merged primary10% contours are ineligible for the primary B gap claim. One moving/merged region is a separate outcome.',
            'The finer-mesh Goryachev B ratio1 example has about4% resolution-time sensitivity despite seed-time and winner agreement.',
            'The native-pattern responses do not uniquely identify a molecular mechanism or provide causal allocation to one kinetic term. No noise, moving cue, dose sweep or deBelly ring enters this figure.'
        ],
        claims=['Separation response depends on the model and on whether competition is externally driven or cue-free.',
            'Goryachev has a stronger continuous-cue handover distance response than Holmes3 over valid external two-region placements.',
            'Measured global heading qualifies the proposed geometry-only persistence translation; complete-region handover is typically near half the requested angle.'],
        mechanistic_scope='Fixed-size rotation isolates separation from changes in size, incumbent extent and pool size. The observed responses constrain mechanism classes, but are not a kinetic-term ablation, molecular-transport measurement or unique causal identification.',
        actual_horizons=dict(chemical_A1_initial=5000,chemical_A1_extended_cap=20000,spring_A1=2000,
            counts_including_A1_mirrors=[dict(model=model,native_horizon=float(horizon),cases=int(count)) for (model,horizon),count in A1.groupby(['model','original_native_horizon']).size().items()],
            original_A1_continuation_rule='Chemical5,000-unit histories extend to20,000 unless the retained fixed-footprint incumbent excess is at most10% and a cue-associated single threshold region exists. This original continuation gate does not define the corrected full-region handover.'),
        full_region_definition=dict(signal='Native active species; active Rac for Holmes3',floor='Common instantaneous profile minimum',
            region='Two complete periodic watersheds bounded by intervening subcell parabolic saddles',amount='Piecewise-linear integral of concentration minus common floor, including fractional boundary cells',
            handover='First continuous EC>EI crossing in genuine two-region topology; D=0 crossing, distinct from D ≥ 0.5 replacement',
            ambiguities='Active flat/nonunique saddles, merged/single-region movement and unresolved spatial identity remain explicit unavailable/alternate outcomes'),
        direct_heading_check=dict(median_handover_angle_fraction=angle,coherence='Instantaneous |M₁|/M₀',
            exact_direction_floor=1e-10,reporter_floor=.05,sensitivity_floors=[.1,.2],target_fractions=[1.,.5,.9,.99]),
        observed_A1_endpoints=evidence,
        coverage=dict(fresh_basals=12,main_A1_placements=152,A1_mirrors=24,valid_handover_refinements=133,
            unavailable_external_or_borderline_handover=34,A2_native_runs=545,A2_shortest_tested_switch_brackets=72,A2_no_observed_switch_at_tested_durations=37,
            primary_B_runs=357,B_mirrors=24,B_total_eligible_withdrawals=307,B_cue_phase_single_region_takeovers=39,
            B_ineligible_primary_contours=28,B_unreached_seed_caps=7,B_resolution_lower_bounds=52,native_history_images=1150),
        controls=dict(no_cue_duration_chemical=53000,no_cue_duration_spring=5300,fine_mesh_basals=2,A1_numeric_runs=12,A1_control_roots=10,
            B_numeric_runs=36,B_status_eligibility_winner_agreement='36/36',B_resolution_bracket_overlap='35/36',
            B_fine_mesh_exception='Goryachev C37.8, initial gap4, ratio1:1040–1060 versus1080–1100'),
        sources=[dict(label=label,path=str(TOPIC/path)) for label,path in [
            ('Complete corrected A1 evidence','supporting/data/A1_complete_region_evidence.csv'),
            ('All A2 switching brackets and unavailable calibrations','supporting/data/A2_all_switching_brackets.csv'),
            ('All B gap, fate, censoring and actual-heading measurements','supporting/data/B_complete_region_evidence_and_geometry.csv'),
            ('Support-limited measured geometry remaps','supporting/data/A1_fixed_angle_geometry_remap.csv'),
            ('Direct global-heading target checks','supporting/data/A1_heading_target_time_brackets.csv'),
            ('B numerical convergence including the material exception','supporting/data/B_control_convergence.csv'),
            ('Final artifacts and source-link audit','supporting/metadata/final_artifact_audit.json'),
            ('Lossless full native-history pixel audit','supporting/metadata/lossless_pixel_audit.json')]],
        provenance=dict(generated_utc=datetime.now(timezone.utc).isoformat(),generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_csv_sha256={name:hashlib.sha256((OUT/'supporting/data'/name).read_bytes()).hexdigest() for name in ['A1_complete_region_evidence.csv','B_complete_region_evidence_and_geometry.csv','A1_direct_heading_check.csv']}))
    # Readable prose spacing while preserving paths/measurement identifiers.
    import re
    for key in ['legend','narrative','measurements','limitations','claims']:
        content[key]=[re.sub(r'(?<=[a-z])(?=\d)',' ',re.sub(r'(?<=\d),(?=\d)',', ',text)) for text in content[key]]
    path=OUT/'supporting/metadata/content.json';path.write_text(json.dumps(content,indent=2,allow_nan=False)+'\n')
    audit=dict(content_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),**content['provenance'])
    (OUT/'supporting/metadata/content_generation.json').write_text(json.dumps(audit,indent=2)+'\n');print(path)


if __name__=='__main__':main()
