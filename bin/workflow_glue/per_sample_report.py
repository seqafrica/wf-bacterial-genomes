"""Create a report for each sample."""
from enum import Enum
import json
import os

import dacite
from dominate import tags as html_tags
from ezcharts.components.reports import labs
import pandas as pd
from workflow_glue.models.custom import Sample
from workflow_glue.report_utils import convert_bp, fill_none

from .collect_results import (  # noqa: ABS101
    gather_sample_files
)
from .util import get_named_logger, wf_parser  # noqa: ABS101


def lead_section(text):
    """From a dict of key value pairs make a text section."""
    _div = html_tags.div()
    for key, value in text.items():
        if key == 'header':
            _div.add(html_tags.h3(value))
        else:
            _div.add(
                html_tags.small(key.replace("_", " ").upper(), cls="text-muted"),
                html_tags.h5(value if value is not None else "None")
            )
    return _div


def flye_section(assembly):
    """Extract information on denovo assembly."""
    summary = assembly.get_flye_summary()

    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(
        html_tags.tr(
            html_tags.th("# Contigs"),
            html_tags.th("Total length"),
            html_tags.th("# Circular contigs"),
        )
    )
    _table.add(_thead)
    _tr = html_tags.tr()
    _tr.add(html_tags.td(summary['contig_num']))
    _tr.add(html_tags.td(summary['total_yield']))
    _tr.add(html_tags.td(summary['circular_num']))
    _table.add(_tr)
    _div.add(_table)
    return _div


def ref_section(total):
    """Return coverage information on reference based assembly."""
    depth_df = pd.read_csv(total, sep="\t", names=["ref", "start", "end", "depth"])
    coverage_dict = dict()
    thresholds = [30, 50]
    for ref, depths in depth_df.groupby("ref"):
        ref_dict = dict()
        for t in thresholds:
            total_length = 0
            above_threshold = 0
            for index, row in depths.iterrows():
                total_length += row["end"] - row["start"]
                if row["depth"] >= t:
                    above_threshold += row["end"] - row["start"]
            above_threshold_pct = above_threshold / total_length * 100
            if "total_len" not in ref_dict:
                ref_dict["total_len"] = convert_bp(total_length)
            ref_dict[t] = round(above_threshold_pct, 3)
        coverage_dict[ref] = ref_dict

    # HTML table
    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(
        html_tags.tr(
            html_tags.th("Contig"),
            html_tags.th("Total length"),
            [html_tags.th(f"% Coverage at {t}x") for t in thresholds]
        )
    )
    _table.add(_thead)
    for k, v in coverage_dict.items():
        _tr = html_tags.tr()
        _tr.add(html_tags.td(k))
        _tr.add(html_tags.td(v["total_len"]))
        _tr.add([html_tags.td(v[t]) for t in thresholds])
        _table.add(_tr)
    _div.add(_table)
    return _div


def amr_section(acquired_data, point_data, html_id):
    """Parse resfinder JSON for accordion style table."""
    if (not acquired_data and not point_data):
        return html_tags.b("No AMR genes detected in sample.")
    # Check if any point mutations have FungAMR entries and add note
    has_fungamr = False
    if point_data:
        for gene, evidence in point_data.items():
            if (
                evidence and
                any(mutation.get("has_fungamr_notes", False) for mutation in evidence)
            ):
                has_fungamr = True
                break
    _container = html_tags.div()
    if has_fungamr:
        _note = html_tags.p(
            html_tags.a(
                "FungAMR",
                href="https://doi.org/10.1038/s41564-025-02084-7"),
            " database entries may indicate antimicrobial resistance or "
            "increased susceptibility. Prediction confidence scores shown"
            " in notes: positive scores (1-8) indicate resistance (1="
            "strongest evidence), negative scores (-1 to -8) indicate "
            "increased susceptibility (-1=strongest evidence).",
            cls="small"
        )
        _container.add(_note)

    _div = html_tags.div(cls="accordion-item")
    # Point mutations
    row = 0
    for gene, evidence in point_data.items():
        if not evidence:
            continue
        row += 1
        show_notes_column = any(
            mutation.get("has_fungamr_notes", False) for mutation in evidence)
        drugs = {drug.capitalize() for mut in evidence for drug in mut["drugs"]}
        _head = html_tags.h2(id=f"{row}", style="border: 1px solid rgba(0,0,0,.125);\
                            border-collapse: collapse;\
                            padding:0;\
                            margin-bottom:0")
        _button = html_tags.button(
            html_tags.span(html_tags.b(gene)),
            html_tags.span(
                [html_tags.span(d, cls="badge bg-dark me-1") for d in drugs]
            ),
            html_tags.span(
                "PointFinder",
                html_tags.span(
                    len(evidence),
                    cls="position-absolute badge rounded-pill bg-danger",
                    style="top:8px"
                    ),
                cls="badge bg-primary"),
            cls="accordion-button collapsed",
            type="button",
            data_bs_toggle="collapse",
            data_bs_target=f"#collapse{row}",
            aria_expanded="false",
            aria_controls=f"collapse{row}",
            style="display: grid; \
                    align-items: center;\
                    grid-template-columns: 150px 1fr max-content max-content;\
                    grid-gap: 25px"
            )
        _head.add(_button)
        _div.add(_head)
        _div1 = html_tags.div(
            id=f"collapse{row}",
            cls="accordion-collapse collapse",
            fr=f"{row}",
            aria_labelledby=f"{row}",
            data_bs_parent=f"#{html_id}")
        _div2 = html_tags.div(cls="accordion body")
        _table = html_tags.table(cls="table table-striped")
        _thead = html_tags.thead()
        header_cells = [
            html_tags.th("Antimicrobial resistance"),
            html_tags.th("Amino Acid"),
            html_tags.th("Nucleotide"),
            html_tags.th("PMID/ Database")
        ]
        if show_notes_column:
            header_cells.append(html_tags.th("Notes"))
        _thead.add(html_tags.tr(header_cells))
        _table.add(_thead)
        for mutation in evidence:
            _tr = html_tags.tr()
            row_cells = [
                html_tags.td("; ".join(d.capitalize() for d in mutation["drugs"])),
                html_tags.td(mutation["aa"]),
                html_tags.td(mutation["nuc"].upper()),
                html_tags.td(mutation["pmids"])
            ]
            if show_notes_column:
                row_cells.append(html_tags.td(mutation.get("notes", "-")))
            _tr.add(row_cells)
            _table.add(_tr)
        _div2.add(_table)
        _div1.add(_div2)
        _div.add(_div1)

    # Acquired resistance
    for gene, evidence in acquired_data.items():
        row += 1
        _head = html_tags.h2(id="f{row}", style="border: 1px solid rgba(0,0,0,.125);\
                            border-collapse: collapse;\
                            padding:0;\
                            margin-bottom:0")
        _button = html_tags.button(
            html_tags.span(html_tags.b(gene)),
            html_tags.span(
                {html_tags.span(
                    d.capitalize(),
                    cls="badge bg-dark me-1"
                    ) for d in evidence["drugs"]}
                ),
            html_tags.span(
                "ResFinder",
                html_tags.span(
                    1,  # TODO add evidence as list for multiple hits of gene
                    cls="position-absolute badge rounded-pill bg-danger",
                    style="top:8px"
                    ),
                cls="badge bg-secondary"),
            cls="accordion-button collapsed",
            type="button",
            data_bs_toggle="collapse",
            data_bs_target=f"#collapse{row}",
            aria_expanded="false",
            aria_controls=f"collapse{row}",
            style="display: grid;\
                grid-template-columns: 150px 1fr max-content max-content;\
                align-items: center;\
                grid-gap: 25px")
        _head.add(_button)
        _div.add(_head)
        _div1 = html_tags.div(
            id=f"collapse{row}",
            cls="accordion-collapse collapse",
            fr=f"{row}",
            aria_labelledby=f"{row}",
            data_bs_parent=f"#{html_id}"
            )
        _div2 = html_tags.div(cls="accordion body")
        _table = html_tags.table(cls="table table-striped")
        _thead = html_tags.thead()
        _thead.add(
            html_tags.tr(
                html_tags.th("Antimicrobial resistance"),
                html_tags.th("Identity"),
                html_tags.th("Coverage")
            )
        )
        _table.add(_thead)
        _tr = html_tags.tr()
        _tr.add(
            html_tags.td("; ".join(d.capitalize() for d in evidence["drugs"])),
            html_tags.td(evidence["identity"]),
            html_tags.td(evidence["coverage"])
        )
        _table.add(_tr)
        _div2.add(_table)
        _div1.add(_div2)
        _div.add(_div1)
    return _div


def mlst_section(mlst_data):
    """Extract mlst results."""
    # Build column names and row data
    col_names = ["Id", "Scheme", "Sequence type"]
    row_data = ["-", mlst_data.detected_species, mlst_data.sequence_type]
    # Add allele columns
    for schema in mlst_data.typing_schema:
        col_names.append(schema.schema_identifier)
        row_data.append(schema.allele_variant)
    # HTML section
    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(html_tags.tr([html_tags.th(c) for c in col_names]))
    _table.add(_thead)
    _tr = html_tags.tr()
    for cell in row_data:
        _tr.add(html_tags.td(fill_none(cell, fill_na='-')))
    _table.add(_tr)
    _div.add(_table)
    return _div


def species_id_section(
        species_identification, sourmash_excluded, fill_na="-"
):
    """Extract species identification from Sourmash results."""
    if not species_identification.detected_matches:
        return html_tags.div(html_tags.p(
            "Species identification was unable to identify the sample. "
            "Please check coverage of sample."
        ))
    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(html_tags.tr([
        html_tags.th("Match rank"),
        html_tags.th("Species call"),
        html_tags.th("Query containment ANI"),
        html_tags.th("Intersect (bp)"),
        html_tags.th("Remaining (bp)"),
        html_tags.th("Fraction original query"),
        html_tags.th("Fraction unique"),
        html_tags.th("Reference Genome"),
    ]))
    _table.add(_thead)
    for match in species_identification.detected_matches:
        _tr = html_tags.tr()
        _tr.add(html_tags.td(match.rank))
        _tr.add(html_tags.td(fill_none(match.species_call, fill_na)))
        ani_value = (
            f"{match.query_containment_ani: .3f}"
            if match.query_containment_ani is not None else fill_na
        )
        _tr.add(html_tags.td(ani_value))
        intersect_value = (
            f"{match.intersect_bp}"
            if match.intersect_bp is not None else fill_na
        )
        _tr.add(html_tags.td(intersect_value))
        remaining_value = (
            f"{match.remaining_bp}"
            if match.remaining_bp is not None else fill_na
        )
        _tr.add(html_tags.td(remaining_value))
        f_orig_value = (
            f"{match.f_orig_query: .3f}"
            if match.f_orig_query is not None else fill_na
        )
        _tr.add(html_tags.td(f_orig_value))
        f_unique_value = (
            f"{match.f_unique_to_query: .3f}"
            if match.f_unique_to_query is not None else fill_na
        )
        _tr.add(html_tags.td(f_unique_value))
        _tr.add(html_tags.td(match.reference_name or fill_na))
        _table.add(_tr)
    _div.add(_table)

    if sourmash_excluded and sourmash_excluded.excluded_assemblies:
        _div.add(html_tags.small(
            html_tags.strong(
                f"Assemblies excluded from default Sourmash database "
                f"({sourmash_excluded.total_excluded}): "
            ),
            ", ".join(sourmash_excluded.excluded_assemblies),
        ))
    else:
        _div.add(html_tags.small(
            "No assemblies were excluded from the Sourmash database.",
        ))

    return _div


def serotype_section(serotype_data):
    """Extract serotyping results."""
    columns = [
        "Predicted serotype",
        "Predicted antigenic profile",
        "Predicted identification",
        "O antigen prediction",
        "H1 antigen prediction(fliC)",
        "H2 antigen prediction(fljB)",
        "Note"
    ]
    row_data = [
        serotype_data.predicted_serotype,
        serotype_data.predicted_antigenic_profile,
        serotype_data.predicted_identification,
        serotype_data.o_antigen_prediction,
        serotype_data.h1_antigen_prediction,
        serotype_data.h2_antigen_prediction,
        serotype_data.note
    ]
    _div = html_tags.div()
    _table = html_tags.table(cls="table table-striped")
    _thead = html_tags.thead()
    _thead.add(html_tags.tr([html_tags.th(c) for c in columns]))
    _table.add(_thead)
    _tr = html_tags.tr()
    for cell in row_data:
        _tr.add(html_tags.td(fill_none(cell, fill_na='-')))
    _table.add(_tr)
    _div.add(_table)
    return _div


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")

    sample_json = os.path.join(args.data_dir, f"{args.sample_alias}.json")
    with open(sample_json) as f:
        sample_dict = json.load(f)
        sample = dacite.from_dict(
            data_class=Sample,
            data=sample_dict,
            config=dacite.Config(cast=[Enum])
        )

    report = labs.LabsReport(
        f"{args.sample_alias} | Isolate Sequencing Report",
        "wf-bacterial-genomes",
        args.params,
        args.versions,
        args.wf_version
    )

    files = gather_sample_files(args.sample_alias, args.data_dir)

    lead_summary = lead_section(dict(
        sample=args.sample_alias,
        barcode=args.sample_barcode,
        run=args.wf_session,
        version=args.wf_version,
    ))
    with open(args.params) as fh:
        params_data = json.load(fh)

    run_summary_dict = sample.get_run_summary(
        reference=params_data.get("reference")
    )

    run_summary = lead_section(run_summary_dict)

    with report.add_section('Workflow details', 'Workflow', True):
        with html_tags.div(cls="row"):
            html_tags.div(lead_summary, cls="col-sm-12", style="float: left; width:50%")
            html_tags.div(run_summary, cls="col-sm-12", style="float: right; width:50%")

    with report.add_section("Multilocus sequencing typing (MLST)", "MLST"):
        with html_tags.div(cls="row"):
            if sample.has_mlst():
                mlst_section(sample.results.sequence_typing)
            else:
                html_tags.b(
                    "MLST was unable to identify scheme for "
                    "this sample. Please check coverage of sample."
                )

    with report.add_section("Species identification", "SpeciesID"):
        with html_tags.div(cls="row"):
            if sample.has_species_identification():
                species_data = sample.results.species_identification
                excluded_data = sample.results.sourmash_excluded_genomes
                html_tags.p(
                    "Taxonomic identification was done using Sourmash."
                )
                species_id_section(species_data, excluded_data)
            else:
                html_tags.b(
                    "Sourmash was unable to identify species for this sample."
                )

    if sample.has_serotype():
        with report.add_section("Salmonella serotyping", "Serotype"):
            with html_tags.div(cls="row"):
                serotype_section(sample.results.serotyping)

    with report.add_section('Antimicrobial resistance prediction', 'AMR', True):
        with html_tags.div(cls="accordion", id="accordionTable"):
            if sample.has_amr():
                html_tags.p(
                    "Analysis was performed using ResFinder. "
                    "Click on each gene for more information."
                )
                acquired = (
                    sample.results.antimicrobial_resistance
                    .get_acquired_resistance_display()
                )
                mutations = (
                    sample.results.antimicrobial_resistance
                    .get_point_mutations_display()
                )
                amr_section(acquired, mutations, "accordionTable")
            else:
                html_tags.b("No AMR genes detected for this sample.")
            if params_data.get("pointfinder_ignore_indels", False):
                html_tags.small(
                    "PointFinder analysis was run with --ignore_indels "
                    "flag enabled."
                )
                html_tags.br()
            if params_data.get("pointfinder_ignore_stop_codons", False):
                html_tags.small(
                    "PointFinder analysis was run with --ignore_stop_codons "
                    "flag enabled."
                )

    with report.add_section("Assembly QC", "Assembly", True):
        with html_tags.div(cls="row"):
            if not args.denovo:
                html_tags.p(
                    "Analysis was completed using an alignment with the provided "
                    "reference, and Medaka was used for variant calling."
                )
                if files["depth"]:
                    ref_section(files["depth"])
                else:
                    html_tags.b(
                        "No coverage information, please check input data quality."
                        )
            else:
                html_tags.p(
                    "As no reference was provided the reads were assembled"
                    " and corrected using Flye and Medaka."
                )
                if sample.has_assembly():
                    flye_section(sample.results.assembly)
                else:
                    html_tags.b(
                        "No denovo assembly was produced for this sample. "
                        "Please check input data quality."
                    )

    report.write(args.output)
    logger.info(f"Report written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--denovo",
        action="store_true",
        help="Analysis performed de-novo assembly (instead of variant calling).",
    )
    parser.add_argument(
        "--versions",
        required=True,
        help="directory containing CSVs containing name,version.",
    )
    parser.add_argument(
        "--params",
        default=None,
        required=True,
        help="A JSON file containing the workflow parameter key/values",
    )
    parser.add_argument("--output", help="Report output filename")
    parser.add_argument("--sample_alias", required=True)
    parser.add_argument("--sample_barcode", required=True)
    parser.add_argument("--data_dir", required=True, help="Analysis results directory")
    parser.add_argument("--wf_session", required=True)
    parser.add_argument(
        "--wf_version", default="unknown",
        help="version of the executed workflow")
    return parser
