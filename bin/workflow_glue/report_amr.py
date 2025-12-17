"""Functions for creation of AMR report section."""
from dominate import tags as html_tags
from workflow_glue.report_utils import capitalize_name, fill_none


# Special AMR classes that should appear at the end of list
SPECIAL_AMR_CLASSES = {
    'unspecified', 'under_development', 'other', 'unknown',
    'various', 'various antifungals'
}


def sort_amr_classes(class_name):
    """Sort AMR classes: alphabetically, with special classes at the end."""
    class_lower = class_name.lower()
    if class_lower in SPECIAL_AMR_CLASSES:
        return (1, class_lower)  # Special classes at the end
    else:
        return (0, class_lower)  # Regular classes alphabetically


def format_gene_name(gene_name):
    """Format gene name with lowercase first letter."""
    if not gene_name:
        return gene_name
    return gene_name[0].lower() + gene_name[1:]


def has_fungamr_pmid(variation):
    """Check if variation has FungAMR database reference."""
    return variation.pmids and any(
        "fungamr-db" in str(pmid).lower()
        for pmid in variation.pmids
    )


def get_badge_style(status):
    """Get badge background color style for a given status."""
    class_map = {
        'detected_with_resistance': 'bg-primary',
        'detected_no_resistance': 'bg-secondary'
    }
    return class_map.get(status, None)


def _is_pointfinder(item, databases):
    """Check if a region or variation is from PointFinder database."""
    if not (item.ref_database and databases):
        return False
    ref_dbs = (
        [item.ref_database] if isinstance(item.ref_database, str)
        else item.ref_database
    )
    db_keys = {db.key for db in databases if db.database_name.lower() == 'pointfinder'}
    return any(db_key in db_keys for db_key in ref_dbs)


def _add_genes_from_regions(
    phenotype, seq_regions, databases, genes_with_resistance
):
    """Add genes with resistance from seq_regions to the set."""
    if not (phenotype.seq_regions and seq_regions):
        return
    for region in seq_regions:
        if region.key in phenotype.seq_regions:
            if _is_pointfinder(region, databases):
                genes_with_resistance.add(region.get_gene_name())


def _add_parent_genes(variation, seq_regions, genes_with_resistance):
    """Add parent genes for a variation to the set."""
    if not (variation.seq_regions and seq_regions):
        return
    for region in seq_regions:
        if region.key in variation.seq_regions:
            genes_with_resistance.add(region.get_gene_name())


def _add_genes_from_variations(
    phenotype, seq_variations, seq_regions, databases, genes_with_resistance
):
    """Add genes with resistance from seq_variations to the set."""
    if not (phenotype.seq_variations and seq_variations):
        return
    for variation in seq_variations:
        if variation.key in phenotype.seq_variations:
            if _is_pointfinder(variation, databases):
                _add_parent_genes(variation, seq_regions, genes_with_resistance)


def get_pointfinder_genes_with_resistance(amr_data):
    """Get all PointFinder genes found and identify which have resistance phenotypes."""
    pointfinder_genes = {}
    genes_with_resistance = set()
    if not amr_data.phenotypes:
        return pointfinder_genes
    # Collect genes with resistance phenotypes
    for phenotype in amr_data.phenotypes:
        _add_genes_from_regions(
            phenotype, amr_data.seq_regions, amr_data.databases, genes_with_resistance
        )
        _add_genes_from_variations(
            phenotype, amr_data.seq_variations, amr_data.seq_regions,
            amr_data.databases, genes_with_resistance
        )
    # Collect detected PointFinder genes and mark resistance status
    detected_genes = set()
    if amr_data.seq_regions:
        for region in amr_data.seq_regions:
            if _is_pointfinder(region, amr_data.databases):
                gene_name = region.get_gene_name()
                detected_genes.add(gene_name)
                has_resistance = gene_name in genes_with_resistance
                pointfinder_genes[gene_name] = (
                    'detected_with_resistance' if has_resistance
                    else 'detected_no_resistance'
                )
    # Add genes from database that were not detected
    if amr_data.pointfinder_available_genes:
        for gene in amr_data.pointfinder_available_genes:
            if gene not in detected_genes:
                pointfinder_genes[gene] = 'not_detected'
    return pointfinder_genes


def create_badge_header(resistance_by_class, pointfinder_genes):
    """Header with antimicrobial classes and PointFinder gene badges."""
    with html_tags.div(cls="row mb-3 mt-0"):
        # Left column - Antibiotic classes
        with html_tags.div(cls="col-md-6"):
            html_tags.p(
                html_tags.b(
                    "Resistance to the following antimicrobial and "
                    "disinfectant classes was found:"
                ),
                cls="mb-2"
            )
            with html_tags.div(cls="mb-3"):
                for amr_class in sorted(
                    resistance_by_class.keys(), key=sort_amr_classes
                ):
                    html_tags.span(
                        capitalize_name(amr_class),
                        cls="badge bg-primary me-2 mb-2"
                    )
        # Right column - PointFinder genes
        if pointfinder_genes:
            with html_tags.div(cls="col-md-6"):
                html_tags.p(
                    html_tags.b("Genes analysed for AMR point mutations:"),
                    cls="mb-2"
                )
                with html_tags.div(cls="mb-3"):
                    for gene_name in sorted(
                        pointfinder_genes.keys(), key=lambda x: x.lower()
                    ):
                        status = pointfinder_genes[gene_name]
                        badge_class = get_badge_style(status)
                        display_gene_name = format_gene_name(gene_name)
                        if badge_class:
                            html_tags.span(
                                display_gene_name,
                                cls=f"badge {badge_class} me-2 mb-2"
                            )
                        else:
                            html_tags.span(
                                display_gene_name,
                                cls="badge me-2 mb-2",
                                style="background-color: #a0a0a0; color: white;"
                            )
                # Add legend
                with html_tags.div(cls="mt-2"):
                    html_tags.small(
                        html_tags.span(
                            "",
                            cls="badge bg-primary me-1",
                            style=(
                                "display: inline-block; width: 10px; height: 10px; "
                                "vertical-align: middle; padding: 0;"
                            )
                        ),
                        "Known AMR mutations",
                        html_tags.span(
                            "",
                            cls="badge bg-secondary me-1 ms-3",
                            style=(
                                "display: inline-block; width: 10px; height: 10px; "
                                "vertical-align: middle; padding: 0;"
                            )
                        ),
                        "No known AMR mutations",
                        html_tags.span(
                            "",
                            cls="badge me-1 ms-3",
                            style=(
                                "display: inline-block; width: 10px; height: 10px; "
                                "vertical-align: middle; padding: 0; "
                                "background-color: #a0a0a0; color: white;"
                            )
                        ),
                        "Region of interest not detected",
                        cls="text-muted"
                    )


def _create_database_badges(database_counts):
    """Create database badges with proper capitalization."""
    if not database_counts:
        return []
    db_name_map = {
        'pointfinder': 'PointFinder',
        'resfinder': 'ResFinder',
        'disinfinder': 'DisinFinder'
    }
    db_badges = []
    for db_name in sorted(database_counts.keys()):
        db_count = database_counts[db_name]
        if db_count == 0:
            continue
        db_display = db_name_map.get(db_name.lower(), db_name.capitalize())
        db_badges.append(
            html_tags.span(
                db_display,
                html_tags.span(
                    str(db_count),
                    cls="position-absolute badge rounded-pill bg-danger",
                    style="top:-8px; right:-8px;"
                ),
                cls="badge bg-primary position-relative me-2"
            )
        )
    return db_badges


def create_accordion_button(amr_class, class_info, row):
    """Create accordion button wit antibiotic class name and badges."""
    db_badges = _create_database_badges(class_info['database_counts'])
    _button = html_tags.button(
        html_tags.span(
            html_tags.b(capitalize_name(amr_class)),
            style="flex: 0 0 auto; max-width: 230px;"
        ),
        html_tags.span(
            f"{len(class_info['unique_drugs'])} drug(s)",
            cls="badge bg-secondary",
            style="position: absolute; left: 250px;"
        ),
        html_tags.span(style="flex: 1 1 auto;"),  # Spacer
        html_tags.span(
            *db_badges,
            style="flex: 0 0 auto; white-space: nowrap; margin-right: 2rem;"
        ),
        cls="accordion-button collapsed",
        type="button",
        data_bs_toggle="collapse",
        data_bs_target=f"#collapse{row}",
        aria_expanded="false",
        aria_controls=f"collapse{row}",
        style="display: flex; align-items: center;"
    )
    return _button


def _get_region_databases(region, phenotypes_for_class, databases):
    """Extract database sources for a region."""
    db_sources = set()
    for phenotype in phenotypes_for_class:
        if not (phenotype.seq_regions and region.key in phenotype.seq_regions):
            continue
        if not phenotype.ref_database:
            continue
        for db_key in phenotype.ref_database:
            for database in databases:
                if database.key == db_key:
                    db_sources.add(database.database_name)
    return db_sources


def _get_variation_databases(variation, phenotypes_for_class, databases):
    """Extract database sources for a variant."""
    db_sources = set()
    for phenotype in phenotypes_for_class:
        if not (phenotype.seq_variations and variation.key in phenotype.seq_variations):
            continue
        for database in databases:
            if database.key == variation.ref_database:
                db_sources.add(database.database_name)
    return db_sources


def _get_database_sources(region, variations, databases, phenotypes_for_class):
    """Get database sources relevant to this AMR class only."""
    db_sources = set()
    if not (databases and phenotypes_for_class):
        return db_sources
    if region.phenotypes:
        db_sources.update(
            _get_region_databases(region, phenotypes_for_class, databases))
    for var_data in variations.values():
        variation = var_data['variation']
        if variation.ref_database:
            db_sources.update(
                _get_variation_databases(variation, phenotypes_for_class, databases))
    return db_sources


def _create_gene_region_row(
    gene_name, gene_info, has_mutations, row_id, amr_data, fill_na,
    visual_row_num, phenotypes_for_class
):
    """Create a single gene region table row."""
    region = gene_info['region']
    # "Gene" cell content
    display_gene_name = (
        gene_name[0].lower() + gene_name[1:] if gene_name
        else gene_name
    )
    gene_cell_content = display_gene_name
    if has_mutations:
        gene_cell_content = html_tags.span(
            "▸ ",
            display_gene_name,
            html_tags.span(
                f" ({len(gene_info['variations'])})",
                cls="text-muted small"
            )
        )
    # List of antimicrobials
    drugs_str = ", ".join(
        sorted([capitalize_name(d) for d in gene_info['drugs']])
    ) if gene_info['drugs'] else fill_na
    # Databases cell - only show databases relevant to this abx class
    db_sources = _get_database_sources(
        region, gene_info['variations'], amr_data.databases, phenotypes_for_class)
    db_sources_str = ", ".join(sorted(db_sources)) if db_sources else fill_na

    # Determine background color for striping
    bg_color = "#f2f2f2" if (visual_row_num % 2) == 1 else "#ffffff"
    _tr = html_tags.tr(
        data_bs_toggle="collapse" if has_mutations else None,
        data_bs_target=f"#{row_id}" if has_mutations else None,
        style=(
            f"cursor: pointer; background-color: {bg_color}; "  # noqa: E702
            if has_mutations
            else f"background-color: {bg_color}; "  # noqa: E702
        )
    )
    _tr.add(
        html_tags.td(gene_cell_content),
        html_tags.td(drugs_str),
        html_tags.td(fill_none(region.query_contig, fill_na)),
        html_tags.td(
            f"{fill_none(region.query_start_pos, fill_na)} - "
            f"{fill_none(region.query_end_pos, fill_na)}"
        ),
        html_tags.td(
            f"{region.identity: .2f}%" if region.identity is not None else fill_na
        ),
        html_tags.td(
            f"{region.coverage: .2f}%" if region.coverage is not None else fill_na
        ),
        html_tags.td(fill_none(region.ref_acc, fill_na)),
        html_tags.td(", ".join(region.pmids) if region.pmids else fill_na),
        html_tags.td(db_sources_str)
    )
    return _tr


def _create_mutation_dropdown(variations, all_regions, row_id, fill_na):
    """Create expandable mutation table for a gene."""
    # Check FungAMR entries are present
    has_fungamr = any(
        has_fungamr_pmid(var_data['variation'])
        for var_data in variations.values()
    )
    _mutation_row = html_tags.tr()
    _mutation_cell = html_tags.td(
        colspan="9",
        cls="p-0",
        style="background-color: #ffffff;"
    )
    _mutation_collapse = html_tags.div(
        id=row_id,
        cls="collapse",
        style="padding: 0rem 1rem;"
        )

    _mut_table = html_tags.table(
        cls="table table-sm mb-0",
        style="background-color: #ffffff;"
    )
    _mut_thead = html_tags.thead(cls="table-light")
    _header_tr = html_tags.tr(
        html_tags.th("Mutation"),
        html_tags.th("Drug(s)"),
        html_tags.th("Position"),
        html_tags.th("Nucleotide change"),
        html_tags.th("Amino acid change"),
        html_tags.th("PMID")
    )
    if has_fungamr:
        _header_tr.add(html_tags.th("Notes"))
    _mut_thead.add(_header_tr)
    _mut_table.add(_mut_thead)
    _mut_tbody = html_tags.tbody()

    # Sort by position
    sorted_variations = sorted(
        variations.values(),
        key=lambda x: x['variation'].ref_start_pos
        if x['variation'].ref_start_pos is not None else 0
    )
    for var_data in sorted_variations:
        variation = var_data['variation']
        drugs = var_data['drugs']
        contig, start, end = variation.get_contig_and_positions(all_regions)
        nuc_change = (
            f"{variation.ref_codon.upper()} → {variation.var_codon.upper()}"
            if variation.ref_codon and variation.var_codon
            else fill_na
        )
        if variation.ref_aa and variation.var_aa:
            aa_change = f"{variation.ref_aa.upper()} → {variation.var_aa.upper()}"
        elif variation.is_promoter_mutation():
            aa_change = "Promoter mutation"
        else:
            aa_change = fill_na
        # Format PMIDs and replace fungamr-db
        pmid_display = fill_na
        if variation.pmids:
            formatted_pmids = []
            for pmid in variation.pmids:
                pmid_str = str(pmid)
                if "fungamr-db" in pmid_str.lower():
                    formatted_pmids.append("FungAMR-DB")
                else:
                    formatted_pmids.append(pmid_str)
            pmid_display = ", ".join(formatted_pmids)
        # Format drugs - sort and capitalize
        drugs_display = fill_na
        if drugs:
            drugs_display = ", ".join(
                sorted([capitalize_name(d) for d in drugs])
            )
        _mut_tr = html_tags.tr(style="background-color: #ffffff;")
        _mut_tr.add(
            html_tags.td(variation.seq_var or fill_na),
            html_tags.td(drugs_display),
            html_tags.td(fill_none(start, fill_na)),
            html_tags.td(nuc_change),
            html_tags.td(aa_change),
            html_tags.td(pmid_display)
        )
        # Add Notes column if FungAMR present
        if has_fungamr:
            notes_display = fill_na
            if variation.notes:
                notes_display = (
                    "; ".join(str(note) for note in variation.notes)
                    if isinstance(variation.notes, list)
                    else str(variation.notes)
                )
            _mut_tr.add(html_tags.td(notes_display))
        _mut_tbody.add(_mut_tr)
    _mut_table.add(_mut_tbody)
    _mutation_collapse.add(_mut_table)
    _mutation_cell.add(_mutation_collapse)
    _mutation_row.add(_mutation_cell)
    return _mutation_row


def create_gene_region_table(
    gene_data, amr_data, row, fill_na, phenotypes_for_class
):
    """Create gene region table with expandable mutation rows."""
    _region_table = html_tags.table(cls="table")
    _thead = html_tags.thead()
    _thead.add(
        html_tags.tr(
            html_tags.th("Gene"),
            html_tags.th("Drug(s)"),
            html_tags.th("Contig"),
            html_tags.th("Gene position"),
            html_tags.th("Identity"),
            html_tags.th("Coverage"),
            html_tags.th("Accession"),
            html_tags.th("PMID"),
            html_tags.th("Database")
        )
    )
    _region_table.add(_thead)
    _tbody = html_tags.tbody()
    for row_num, gene_name in enumerate(
        sorted(gene_data.keys(), key=lambda x: x.lower()), start=1
    ):
        gene_info = gene_data[gene_name]
        has_mutations = len(gene_info['variations']) > 0
        mutation_row_id = f"mutations{row}_{row_num}"
        gene_row = _create_gene_region_row(
            gene_name, gene_info, has_mutations, mutation_row_id,
            amr_data, fill_na, row_num, phenotypes_for_class
        )
        _tbody.add(gene_row)
        # Add mutation dropdown if applicable
        if has_mutations:
            mutation_row = _create_mutation_dropdown(
                gene_info['variations'],
                amr_data.seq_regions or [],
                mutation_row_id,
                fill_na
            )
            _tbody.add(mutation_row)
    _region_table.add(_tbody)
    return _region_table


def antibiotic_section(amr_data, fill_na="-"):
    """Create antibiotic resistance section organized by class."""
    if not amr_data:
        html_tags.b("No AMR data available for this sample.")
        return
    # Get resistance data grouped by class
    resistance_by_class = amr_data.get_resistance_by_class()
    if not resistance_by_class:
        html_tags.b("No antibiotic class data available.")
        return

    # Get PointFinder genes with resistance status
    pointfinder_genes = get_pointfinder_genes_with_resistance(amr_data)
    # Create header with badges
    create_badge_header(resistance_by_class, pointfinder_genes)

    # Check if any FungAMR entries are present
    has_fungamr_in_data = False
    if amr_data.seq_variations:
        for variation in amr_data.seq_variations:
            if variation.pmids and any(
                "fungamr-db" in str(pmid).lower() for pmid in variation.pmids
            ):
                has_fungamr_in_data = True
                break

    # Create accordion table
    with html_tags.div(
        cls="accordion mt-3",
        id="antibioticAccordion",
        style="border-radius: 0.375rem; overflow: hidden;"
    ):
        for row, amr_class in enumerate(
            sorted(resistance_by_class.keys(), key=sort_amr_classes), start=1
        ):
            class_info = resistance_by_class[amr_class]
            gene_data = amr_data.get_genes_by_class(class_info['phenotypes'])
            with html_tags.div(cls="accordion-item"):
                # Header
                _head = html_tags.h2(id=f"heading{row}", cls="accordion-header")
                _button = create_accordion_button(amr_class, class_info, row)
                _head.add(_button)
                # Body
                with html_tags.div(
                    id=f"collapse{row}",
                    cls="accordion-collapse collapse",
                    fr=f"{row}",
                    aria_labelledby=f"heading{row}",
                    data_bs_parent="#antibioticAccordion"
                ):
                    with html_tags.div(cls="accordion-body"):
                        if gene_data:
                            create_gene_region_table(
                                gene_data, amr_data, row, fill_na,
                                class_info['phenotypes']
                            )
    if has_fungamr_in_data:
        html_tags.p(
            html_tags.a(
                "FungAMR",
                href="https://doi.org/10.1038/s41564-025-02084-7"
            ),
            " database entries may indicate antimicrobial resistance or "
            "increased susceptibility.",
            html_tags.br(),
            "Prediction confidence scores shown"
            " in notes: positive scores (1-8) indicate resistance (1="
            "strongest evidence), negative scores (-1 to -8) indicate "
            "increased susceptibility (-1=strongest evidence).",
            html_tags.br(),
            cls="small mt-4"
        )
