#!/usr/bin/env python
"""Process resfinder information from output JSON for isolate report sequencing."""
import os

import pandas as pd
from workflow_glue.models.custom import (
    AMRDatabase, AMRPhenotype, AMRSeqRegion, AMRSeqVariation, MutationType
)


def get_mutation_type(variant_data):
    """Determine mutation type from ResFinder boolean flags."""
    if variant_data.get("substitution", False):
        return MutationType.substitution
    elif variant_data.get("deletion", False):
        return MutationType.deletion
    elif variant_data.get("insertion", False):
        return MutationType.insertion
    else:
        return None


def parse_databases(databases_dict):
    """Parse database information from ResFinder JSON."""
    databases = []
    for key, data in databases_dict.items():
        databases.append(AMRDatabase(
            key=key,
            database_name=data["database_name"],
            database_version=data["database_version"]
        ))
    return databases


def parse_seq_regions(seq_regions_dict):
    """Parse sequence regions from ResFinder JSON."""
    seq_regions = []
    for key, data in seq_regions_dict.items():
        seq_regions.append(AMRSeqRegion(
            key=key,
            phenotypes=data["phenotypes"],
            ref_database=data["ref_database"],
            gene_name=data["name"],
            ref_acc=data.get("ref_acc"),
            identity=data["identity"],
            coverage=data["coverage"],
            alignment_length=data["alignment_length"],
            ref_seq_length=data["ref_seq_length"],
            ref_start_pos=data["ref_start_pos"],
            ref_end_pos=data["ref_end_pos"],
            query_contig=data["query_id"],
            query_start_pos=data["query_start_pos"],
            query_end_pos=data["query_end_pos"],
            pmids=data["pmids"],
            notes=data["notes"],
            grade=data["grade"]
        ))
    return seq_regions


def parse_seq_variations(seq_variations_dict):
    """Parse sequence variations from ResFinder JSON."""
    seq_variations = []
    for key, data in seq_variations_dict.items():
        # Skip variations without phenotypes
        # These are point mutations, that are not present in the DB
        if not data["phenotypes"]:
            continue

        seq_variations.append(AMRSeqVariation(
            key=key,
            ref_database=data["ref_database"],
            seq_regions=data["seq_regions"],
            phenotypes=data["phenotypes"],
            seq_var=data["seq_var"],
            ref_codon=data["ref_codon"],
            var_codon=data["var_codon"],
            ref_aa=data.get("ref_aa"),
            var_aa=data.get("var_aa"),
            codon_change=data.get("codon_change"),
            nuc_change=data.get("nuc_change"),
            ref_start_pos=data["ref_start_pos"],
            ref_end_pos=data["ref_end_pos"],
            mut_type=get_mutation_type(data),
            pmids=data["pmids"],
            notes=data["notes"]
        ))
    return seq_variations


def parse_phenotypes(phenotypes_dict):
    """Parse phenotype predictions from ResFinder JSON."""
    phenotypes = []
    for drug, data in phenotypes_dict.items():
        phenotypes.append(AMRPhenotype(
            drug=drug,
            amr_classes=data["amr_classes"],
            seq_regions=data["seq_regions"],
            seq_variations=data["seq_variations"],
            ref_database=data["ref_database"],
            amr_resistant=data["amr_resistant"],
            amr_species_relevant=data["amr_species_relevant"]
        ))
    return phenotypes

# ---------------------
# The following functions provide a workaround for an issue with ResFinder.
# There appears to be a bug (or possibly a missing implementation) in ResFinder
# that causes resfinder.json to be incomplete when a seq_region is associated
# with both ResFinder and Disinfinder phenotypes.
# In such cases, the phenotypes field includes only the ResFinder drugs,
# while the disinfectant phenotypes lack entries in the seq_region field.
# This missing information can instead be extracted from disinfinder_results_table.txt.


def parse_disinfinder_table(disinfinder_table_file):
    """Parse DisinFinder results table to extract missing phenotype information."""
    if not disinfinder_table_file or not os.path.exists(disinfinder_table_file):
        return {}
    with open(disinfinder_table_file, 'r') as f:
        lines = f.readlines()
    # Find the header line (starts with "Resistance gene")
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('Resistance gene'):
            header_idx = i
            break
    if header_idx is None:
        return {}
    df = pd.read_csv(
        disinfinder_table_file,
        sep='\t',
        skiprows=header_idx
    )
    if df.empty:
        return {}
    # Create mapping of gene name to disinfectant phenotypes
    gene_phenotypes = {}
    for _, row in df.iterrows():
        gene_name = str(row['Resistance gene']).strip()
        phenotype_str = str(row['Phenotype']).strip()
        # Parse comma-separated phenotypes and convert to lowercase
        phenotypes = [
            p.strip().lower() for p in phenotype_str.split(',')
        ]
        gene_phenotypes[gene_name] = phenotypes
    return gene_phenotypes


def update_regions_with_disinfinder(seq_regions, gene_phenotypes, databases):
    """Add DisinFinder phenotypes to seq_regions."""
    for region in seq_regions:
        # Check if this region is associated with DisinFinder
        if not region.ref_database:
            continue
        has_disinfinder = any(
            db.database_name.lower() == 'disinfinder'
            for db in databases
            if db.key in region.ref_database
        )
        if not has_disinfinder:
            continue
        gene_name = region.gene_name or (
            region.key.split(";;")[0] if region.key else None
        )
        if gene_name and gene_name in gene_phenotypes:
            # Add disinfectant phenotypes if they're not already present
            disinfectant_phenotypes = gene_phenotypes[gene_name]
            if region.phenotypes is None:
                region.phenotypes = []
            for phenotype in disinfectant_phenotypes:
                if phenotype not in region.phenotypes:
                    region.phenotypes.append(phenotype)


def update_phenotypes_with_regions(phenotypes, seq_regions):
    """Update phenotype seq_regions with region keys for matching drugs."""
    for phenotype in phenotypes:
        if not phenotype.drug or phenotype.seq_regions:
            continue
        for region in seq_regions:
            if not region.phenotypes or phenotype.drug not in region.phenotypes:
                continue
            if phenotype.seq_regions is None:
                phenotype.seq_regions = []
            if region.key not in phenotype.seq_regions:
                phenotype.seq_regions.append(region.key)
            # Add database keys from region (DisinFinder phenotypes)
            if not phenotype.ref_database and region.ref_database:
                phenotype.ref_database = []
                for db_key in region.ref_database:
                    if 'DisinFinder' in db_key:
                        phenotype.ref_database.append(db_key)
# -------------------------


def parse_resfinder(json_data, disinfinder_table_file=None):
    """Parse complete ResFinder JSON into structured data."""
    databases = parse_databases(json_data.get("databases", {}))
    seq_regions = parse_seq_regions(json_data.get("seq_regions", {}))
    seq_variations = parse_seq_variations(json_data.get("seq_variations", {}))
    phenotypes = parse_phenotypes(json_data.get("phenotypes", {}))
    # If DisinFinder table provided, use it to supplement missing data
    if disinfinder_table_file:
        gene_phenotypes = parse_disinfinder_table(disinfinder_table_file)
        if gene_phenotypes:
            # Add disinfectant phenotypes to regions
            update_regions_with_disinfinder(seq_regions, gene_phenotypes, databases)
            # Update phenotype objects with region associations AND databases
            update_phenotypes_with_regions(phenotypes, seq_regions)
    return {
        "databases": databases,
        "seq_regions": seq_regions,
        "seq_variations": seq_variations,
        "phenotypes": phenotypes,
        "provided_species": json_data.get("provided_species")
    }
