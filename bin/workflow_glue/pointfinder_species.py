#!/usr/bin/env python
"""Obtain pointfinder species option from mlst or Sourmash results."""

import json
import os
import re
import sys

import pandas as pd

from .collect_results import extract_species_from_lineage  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


def load_sourmash_pointfinder_mapping(mapping_csv_path, logger):
    """Load sourmash to pointfinder mapping from CSV file."""
    if not mapping_csv_path or not os.path.exists(mapping_csv_path):
        return {}
    df = pd.read_csv(mapping_csv_path)

    if 'Sourmash species' not in df.columns or 'PointFinder species' not in df.columns:
        if logger:
            logger.warning(
                f"CSV file {mapping_csv_path} must contain 'Sourmash species' "
                f"and 'PointFinder species' columns, proceeding without custom mapping"
            )
        return {}

    if df.empty:
        if logger:
            logger.warning(
                f"CSV file {mapping_csv_path} is empty, proceeding without custom "
                f"mapping"
            )
        return {}

    mapping = dict(
        zip(
            df['Sourmash species'].str.strip(),
            df['PointFinder species'].str.strip(),
        )
    )

    return mapping


def sourmash_pointfinder_species(sourmash_csv_path, mapping_csv_path, logger):
    """Try to get species from sourmash if available."""
    df = pd.read_csv(sourmash_csv_path)
    if df.empty:
        return "other"

    # Extract species from lineage
    lineage = df.iloc[0].get('lineage', '')
    species = extract_species_from_lineage(lineage)

    if not species:
        return "other"

    if not mapping_csv_path:
        return "other"

    sourmash_to_pointfinder = load_sourmash_pointfinder_mapping(
        mapping_csv_path,
        logger,
    )

    # Allow exact matches and matches with suffixes after species names
    for mapping_species, pointfinder_name in sourmash_to_pointfinder.items():
        pattern = re.escape(mapping_species) + r"(?:[_\s.\-]|$)"
        if re.match(pattern, species):
            return pointfinder_name

    return "other"


def main(args):
    """Run entry point."""
    logger = get_named_logger("pointfinder_species")
    """Extract mlst scheme and assign pointfinder species."""
    pointfinder_dict = {
        "campylobacter_nonjejuni_7": "campylobacter",
        "campylobacter_nonjejuni_8": "campylobacter",
        "campylobacter_nonjejuni": "campylobacter",
        "campylobacter_nonjejuni_6": "campylobacter",
        "campylobacter_nonjejuni_3": "campylobacter",
        "campylobacter_nonjejuni_5": "campylobacter",
        "campylobacter": "campylobacter",
        "campylobacter_nonjejuni_4": "campylobacter",
        "campylobacter_nonjejuni_2": "campylobacter",
        "efaecium": "enterococcus faecium",
        "efaecalis": "enterococcus faecalis",
        "neisseria": "neisseria gonorrhoeae",
        "senterica_achtman_2": "salmonella",
        "ecoli": 'escherichia_coli',
        "klebsiella": "klebsiella",
        "koxytoca": "klebsiella",
        "kaerogenes": "klebsiella",
        "saureus": "staphylococcus aureus",
        "helicobacter": "helicobacter pylori",
        "mycobacteria_2": "mycobacterium tuberculosis",
    }
    with open(args.mlst_json) as f:
        data = json.load(f)
    pointfinder_species = pointfinder_dict.get(data[0]["scheme"], "other")

    if (
        pointfinder_species == "other"
        and hasattr(args, 'sourmash_csv')
        and args.sourmash_csv
    ):
        mapping_csv = getattr(args, 'mapping_csv', None)
        pointfinder_species = sourmash_pointfinder_species(
            args.sourmash_csv,
            mapping_csv,
            logger,
        )

    logger.info("Pointfinder species identified.")
    sys.stdout.write(pointfinder_species)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("pointfinder_species")
    parser.add_argument(
        "--mlst_json",
        help="MLST json results file")
    parser.add_argument(
        "--sourmash_csv",
        help="Sourmash taxonomy CSV file (optional)")
    parser.add_argument(
        "--mapping_csv",
        help="CSV file with Sourmash to PointFinder species mappings (optional)")
    return parser
