"""Create a sourmash picklist from combined databases, exclude specified assemblies."""

import os
import re
import zipfile

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def is_valid_assembly_id(assembly_id):
    """Check if assembly ID matches GCA_XXXXXXXXX.V or GCF_XXXXXXXXX.V format."""
    pattern = r'^GC[AF]_\d{9}\.\d+$'
    return re.match(pattern, assembly_id) is not None


def parse_exclude_file(file_path, logger=None):
    """Parse exclude file and return set of valid assembly IDs to exclude."""
    assemblies_to_remove = set()

    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            # Check for multiple entries on one line
            entries = line.split()
            if len(entries) > 1:
                if logger:
                    logger.warning(
                        f"Line {line_num} contains multiple entries, "
                        f"skipping: {line}"
                    )
                continue
            # Check format with regex validation
            if is_valid_assembly_id(line):
                assemblies_to_remove.add(line)
            else:
                if logger:
                    logger.warning(
                        f"Line {line_num} invalid format (expected "
                        f"GCA_XXXXXXXXX.V or GCF_XXXXXXXXX.V), skipping: {line}"
                    )
    return assemblies_to_remove


def main(args):
    """Run entry point."""
    logger = get_named_logger("sourmash_picklist")

    # Read assemblies to exclude
    assemblies_to_remove = set()
    if args.exclude_file and os.path.exists(args.exclude_file):
        assemblies_to_remove = parse_exclude_file(args.exclude_file, logger)
    else:
        if args.exclude_file:
            logger.warning(
                f"Exclude file {args.exclude_file} not found, "
                f"proceeding without exclusions"
            )
        else:
            logger.info("No Sourmash DB exclusion file provided")
    logger.info(
        f"Total assemblies to exclude from Sourmash DB: {len(assemblies_to_remove)}"
    )

    # Extract manifests from both databases
    manifests = []
    for db_path in [args.fungi_db, args.bacteria_db]:
        with zipfile.ZipFile(db_path, 'r') as z:
            with z.open('SOURMASH-MANIFEST.csv') as f:
                manifests.append(pd.read_csv(f, comment='#'))

    # Combine manifests
    combined_df = pd.concat(manifests, ignore_index=True)

    # Filter out assemblies to remove
    excluded_assemblies = set()
    if assemblies_to_remove:
        mask = combined_df['name'].apply(
            lambda x: not any(assembly in str(x) for assembly in assemblies_to_remove)
        )
        # Find which assemblies were actually found and excluded
        excluded_mask = ~mask
        excluded_entries = combined_df[excluded_mask]
        for assembly in assemblies_to_remove:
            if any(assembly in str(name) for name in excluded_entries['name']):
                excluded_assemblies.add(assembly)
        combined_df = combined_df[mask]

    # Create picklist with required columns: md5, name
    picklist = combined_df[['md5', 'name']]
    picklist.to_csv(args.output, index=False)

    excluded_file = "sourmash_picklist_excluded.txt"
    with open(excluded_file, 'w') as f:
        for assembly in sorted(excluded_assemblies):
            f.write(f"{assembly}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("sourmash_picklist")
    parser.add_argument(
        "--fungi-db",
        required=True,
        help="Path to fungal database zip file"
    )
    parser.add_argument(
        "--bacteria-db",
        required=True,
        help="Path to bacterial database zip file"
    )
    parser.add_argument(
        "--exclude-file",
        help="File containing assembly IDs to exclude (one per line)"
    )
    parser.add_argument(
        "--output",
        default="sourmash_picklist.csv",
        help="Output picklist CSV file (default: sourmash_picklist.csv)"
    )
    return parser
