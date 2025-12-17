"""Pytests for workflow_glue modules."""

import math
from pathlib import Path

import pandas as pd
import pytest
from workflow_glue.collect_results import (
    extract_species_from_lineage,
    parse_sourmash_taxonomy,
)
from workflow_glue.models.custom import (
    AMRDatabase, AMRPhenotype, AMRSeqRegion
)
from workflow_glue.pointfinder_species import (
    load_sourmash_pointfinder_mapping,
    sourmash_pointfinder_species,
)
from workflow_glue.process_resfinder_iso import (
    update_phenotypes_with_regions,
    update_regions_with_disinfinder,
)
from workflow_glue.sourmash_picklist import (
    is_valid_assembly_id,
    parse_exclude_file
)


@pytest.fixture
def test_data_location(request):
    """Define test_data location fixture."""
    return request.config.getoption("--test_data_location")


# extract_species_from_lineage
@pytest.mark.parametrize(
    "lineage,expected",
    [
        ("k__Bacte;p__Firm;g__Bacius;s__Bacillus subtilis", "Bacillus subtilis"),
        ("k__Bacte;p__Prot;g__Escher;s__Escherichia coli", "Escherichia coli"),
        ("", "Unknown"),
        (None, "Unknown"),
        (float("nan"), "Unknown"),
    ],
)
def test_extract_species_from_lineage_expected_cases(lineage, expected):
    """Test extract_species_from_lineage with expected cases."""
    if isinstance(lineage, float) and math.isnan(lineage):
        lineage = pd.NA
    assert extract_species_from_lineage(lineage) == expected


def test_extract_species_from_lineage_without_species_prefix_returns_genus():
    """Returns genus when there is no `s__` in the lineage but `g__` exists."""
    lineage = "k__Bacteria;p__Proteobacteria;g__Escherichia"
    assert extract_species_from_lineage(lineage) == "Escherichia"


# parse_sourmash_taxonomy
def test_parse_sourmash_taxonomy(tmp_path: Path):   # noqa: NT001
    """Test parse_sourmash_taxonomy function with various data scenarios."""
    df = pd.DataFrame([
        {
            "gather_result_rank": 0,
            "name": "ref1",
            "lineage": "k__Bacteria;g__Escherichia;s__Escherichia coli",
            "intersect_bp": 100,
            "f_orig_query": 0.9,
            "f_unique_to_query": 0.8,
            "query_containment_ani": 97.5,
            "remaining_bp": 50,
        },
        {
            # NaN fields should convert to None
            "gather_result_rank": 1,
            "name": "ref1",
            "lineage": pd.NA,
            "intersect_bp": pd.NA,
            "f_orig_query": pd.NA,
            "f_unique_to_query": pd.NA,
            "query_containment_ani": pd.NA,
            "remaining_bp": pd.NA,
        },
        {
            # No species prefix
            "gather_result_rank": 2,
            "name": "ref1",
            "lineage": "k__Bacteria;g__Escherichia",
            "intersect_bp": 100,
            "f_orig_query": 0.9,
            "f_unique_to_query": 0.8,
            "query_containment_ani": 97.5,
            "remaining_bp": 50,
        },
        {
            # Species not last in lineage
            "gather_result_rank": 3,
            "name": "ref1",
            "lineage": "k__Bac;g__Esc;s__Escherichia coli;ss__Escherichia coli c",
            "intersect_bp": 100,
            "f_orig_query": 0.9,
            "f_unique_to_query": 0.8,
            "query_containment_ani": 97.5,
            "remaining_bp": 50,
        },
        {
            # No species or genus prefix
            "gather_result_rank": 4,
            "name": "ref1",
            "lineage": "k__Bacteria;o__Escherichia",
            "intersect_bp": 100,
            "f_orig_query": 0.9,
            "f_unique_to_query": 0.8,
            "query_containment_ani": 97.5,
            "remaining_bp": 50,
        },
    ])
    csv_path = tmp_path / "taxonomy.csv"
    df.to_csv(csv_path, index=False)

    species_id = parse_sourmash_taxonomy(csv_path)
    matches = species_id.detected_matches

    assert len(matches) == 5

    # Row 1 — normal
    m1 = matches[0]
    assert m1.rank == 0
    assert m1.species_call == "Escherichia coli"

    # Row 2 — NaNs -> None
    m2 = matches[1]
    assert m2.rank == 1
    assert m2.species_call == "Unknown"
    assert m2.intersect_bp is None

    # Row 3 — no 's__' but has 'g__'
    m3 = matches[2]
    assert m3.rank == 2
    assert m3.species_call == "Escherichia"

    # Row 4 — species not last
    m4 = matches[3]
    assert m4.rank == 3
    assert m4.species_call == "Escherichia coli"

    # Row 5 — no 's__' or 'g__' -> "Unknown"
    m5 = matches[4]
    assert m5.rank == 4
    assert m5.species_call == "Unknown"


@pytest.mark.parametrize(
    "lineage,expected",
    [
        # species not last and 'ss__'
        (
            "k__Bacteria;g__Escherichia;s__Escherichia coli;"
            "ss__Escherichia coli coli",
            "Escherichia coli"
        ),
        # multiple s__ markers -> function currently takes the LAST 's__' segment
        (
            "k__Bacteria;g__Escherichia;s__Escherichia coli;"
            "s__Escherichia notcoli",
            "Escherichia notcoli"
        ),
        # no species prefix but has genus
        ("k__Bacteria;g__Escherichia", "Escherichia"),
        # no species or genus prefix
        ("k__Bacteria;o__SomeOrder", "Unknown"),
        # NaN / empty -> "Unknown"
        (pd.NA, "Unknown"),
        ("", "Unknown"),
    ],
)
def test_extract_species_from_lineage_edge_cases(lineage, expected):
    """Test extract_species_from_lineage with edge cases."""
    assert extract_species_from_lineage(lineage) == expected


def test_parse_sourmash_taxonomy_empty_file(tmp_path: Path):  # noqa: NT001
    """Test parse_sourmash_taxonomy with empty file."""
    empty = pd.DataFrame(columns=[
        "gather_result_rank", "name", "lineage", "intersect_bp",
        "f_orig_query", "f_unique_to_query", "query_containment_ani", "remaining_bp"
    ])
    csv_path = tmp_path / "empty.csv"
    empty.to_csv(csv_path, index=False)

    species_id = parse_sourmash_taxonomy(csv_path)
    assert species_id.detected_matches == []


@pytest.mark.parametrize(
    "assembly_id,expected",
    [
        ("GCA_123456789.1", True),
        ("GCF_987654321.2", True),
        ("GCA_000005825.15", True),
        ("GCA_123456789", False),  # Missing version
        ("GCA_12345678.1", False),  # Too few digits
        ("GCA_1234567890.1", False),  # Too many digits
        ("GCA_abcdefghi.1", False),  # Non-numeric
        ("GCX_123456789.1", False),  # Wrong prefix
        ("", False),  # Empty string
    ],
)
def test_is_valid_assembly_id(assembly_id, expected):
    """Test is_valid_assembly_id with various assembly ID formats."""
    assert is_valid_assembly_id(assembly_id) == expected


def test_is_valid_assembly_id_edge_cases():
    """Test is_valid_assembly_id with edge cases."""
    # Multiple dots
    assert is_valid_assembly_id("GCA_123456789.1.2") is False

    # Missing version number after dot
    assert is_valid_assembly_id("GCA_123456789.") is False

    # Version with leading zeros
    assert is_valid_assembly_id("GCA_123456789.01") is True


@pytest.mark.parametrize(
    "exclude_content,expected_count",
    [
        # Two valid entries
        ("GCA_123456789.1\nGCF_987654321.2", 2),
        # Multiple on first line
        ("GCA_123456789.1 GCF_987654321.2\nGCF_111111111.1", 1),
        # Empty line in middle
        ("GCA_123456789.1\n\nGCF_987654321.2", 2),
        # Empty file
        ("", 0),
        # Only empty lines
        ("\n\n\n", 0),
        # Whitespace line in middle
        ("GCA_123456789.1\n   \nGCF_987654321.2", 2),
    ],
)
def test_parse_exclude_file(tmp_path, exclude_content, expected_count):
    """Test parse_exclude_file with various file contents."""
    exclude_file = tmp_path / "exclude.txt"
    exclude_file.write_text(exclude_content)

    result = parse_exclude_file(str(exclude_file))
    assert len(result) == expected_count


@pytest.mark.parametrize(
    "csv_content,expected",
    [
        # Valid mapping
        ("Sourmash species,PointFinder species\nKlebsiella pneumoniae,klebsiella",
         {"Klebsiella pneumoniae": "klebsiella"}),
        # Empty file
        ("Sourmash species,PointFinder species", {}),
        # Wrong columns
        ("Species,Database\nKlebsiella pneumoniae,klebsiella", {}),
        # Whitespace trimming
        ("Sourmash species,PointFinder species\n  Klebsiella pneumoniae , klebsiella  ",
         {"Klebsiella pneumoniae": "klebsiella"}),
    ],
)
def test_load_sourmash_pointfinder_mapping(
    tmp_path: Path, csv_content, expected  # noqa: NT001
):
    """Test load_sourmash_pointfinder_mapping with various CSV formats."""
    csv_file = tmp_path / "mapping.csv"
    csv_file.write_text(csv_content)
    result = load_sourmash_pointfinder_mapping(str(csv_file), logger=None)
    assert result == expected


@pytest.mark.parametrize(
    "path_input,expected",
    [
        (None, {}),
        ("/nonexistent/file.csv", {}),
    ],
)
def test_load_sourmash_pointfinder_mapping_invalid_paths(path_input, expected):
    """Test load_sourmash_pointfinder_mapping with invalid paths."""
    result = load_sourmash_pointfinder_mapping(path_input, logger=None)
    assert result == expected


@pytest.mark.parametrize(
    "lineage,mapping_content,expected",
    [
        # Successful mapping
        ("k__Bacteria;s__Klebsiella pneumoniae",
         "Sourmash species,PointFinder species\nKlebsiella pneumoniae,klebsiella",
         "klebsiella"),
        # Species not in mapping
        ("k__Bacteria;s__Pseudomonas aeruginosa",
         "Sourmash species,PointFinder species\nKlebsiella pneumoniae,klebsiella",
         "other"),
        # No mapping file
        ("k__Bacteria;s__Klebsiella pneumoniae", None, "other"),
        # Genus only (returns "Unknown" from extract_species_from_lineage)
        ("k__Bacteria;g__Klebsiella",
         "Sourmash species,PointFinder species\nKlebsiella pneumoniae,klebsiella",
         "other"),
    ],
)
def test_sourmash_pointfinder_species(
    tmp_path: Path, lineage, mapping_content, expected  # noqa: NT001
):
    """Test sourmash_pointfinder_species with various scenarios."""
    # Create sourmash taxonomy CSV
    sourmash_df = pd.DataFrame([{"lineage": lineage}])
    sourmash_file = tmp_path / "taxonomy.csv"
    sourmash_df.to_csv(sourmash_file, index=False)

    # Handle mapping file
    mapping_file = None
    if mapping_content:
        mapping_file = tmp_path / "mapping.csv"
        mapping_file.write_text(mapping_content)
        mapping_file = str(mapping_file)

    result = sourmash_pointfinder_species(str(sourmash_file), mapping_file, logger=None)
    assert result == expected


def test_sourmash_pointfinder_species_empty_file(tmp_path: Path):  # noqa: NT001
    """Test sourmash_pointfinder_species with empty sourmash file."""
    sourmash_df = pd.DataFrame(columns=["lineage"])
    sourmash_file = tmp_path / "empty.csv"
    sourmash_df.to_csv(sourmash_file, index=False)

    result = sourmash_pointfinder_species(str(sourmash_file), None, logger=None)
    assert result == "other"


@pytest.mark.parametrize(
    "lineage,mapping_content,expected",
    [
        # Underscore suffix
        (
            "k__Bacteria;s__Helicobacter pylori_A",
            "Sourmash species,PointFinder species\n"
            "Helicobacter pylori,helicobacter_pylori",
            "helicobacter_pylori"
        ),
        # Genus-level mapping
        (
            "k__Bacteria;s__Helicobacter pylori_A",
            "Sourmash species,PointFinder species\n"
            "Helicobacter,helicobacter",
            "helicobacter"
        ),
        # No delimiter
        (
            "k__Bacteria;s__Helicobacter pyloriA",
            "Sourmash species,PointFinder species\n"
            "Helicobacter pylori,helicobacter_pylori",
            "other"
        ),
        # Genus-only species call without genus-level mapping
        (
            "k__Bacteria;g__Helicobacter",
            "Sourmash species,PointFinder species\n"
            "Helicobacter pylori,helicobacter",
            "other"
        ),
    ],
)
def test_sourmash_pointfinder_species_suffix_matching(
    tmp_path: Path, lineage, mapping_content, expected  # noqa: NT001
):
    """Test species suffix matching."""
    sourmash_df = pd.DataFrame([{"lineage": lineage}])
    sourmash_file = tmp_path / "taxonomy.csv"
    sourmash_df.to_csv(sourmash_file, index=False)

    mapping_file = tmp_path / "mapping.csv"
    mapping_file.write_text(mapping_content)

    result = sourmash_pointfinder_species(
        str(sourmash_file), str(mapping_file), logger=None
    )
    assert result == expected


def test_sourmash_pointfinder_species_regex_escape(
    tmp_path: Path  # noqa: NT001
):
    """Test special regex characters in species names are escaped."""
    lineage = "k__Bacteria;s__Helicobacter (species) group_strain1"
    mapping_content = (
        "Sourmash species,PointFinder species\n"
        "Helicobacter (species) group,helicobacter"
    )
    sourmash_df = pd.DataFrame([{"lineage": lineage}])
    sourmash_file = tmp_path / "taxonomy.csv"
    sourmash_df.to_csv(sourmash_file, index=False)

    mapping_file = tmp_path / "mapping.csv"
    mapping_file.write_text(mapping_content)

    result = sourmash_pointfinder_species(
        str(sourmash_file), str(mapping_file), logger=None
    )
    assert result == "helicobacter"


# ------------ update_regions_with_disinfinder ----------------
def test_update_regions_with_disinfinder_adds_phenotypes():
    """Test that DisinFinder phenotypes are added to matching regions."""
    databases = [AMRDatabase(
        key="DisinFinder-1.0",
        database_name="DisinFinder",
        database_version="1.0")]
    region = AMRSeqRegion(
        key="geneA;;1;;ABC123",
        gene_name="geneA",
        ref_database=["DisinFinder-1.0"],
        phenotypes=None)
    gene_phenotypes = {"geneA": ["compound_one", "compound_two"]}
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    assert region.phenotypes == ["compound_one", "compound_two"]


def test_update_regions_with_disinfinder_no_duplicates():
    """Test that duplicate phenotypes are not added."""
    databases = [AMRDatabase(
        key="DisinFinder-1.0",
        database_name="DisinFinder",
        database_version="1.0")]
    region = AMRSeqRegion(
        key="geneA;;1;;ABC123",
        gene_name="geneA",
        ref_database=["DisinFinder-1.0"],
        phenotypes=["compound_one"])
    gene_phenotypes = {"geneA": ["compound_one", "compound_two"]}
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    assert region.phenotypes.count("compound_one") == 1


def test_update_regions_with_disinfinder_partial_overlap():
    """Test that new phenotypes are added when region already has some but not all."""
    databases = [AMRDatabase(
        key="DisinFinder-1.0",
        database_name="DisinFinder",
        database_version="1.0")]
    region = AMRSeqRegion(
        key="geneA;;1;;ABC123",
        gene_name="geneA",
        ref_database=["DisinFinder-1.0"],
        phenotypes=["compound_one", "compound_three"]
    )
    gene_phenotypes = {"geneA": ["compound_one", "compound_two", "compound_four"]}
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    assert "compound_one" in region.phenotypes
    assert "compound_two" in region.phenotypes
    assert "compound_three" in region.phenotypes
    assert "compound_four" in region.phenotypes
    assert region.phenotypes.count("compound_one") == 1


def test_update_regions_with_disinfinder_gene_not_in_mapping():
    """Test that regions are unchanged when gene not in phenotypes mapping."""
    databases = [AMRDatabase(
        key="DisinFinder-1.0",
        database_name="DisinFinder",
        database_version="1.0")]
    region = AMRSeqRegion(
        key="geneA;;1;;ABC123",
        gene_name="geneA",
        ref_database=["DisinFinder-1.0"],
        phenotypes=["existing_compound"]
    )
    gene_phenotypes = {"geneB": ["compound_one"]}  # Different gene
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    # Should remain unchanged
    assert region.phenotypes == ["existing_compound"]


def test_update_regions_with_disinfinder_empty_phenotypes_list():
    """Test that empty phenotypes list in mapping doesn't break function."""
    databases = [AMRDatabase(
        key="DisinFinder-1.0",
        database_name="DisinFinder",
        database_version="1.0")]
    region = AMRSeqRegion(
        key="geneA;;1;;ABC123",
        gene_name="geneA",
        ref_database=["DisinFinder-1.0"],
        phenotypes=["existing_compound"]
    )
    gene_phenotypes = {"geneA": []}  # Empty list
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    # Should remain unchanged with original phenotypes
    assert region.phenotypes == ["existing_compound"]


def test_update_regions_with_disinfinder_extracted_gene_not_in_mapping():
    """Test when gene extracted from key is not in phenotypes mapping."""
    databases = [AMRDatabase(
        key="DisinFinder-1.0",
        database_name="DisinFinder",
        database_version="1.0")]
    region = AMRSeqRegion(
        key="unknownGene;;1;;ABC123",
        gene_name=None,  # Will extract "unknownGene" from key
        ref_database=["DisinFinder-1.0"],
        phenotypes=None
    )
    gene_phenotypes = {"geneA": ["compound_one"]}  # unknownGene not in mapping
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    # Should remain unchanged since extracted gene not in mapping
    assert region.phenotypes is None


@pytest.mark.parametrize(
    "ref_database,database_name,should_update",
    [
        (["ResFinder-1.0"], "ResFinder", False),  # Non-DisinFinder
        (None, "DisinFinder", False),  # No ref_database
        (["DisinFinder-1.0"], "DisinFinder", True),  # Valid DisinFinder
    ],
)
def test_update_regions_with_disinfinder_database_filtering(
    ref_database, database_name, should_update
):
    """Test that only DisinFinder db associated regions are updated."""
    databases = [AMRDatabase(
        key=f"{database_name}-1.0",
        database_name=database_name,
        database_version="1.0")]
    region = AMRSeqRegion(
        key="geneA;;1;;ABC",
        gene_name="geneA",
        ref_database=ref_database,
        phenotypes=None)
    gene_phenotypes = {"geneA": ["phenotype_one"]}
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    if should_update:
        assert region.phenotypes == ["phenotype_one"]
    else:
        assert region.phenotypes is None


def test_update_regions_with_disinfinder_extracts_gene_from_key():
    """Test that gene name is extracted from key when gene_name is None."""
    databases = [AMRDatabase(
        key="DisinFinder-1.0",
        database_name="DisinFinder",
        database_version="1.0")]
    region = AMRSeqRegion(
        key="geneA;;1;;ABC123",
        gene_name=None,
        ref_database=["DisinFinder-1.0"],
        phenotypes=None)
    gene_phenotypes = {"geneA": ["compound_one"]}
    update_regions_with_disinfinder([region], gene_phenotypes, databases)
    assert region.phenotypes == ["compound_one"]


# ---------- update_phenotypes_with_regions -----------------
def test_update_phenotypes_with_regions_adds_region_keys():
    """Test that region keys are added to matching phenotypes."""
    phenotype = AMRPhenotype(
        drug="drug_one",
        seq_regions=None,
        ref_database=None)
    region = AMRSeqRegion(key="geneB;;1;;XYZ789", phenotypes=["drug_one"])
    update_phenotypes_with_regions([phenotype], [region])
    assert phenotype.seq_regions == ["geneB;;1;;XYZ789"]


def test_update_phenotypes_with_regions_no_duplicates():
    """Test that duplicate region keys are not added."""
    phenotype = AMRPhenotype(
        drug="drug_one",
        seq_regions=["geneB;;1;;XYZ789"],
        ref_database=None)
    region = AMRSeqRegion(key="geneB;;1;;XYZ789", phenotypes=["drug_one"])
    update_phenotypes_with_regions([phenotype], [region])
    assert len(phenotype.seq_regions) == 1


def test_update_phenotypes_with_regions_phenotype_not_in_region():
    """Test that phenotype is unchanged when drug is not in any region's phenotypes."""
    phenotype = AMRPhenotype(
        drug="drug_nonexistent",
        seq_regions=None,
        ref_database=None)
    region = AMRSeqRegion(key="geneA;;1;;ABC", phenotypes=["drug_one", "drug_two"])
    update_phenotypes_with_regions([phenotype], [region])
    # Should remain unchanged
    assert phenotype.seq_regions is None
    assert phenotype.ref_database is None


@pytest.mark.parametrize(
    "phenotype_drug,phenotype_seq_regions,region_phenotypes,should_update",
    [
        (None, None, ["drug_one"], False),  # No drug
        ("drug_one", ["existing"], ["drug_one"], False),  # Has seq_regions
        ("drug_one", None, None, False),  # Region has no phenotypes
        ("drug_one", None, ["other_drug"], False),  # Drug not in region
        ("drug_one", None, ["drug_one"], True),  # Valid match
    ],
)
def test_update_phenotypes_with_regions_filtering(
    phenotype_drug, phenotype_seq_regions, region_phenotypes, should_update
):
    """Test phenotype filtering logic."""
    phenotype = AMRPhenotype(
        drug=phenotype_drug,
        seq_regions=phenotype_seq_regions,
        ref_database=None)
    region = AMRSeqRegion(key="geneC;;1;;ABC", phenotypes=region_phenotypes)
    update_phenotypes_with_regions([phenotype], [region])
    if should_update:
        assert phenotype.seq_regions == ["geneC;;1;;ABC"]
    else:
        assert phenotype.seq_regions == phenotype_seq_regions


def test_update_phenotypes_with_regions_adds_disinfinder_database():
    """Test that DisinFinder database keys are added to phenotypes."""
    phenotype = AMRPhenotype(
        drug="compound_one",
        seq_regions=None,
        ref_database=None)
    region = AMRSeqRegion(
        key="geneA;;1;;ABC",
        phenotypes=["compound_one"],
        ref_database=["DisinFinder-1.0"])
    update_phenotypes_with_regions([phenotype], [region])
    assert phenotype.ref_database == ["DisinFinder-1.0"]


def test_update_phenotypes_with_regions_skips_non_disinfinder():
    """Test that non-DisinFinder database keys are not added."""
    phenotype = AMRPhenotype(
        drug="drug_one",
        seq_regions=None,
        ref_database=None)
    region = AMRSeqRegion(
        key="geneB;;1;;XYZ",
        phenotypes=["drug_one"],
        ref_database=["ResFinder-1.0"])
    update_phenotypes_with_regions([phenotype], [region])
    # ref_database gets initialized but should not contain ResFinder entries
    assert phenotype.ref_database is not None
    assert len(phenotype.ref_database) == 0
