"""Extended models for the workflow."""
from dataclasses import dataclass, field
from enum import Enum
import os

import pandas as pd
from workflow_glue.models.common import (
    ResultsContents as BaseResultsContents,
    Sample as BaseSample,
    WorkflowResult as BaseWorkflowResult
)
from workflow_glue.report_utils import convert_bp


@dataclass
class FastqStats:
    """Read statistics."""

    n_seqs: int | None = field(
        default=None,
        metadata={
            "title": "Number of reads",
            "description": "The number of sequencing reads"})
    n_bases: int | None = field(
        default=None,
        metadata={
            "title": "Number of bases",
            "description": "The number of bases"})
    min_length: int | None = field(
        default=None,
        metadata={
            "title": "Minimum read length",
            "description": "The minimum read length"})
    max_length: int | None = field(
        default=None,
        metadata={
            "title": "Maximum read length",
            "description": "The maximum read length"})
    mean_quality: float | None = field(
        default=None,
        metadata={
            "title": "Mean read quality",
            "description": "The mean read quality"})
    median_length: float | None = field(
        default=None,
        metadata={
            "title": "Median read length",
            "description": "The median read length"})
    median_quality: float | None = field(
        default=None,
        metadata={
            "title": "Median read quality",
            "description": "The median read quality"})


@dataclass
class AMRSeqRegion:
    """Individual sequence region from ResFinder analysis."""

    key: str | None = field(
        default=None,
        metadata={
            "title": "Region key",
            "description": "Full ResFinder key identifier for this sequence region"})
    phenotypes: list[str] | None = field(
        default=None,
        metadata={
            "title": "Associated phenotypes",
            "description": "List of drug resistance phenotypes associated with "
            "this region"})
    ref_database: list[str] | None = field(
        default=None,
        metadata={
            "title": "Reference databases",
            "description": "List of identifying this sequence region"})
    gene_name: str | None = field(
        default=None,
        metadata={
            "title": "Name",
            "description": "Gene name"})
    ref_acc: str | None = field(
        default=None,
        metadata={
            "title": "Reference accession",
            "description": "Reference sequence accession number"})
    identity: float | None = field(
        default=None,
        metadata={
            "title": "Query to reference identity",
            "description": "Query sequence identity percentage to reference"})
    coverage: float | None = field(
        default=None,
        metadata={
            "title": "Query to reference coverage",
            "description": "Query coverage percentage of reference sequence"})
    alignment_length: int | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Alignment length",
            "description": "Length of sequence alignment in base pairs"})
    ref_seq_length: int | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Reference sequence length",
            "description": "Length of reference sequence in base pairs"})
    ref_start_pos: int | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Reference start position",
            "description": "Start position on reference sequence"})
    ref_end_pos: int | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Reference end position",
            "description": "End position on reference sequence"})
    query_contig: str | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Assembly contig on which the variant was detected"})
    query_start_pos: int | None = field(
        default=None,
        metadata={
            "title": "Query start position",
            "description": "Query start position within contig"})
    query_end_pos: int | None = field(
        default=None,
        metadata={
            "title": "Query end position",
            "description": "Query end position within contig"})
    pmids: list[str] | None = field(
        default=None,
        metadata={
            "title": "PubMed IDs",
            "description": "List of reference citations"})
    notes: list[str] | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Notes",
            "description": "Database notes"})
    grade: int | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Confidence grade",
            "description": "Confidence grade for this resistance prediction (1-3)"})

    def get_gene_name(self):
        """Extract clean gene name."""
        if self.gene_name:
            return self.gene_name
        if self.key:
            return self.key.split(";;")[0]
        return "Unknown"


class MutationType(str, Enum):
    """Types of genetic mutations."""

    substitution = "substitution"
    deletion = "deletion"
    insertion = "insertion"


@dataclass
class AMRSeqVariation:
    """Point mutations from ResFinder analysis."""

    key: str | None = field(
        default=None,
        metadata={
            "title": "Variation key",
            "description": "Full ResFinder key identifier for this sequence variation"})
    ref_database: str | None = field(
        default=None,
        metadata={
            "title": "Reference database",
            "description": "Database used to identify this variation"})
    seq_regions: list[str] | None = field(
        default=None,
        metadata={
            "title": "Associated AMRSeqRegion",
            "description": "List of sequence region keys associated with "
            "this variation"})
    phenotypes: list[str] | None = field(
        default=None,
        metadata={
            "title": "Associated phenotypes",
            "description": "List of drug resistance phenotypes caused by "
            "this variation"})
    seq_var: str | None = field(
        default=None,
        metadata={
            "title": "Sequence variation",
            "description": "Sequence variation in standard notation (e.g. p.P161R)"})
    ref_codon: str | None = field(
        default=None,
        metadata={
            "title": "Reference codon",
            "description": "Reference codon sequence"})
    var_codon: str | None = field(
        default=None,
        metadata={
            "title": "Variant codon",
            "description": "Variant codon sequence"})
    ref_aa: str | None = field(
        default=None,
        metadata={
            "title": "Reference amino acid",
            "description": "Reference amino acid single letter code"})
    var_aa: str | None = field(
        default=None,
        metadata={
            "title": "Variant amino acid",
            "description": "Variant amino acid single letter code"})
    codon_change: str | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Codon change",
            "description": "Codon change notation (e.g. ccg>cgg)"})
    nuc_change: str | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Nucleotide change",
            "description": "Single nucleotide or amino acid change"})
    ref_start_pos: int | None = field(
        default=None,
        metadata={
            "title": "Reference start position",
            "description": "Start position of variation on reference sequence, "
            "relative to the gene start position"})
    ref_end_pos: int | None = field(
        default=None,
        metadata={
            "title": "Reference end position",
            "description": "End position of variation on reference sequence, "
            "relative to the gene start position"})
    mut_type: MutationType | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Mutation type",
            "description": "Type of mutation: substitution, deletion, or insertion"})
    pmids: list[str] | None = field(
        default=None,
        metadata={
            "title": "PubMed IDs",
            "description": "List of PubMed identifiers supporting this variation"})
    notes: list[str] | None = field(
        default=None,
        metadata={
            "title": "Notes",
            "description": "Database notes"})

    def get_gene_name(self):
        """Extract clean gene name from key."""
        if not self.key:
            return "Unknown"
        # Extract: "embA-promoter n.-11C>A;;1;;CP068259.1" -> "embA"
        gene_key = self.key.split(";;")[0]
        gene_name = gene_key.split()[0]  # Remove mutation notation
        gene_name = gene_name.replace("-promoter", "")  # Remove promoter suffix
        return gene_name

    def is_promoter_mutation(self):
        """Check if this is a promoter mutation."""
        if self.key and "-promoter" in self.key:
            return True
        if self.ref_start_pos is not None and self.ref_start_pos < 0:
            return True
        return False

    def has_fungamr_notes(self):
        """Check if this variation has FungAMR database entries."""
        if not self.pmids:
            return False
        return any("fungamr-db" in str(pmid).lower() for pmid in self.pmids)

    def format_pmids_display(self):
        """Format PMIDs for display, replacing fungamr-db entries."""
        if not self.pmids:
            return None
        formatted = []
        for pmid in self.pmids:
            pmid_str = str(pmid)
            if "fungamr-db" in pmid_str.lower():
                formatted.append("FungAMR-DB")
            else:
                formatted.append(pmid_str)
        return ", ".join(formatted)

    def get_contig_and_positions(
        self, all_regions: list['AMRSeqRegion']
    ) -> tuple[str | None, int | None, int | None]:
        """Get contig name and convert gene-relative positions to contig positions."""
        if not self.seq_regions or not all_regions:
            return None, self.ref_start_pos, self.ref_end_pos
        # Find the gene region this mutation belongs to
        for region in all_regions:
            if region.key in self.seq_regions:
                contig = region.query_contig
                # Calculate contig positions if possible
                if (
                    self.ref_start_pos is not None
                    and region.query_start_pos is not None
                    and region.query_end_pos is not None
                ):
                    if region.query_end_pos < region.query_start_pos:
                        # Negative strand
                        start = region.query_start_pos - self.ref_start_pos
                        end = (
                            region.query_start_pos - self.ref_end_pos
                            if self.ref_end_pos is not None else start
                        )
                    else:
                        # Positive strand
                        start = region.query_start_pos + self.ref_start_pos
                        end = (
                            region.query_start_pos + self.ref_end_pos
                            if self.ref_end_pos is not None else start
                        )
                else:
                    start = self.ref_start_pos
                    end = self.ref_end_pos
                return contig, start, end
        # Couldn't find matching region - return None and original positions
        return None, self.ref_start_pos, self.ref_end_pos


@dataclass
class AMRPhenotype:
    """Drug resistance phenotype prediction."""

    drug: str | None = field(
        default=None,
        metadata={
            "title": "Drug name",
            "description": "Name of antimicrobial drug"})
    amr_classes: list[str] | None = field(
        default=None,
        metadata={
            "title": "Antimicrobial classes",
            "description": "List of antimicrobial drug classes"})
    seq_regions: list[str] | None = field(
        default=None,
        metadata={
            "title": "Supporting sequence regions",
            "description": "List of sequence regions supporting this phenotype"})
    seq_variations: list[str] | None = field(
        default=None,
        metadata={
            "title": "Supporting sequence variations",
            "description": "List of sequence variation hits with this phenotype"})
    ref_database: list[str] | None = field(
        default=None,
        metadata={
            "title": "Reference databases",
            "description": "List of databases containing this phenotype prediction"})
    amr_resistant: bool | None = field(
        default=None,
        metadata={
            "title": "Resistance confirmed",
            "description": "Boolean indicating confirmed antimicrobial resistance"})
    amr_species_relevant: bool | None = field(
        default=None,
        metadata={
            "title": "Species relevant",
            "description": "Boolean indicating if this resistance is relevant "
            "for the analyzed species"})


@dataclass
class AMRDatabase:
    """Database version information."""

    key: str | None = field(
        default=None,
        metadata={
            "title": "Database key",
            "description": "Unique database identifier with version"})
    database_name: str | None = field(
        default=None,
        metadata={
            "title": "Database name",
            "description": "Name of the resistance database"})
    database_version: str | None = field(  # Not used in reports
        default=None,
        metadata={
            "title": "Database version",
            "description": "Version number of the database"})


@dataclass
class SequenceTypeSchema:
    """MLST schema and allele variant identified for sample."""

    schema_identifier: str | None = field(
        default=None,
        metadata={
            "title": "Schema identifier",
            "description": "Identifier for the schema used to classify the sequence"})
    allele_variant: str | None = field(
        default=None,
        metadata={
            "title": "Allele variant",
            "description": "Identifier for the specific allele variant detected"})


@dataclass
class Serotype:
    """Salmonella serotyping results."""

    predicted_serotype: str | None = field(
        default=None,
        metadata={
            "title": "Predicted serotype",
            "description": """Predicited serological typing idenitifer for the sample
            (Salmonella only)"""})
    predicted_antigenic_profile: str | None = field(
        default=None,
        metadata={
            "title": "Predicted antigenic profile",
            "description": """Predicted antigenic profile for the sample
            using the O, H1 and H2 antigens identified (Salmonella only)"""})
    predicted_identification: str | None = field(
        default=None,
        metadata={
            "title": "Predicted identification",
            "description": """Predicted identification for the sample
            (Salmonella only)"""})
    o_antigen_prediction: str | None = field(
        default=None,
        metadata={
            "title": "O antigen prediction",
            "description": """Predicted O antigen found in the sample
            (Salmonella only)"""})
    h1_antigen_prediction: str | None = field(
        default=None,
        metadata={
            "title": "H1 antigen prediction",
            "description": """Predicted H1 antigen found in the sample
            (Salmonella only)"""})
    h2_antigen_prediction: str | None = field(
        default=None,
        metadata={
            "title": "H2 antigen prediction",
            "description": """Predicted H2 antigen found in the sample
            (Salmonella only)"""})
    note: str | None = field(
        default=None,
        metadata={
            "title": "Note",
            "description": "Additional notes from serotyping analysis"})


@dataclass
class Annotation:
    """Region of interest identified within assembly."""

    contig: str | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Assembly contig on which region was detected"})
    ID: str | None = field(
        default=None,
        metadata={
            "title": "Identifier",
            "description": "Unique identifier for the annotation from Bakta"})
    start: int | None = field(
        default=None,
        metadata={
            "title": "Start position",
            "description": "Start position of the region on the assembly contig"})
    end: int | None = field(
        default=None,
        metadata={
            "title": "End position",
            "description": "End position of the region on the assembly contig"})
    strand: str | None = field(
        default=None,
        metadata={
            "title": "Strand",
            "description": "Which strand the region was identified on"})
    gene: str | None = field(
        default=None,
        metadata={
            "title": "Gene",
            "description": "Gene name of the region of interest, if applicable"})
    product: str | None = field(
        default=None,
        metadata={
            "title": "Product",
            "description": "Product name of the gene or region, if applicable"})
    ec_number: str | None = field(
        default=None,
        metadata={
            "title": "EC number",
            "description": "Identifier from the enzyme consortium catalogue"})


@dataclass
class Variant:
    """Variants identified in assembly compared to reference."""

    contig: str | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Contig on which the variant was detected"})
    pos: int | None = field(
        default=None,
        metadata={
            "title": "Position",
            "description": ""})
    ref: str | None = field(
        default=None,
        metadata={
            "title": "Reference",
            "description": "The reference nucleotide at the variant site"})
    alt: str | None = field(
        default=None,
        metadata={
            "title": "Alternate allele",
            "description": "The alternate nucleotide at the variant site"})
    depth: float | None = field(
        default=None,
        metadata={
            "title": "Depth",
            "description": "The frequency of the variant in the sample"})


@dataclass
class Coverage:
    """Coverage summary information for each contig in assembly."""

    counts: int | None = field(
        default=None,
        metadata={
            "title": "Counts",
            "description": "Number of reads that map to the contig"})
    median: float | None = field(
        default=None,
        metadata={
            "title": "Median coverage",
            "description": "Median coverage"})
    mean: float | None = field(
        default=None,
        metadata={
            "title": "Mean coverage",
            "description": "Mean coverage"})
    minimum: int | None = field(
        default=None,
        metadata={
            "title": "Minimum coverage",
            "description": "Minimum coverage"})
    maximum: int | None = field(
        default=None,
        metadata={
            "title": "Maximum coverage",
            "description": "Maximum coverage"})


@dataclass
class FlyeContigStats:
    """Flye assembly statistics for a single contig."""

    seq_name: str | None = field(
        default=None,
        metadata={
            "title": "Sequence name",
            "description": "Name of the assembled contig"})
    length: int | None = field(
        default=None,
        metadata={
            "title": "Length",
            "description": "Length of the contig in base pairs"})
    coverage: float | None = field(
        default=None,
        metadata={
            "title": "Coverage",
            "description": "Estimated coverage of the contig"})
    circular: str | None = field(
        default=None,
        metadata={
            "title": "Circular",
            "description": "Whether the contig is circular (Y/N)"})
    repeat: str | None = field(
        default=None,
        metadata={
            "title": "Repeat",
            "description": "Whether the contig is a repeat (Y/N)"})
    multiplicity: int | None = field(
        default=None,
        metadata={
            "title": "Multiplicity",
            "description": "Copy number of the contig"})


@dataclass
class AntimicrobialResistance:
    """Complete ResFinder analysis results for a sample."""

    seq_regions: list[AMRSeqRegion] | None = field(
        default=None,
        metadata={
            "title": "Sequence regions",
            "description": "List of resistance gene sequence regions identified"})
    seq_variations: list[AMRSeqVariation] | None = field(
        default=None,
        metadata={
            "title": "Sequence variations",
            "description": "List of point mutations and sequence variations "
            "identified"})
    phenotypes: list[AMRPhenotype] | None = field(
        default=None,
        metadata={
            "title": "Phenotype predictions",
            "description": "List of predicted antimicrobial resistance phenotypes"})
    databases: list[AMRDatabase] | None = field(
        default=None,
        metadata={
            "title": "Database information",
            "description": "List of databases used in the analysis"})
    provided_species: str | None = field(
        default=None,
        metadata={
            "title": "Species analyzed",
            "description": "Species name provided for analysis"})
    pointfinder_available_genes: list[str] | None = field(
        default=None,
        metadata={
            "title": "Available PointFinder genes",
            "description": "List of all genes available for point mutation analysis"
            " in PointFinder database"
        })

    def get_acquired_resistance_display(self, fill_na="-"):  # used in per_sample_report
        """Convert acquired resistance genes to display format for per-sample report."""
        if not self.seq_regions:
            return {}
        db_lookup = {}
        if self.databases:
            for db in self.databases:
                db_lookup[db.key] = db.database_name
        acquired_dict = {}
        for region in self.seq_regions:
            # Skip regions without phenotypes, these are for point mutations
            if not region.phenotypes:
                continue
            gene = region.get_gene_name()
            db_names = []
            if region.ref_database:
                for db_key in region.ref_database:
                    if db_key in db_lookup:
                        db_names.append(db_lookup[db_key])
            acquired_dict[gene] = {
                "drugs": region.phenotypes,
                "contig": region.query_contig,
                "start": region.query_start_pos,
                "end": region.query_end_pos,
                "database": ", ".join(db_names) if db_names else fill_na,
                "identity": (
                    f"{region.identity: .2f}%" if region.identity else fill_na),
                "coverage": (
                    f"{region.coverage: .2f}%" if region.coverage else fill_na),
                "pmids": ", ".join(region.pmids) if region.pmids else fill_na
            }
        return acquired_dict

    def get_point_mutations_display(self, fill_na="-"):  # used in per_sample_report
        """Convert point mutations to display format for per-sample report."""
        if not self.seq_variations:
            return {}
        db_lookup = {}
        if self.databases:
            for db in self.databases:
                db_lookup[db.key] = db.database_name
        point_dict = {}
        for variation in self.seq_variations:
            gene = variation.get_gene_name()
            if gene not in point_dict:
                point_dict[gene] = []
            nuc_change = (
                f"{variation.ref_codon.upper()}>{variation.var_codon.upper()}"
                if variation.ref_codon and variation.var_codon else None
            )
            pmids_display = variation.format_pmids_display()
            has_fungamr = variation.has_fungamr_notes()
            notes = None
            if has_fungamr and variation.notes:
                notes = (
                    "; ".join(str(note) for note in variation.notes)
                    if isinstance(variation.notes, list)
                    else str(variation.notes)
                )

            point_dict[gene].append({
                "gene": gene,
                "drugs": variation.phenotypes or [],
                "start": variation.ref_start_pos,
                "end": variation.ref_end_pos,
                "database": db_lookup.get(variation.ref_database, ""),
                "aa": variation.seq_var,
                "nuc": nuc_change,
                "pmids": pmids_display or fill_na,
                "has_fungamr_notes": has_fungamr,
                "notes": notes or fill_na
            })
        return point_dict

    def get_resistance_by_class(self):
        """Group resistance data by antibiotic class."""
        class_data = {}
        if not self.phenotypes:
            return class_data
        for phenotype in self.phenotypes:
            if not phenotype.amr_classes:
                continue
            # Skip phenotypes without associated data
            has_data = (phenotype.seq_regions or phenotype.seq_variations)
            if not has_data:
                continue
            for amr_class in phenotype.amr_classes:
                # Replace "Under_development" with "Other"
                display_class = (
                    "Other" if amr_class.lower() == "under_development"
                    else amr_class
                )
                if display_class not in class_data:
                    class_data[display_class] = {
                        'phenotypes': [],
                        'unique_drugs': set(),
                        'database_counts': {}
                    }
                class_data[display_class]['phenotypes'].append(phenotype)
                if phenotype.drug:
                    class_data[display_class]['unique_drugs'].add(phenotype.drug)
        # Count database occurrences - only for this specific abx class
        for amr_class, data in class_data.items():
            database_genes = {}      # db_name -> set of gene names
            database_mutations = {}  # db_name -> mutation count
            for phenotype in data['phenotypes']:
                self._count_genes_from_regions(
                    phenotype, database_genes
                )
                self._count_mutations_from_variations(
                    phenotype, database_mutations
                )

            for db_name, genes in database_genes.items():
                data['database_counts'][db_name] = (
                    data['database_counts'].get(db_name, 0) + len(genes))
            for db_name, mutation_keys in database_mutations.items():
                data['database_counts'][db_name] = (
                    data['database_counts'].get(db_name, 0) + len(mutation_keys))
        return class_data

    def _count_genes_from_regions(self, phenotype, database_genes):
        """Count unique genes from seq_regions by database."""
        if not (phenotype.seq_regions and self.seq_regions and self.databases):
            return
        if not phenotype.ref_database:
            return
        for region in self.seq_regions:
            if region.key not in phenotype.seq_regions:
                continue
            if not region.ref_database:
                continue
            gene_name = region.get_gene_name()
            for db_key in region.ref_database:
                if db_key not in phenotype.ref_database:
                    continue
                db_name = self._get_database_name(db_key)
                if not db_name:
                    continue
                if db_name not in database_genes:
                    database_genes[db_name] = set()
                database_genes[db_name].add(gene_name)

    def _count_mutations_from_variations(self, phenotype, database_mutations):
        """Count mutations from seq_variations by database."""
        if not (phenotype.seq_variations and self.seq_variations and self.databases):
            return
        for variation in self.seq_variations:
            if (
                variation.key not in phenotype.seq_variations
                or not variation.ref_database
            ):
                continue
            db_name = self._get_database_name(variation.ref_database)
            if db_name:
                # track unique mutation keys per database
                if db_name not in database_mutations:
                    database_mutations[db_name] = set()
                database_mutations[db_name].add(variation.key)

    def _get_database_name(self, db_key):
        """Get lowercase database name from database key."""
        for database in self.databases:
            if database.key == db_key:
                return database.database_name.lower()
        return None

    def get_genes_by_class(self, amr_class_phenotypes):
        """Group genes and variations by gene name for a specific class."""
        gene_data = {}
        for phenotype in amr_class_phenotypes:
            self._process_regions_for_genes(phenotype, gene_data)
            self._process_variations_for_genes(phenotype, gene_data)
        return gene_data

    def _process_regions_for_genes(self, phenotype, gene_data):
        """Process seq_regions and add to gene_data."""
        if not (phenotype.seq_regions and self.seq_regions):
            return
        for region in self.seq_regions:
            if region.key not in phenotype.seq_regions:
                continue
            gene_name = region.get_gene_name()
            if gene_name not in gene_data:
                gene_data[gene_name] = {
                    'region': region,
                    'drugs': set(),
                    'variations': {},
                    'has_acquired': False
                }
            if phenotype.drug:
                gene_data[gene_name]['drugs'].add(phenotype.drug)
            if region.phenotypes:
                gene_data[gene_name]['has_acquired'] = True

    def _process_variations_for_genes(self, phenotype, gene_data):
        """Process seq_variations and add to gene_data."""
        if not (phenotype.seq_variations and self.seq_variations):
            return
        for variation in self.seq_variations:
            if variation.key not in phenotype.seq_variations:
                continue
            if not (variation.seq_regions and self.seq_regions):
                continue
            for region in self.seq_regions:
                if region.key not in variation.seq_regions:
                    continue
                gene_name = region.get_gene_name()
                if gene_name not in gene_data:
                    gene_data[gene_name] = {
                        'region': region,
                        'drugs': set(),
                        'variations': {},
                        'has_acquired': False
                    }
                if phenotype.drug:
                    gene_data[gene_name]['drugs'].add(phenotype.drug)
                # Group variations by variation.key
                if variation.key not in gene_data[gene_name]['variations']:
                    gene_data[gene_name]['variations'][variation.key] = {
                        'variation': variation,
                        'drugs': set()
                    }
                if phenotype.drug:
                    gene_data[gene_name]['variations'][variation.key]['drugs'].add(
                        phenotype.drug
                    )


@dataclass
class MLST:
    """Multi-locus sequence typing results for the sample."""

    detected_species: str | None = field(
        default=None,
        metadata={
            "title": "Detected species",
            "description": """The detected species and MLST scheme
            used to classify the sample"""})
    sequence_type: str | None = field(
        default=None,
        metadata={
            "title": "Sequence type",
            "description": "The sequence type assigned to the sample"})
    typing_schema: list[SequenceTypeSchema] | None = field(
        default=None,
        metadata={
            "title": "Typing schema",
            "description": """The MLST schema alleles and variants
            identified in the sample"""})


@dataclass
class SourmashExcludedAssemblies:
    """Genomes excluded from sourmash database."""

    excluded_assemblies: list[str] | None = field(
        default_factory=list,
        metadata={
            "title": "Excluded assemblies",
            "description": "List of assembly IDs excluded from sourmash DB"
        })
    total_excluded: int | None = field(
        default=None,
        metadata={
            "title": "Total excluded count",
            "description": "Number of genomes excluded from the sourmash database"
        })


@dataclass
class SourmashMatch:
    """Individual sourmash match result."""

    rank: int | None = field(
        default=None,
        metadata={
            "title": "Gather result rank",
            "description": "Ranking of this match from sourmash gather"})
    query_containment_ani: float | None = field(
        default=None,
        metadata={
            "title": "Query containment ANI",
            "description": "Average nucleotide identity based on query containment"})
    reference_name: str | None = field(
        default=None,
        metadata={
            "title": "Reference genome name",
            "description": "Name of the reference genome that matched"})
    species_call: str | None = field(
        default=None,
        metadata={
            "title": "Species call",
            "description": "Species name extracted from taxonomic lineage"})
    intersect_bp: int | None = field(
        default=None,
        metadata={
            "title": "Intersect base pairs",
            "description": "Number of base pairs of overlap between query and match"})
    f_orig_query: float | None = field(
        default=None,
        metadata={
            "title": "Fraction of original query",
            "description": "Fraction of the query that matches this reference"})
    f_unique_to_query: float | None = field(
        default=None,
        metadata={
            "title": "Fraction unique to query",
            "description": "Fraction of the query that is unique to this match"})
    remaining_bp: int | None = field(
        default=None,
        metadata={
            "title": "Remaining base pairs",
            "description": "Base pairs remaining in query after this match"})


@dataclass
class SpeciesIdentification:
    """Sourmash species identification results."""

    detected_matches: list[SourmashMatch] | None = field(
        default=None,
        metadata={
            "title": "Detected matches",
            "description": "List of top sourmash matches"})


@dataclass
class Contig:
    """Summary statistics for contig in assembly."""

    name: str | None = field(
        default=None,
        metadata={
            "title": "Name",
            "description": "The name of the contig"})
    length: int | None = field(
        default=None,
        metadata={
            "title": "Contig length",
            "description": "The length of the contig in base pairs"})
    coverage: Coverage | None = field(
        default=None,
        metadata={
            "title": "Coverage",
            "description": "Mean coverage of the contig in the sample"})


@dataclass
class Assembly:
    """Draft genome assembly statistics of the sample."""

    reference: str | None = field(
        default=None,
        metadata={
            "title": "Reference name",
            "description": """Name of the reference used in the assembly
            process. Null for de-novo"""})
    annotations: list[Annotation] | None = field(
        default=None,
        metadata={
            "title": "Annotations",
            "description": """Array of regions of interest identified within
            the assembly"""})
    variants: list[Variant] | None = field(
        default=None,
        metadata={
            "title": "Variants",
            "description": "A list of variants identified in the assembly"})
    contig: list[Contig] | None = field(
        default=None,
        metadata={
            "title": "Contig",
            "description": "Summary statistics for each contig in the assembly"})
    flye_stats: list[FlyeContigStats] | None = field(
        default=None,
        metadata={
            "title": "Flye assembly statistics",
            "description": "Assembly statistics from Flye assembler"})

    def get_flye_summary(self):
        """Calculate summary statistics from Flye assembly data."""
        if not self.flye_stats:
            return {
                'contig_num': 0,
                'total_length': 0,
                'total_yield': '0 bp',
                'circular_num': 0
            }
        contig_num = len(self.flye_stats)
        total_length = sum(stat.length for stat in self.flye_stats if stat.length)
        total_yield = convert_bp(total_length)
        circular_num = sum(1 for stat in self.flye_stats if stat.circular == "Y")
        return {
            'contig_num': contig_num,
            'total_length': total_length,
            'total_yield': total_yield,
            'circular_num': circular_num
        }


class ContigType(str, Enum):
    """Enum for molecule types from MOB-suite."""

    PLASMID = "plasmid"
    CHROMOSOME = "chromosome"


@dataclass
class PlasmidID:
    """Individual plasmid contig information from MOB-suite."""

    contig_id: str | None = field(
        metadata={
            "title": "Contig ID",
            "description": "Identifier for the contig (extracted before space)"})
    contig_type: ContigType | None = field(
        default=None,
        metadata={
            "title": "Contig type",
            "description": "Classification as plasmid or chromosome"})
    primary_cluster_id: str | None = field(
        default=None,
        metadata={
            "title": "Primary cluster ID",
            "description": "Primary MOB cluster identifier"})
    size: int | None = field(
        default=None,
        metadata={
            "title": "Size",
            "description": "Size of the contig in base pairs"})
    rep_type: str | None = field(
        default=None,
        metadata={
            "title": "Replicon types",
            "description": "Replication types identified"})
    relaxase_type: str | None = field(
        default=None,
        metadata={
            "title": "Relaxase types",
            "description": "Relaxase types identified"})
    mpf_type: str | None = field(
        default=None,
        metadata={
            "title": "MPF type",
            "description": "Mating pair formation gene types"})
    orit_type: str | None = field(
        default=None,
        metadata={
            "title": "oriT types",
            "description": "Origin of transfer types"})
    predicted_mobility: str | None = field(
        default=None,
        metadata={
            "title": "Predicted mobility",
            "description": "Predicted mobility status"})
    mash_nearest_neighbor: str | None = field(
        default=None,
        metadata={
            "title": "Mash nearest neighbor",
            "description": "Closest matching reference"})
    mash_neighbor_distance: float | None = field(
        default=None,
        metadata={
            "title": "Mash neighbor distance",
            "description": "Distance to nearest neighbor"})
    mash_neighbor_identification: str | None = field(
        default=None,
        metadata={
            "title": "Mash neighbor identification",
            "description": "Identification of nearest neighbor"})
    predicted_host_range_rank: str | None = field(
        default=None,
        metadata={
            "title": "Predicted host range rank",
            "description": "Taxonomic rank of predicted host range"})
    predicted_host_range_name: str | None = field(
        default=None,
        metadata={
            "title": "Predicted host range name",
            "description": "Taxonomic name of predicted host range"})
    observed_host_range_rank: str | None = field(
        default=None,
        metadata={
            "title": "Observed host range rank",
            "description": "Taxonomic rank of observed host range"})
    observed_host_range_name: str | None = field(
        default=None,
        metadata={
            "title": "Observed host range name",
            "description": "Taxonomic name of observed host range"})


@dataclass
class PlasmidIdentification:
    """MOB-suite plasmid identification results."""

    detected_contigs: list[PlasmidID] = field(
        default_factory=list,
        metadata={
            "title": "Detected contigs",
            "description": "List of contigs identified by MOB-suite"})

    n_plasmid_contigs: int = field(
        default=0,
        metadata={
            "title": "Number of plasmid contigs",
            "description": "Count of contigs classified as plasmids"})

    n_chromosome_contigs: int = field(
        default=0,
        metadata={
            "title": "Number of chromosome contigs",
            "description": "Count of contigs classified as chromosomes"})


@dataclass
class ResultsContents(BaseResultsContents):
    """Results for a sample."""

    fastq: FastqStats | None = field(
        default=None,
        metadata={
            "title": "FASTQ Statistics",
            "description": "Read statistics from FASTQ processing"
        })
    assembly: Assembly | None = field(
        default=None,
        metadata={
            "title": "Assembly",
            "description": "Assembly statistics and information"
        })
    sequence_typing: MLST | None = field(
        default=None,
        metadata={
            "title": "Multi-locus sequence typing",
            "description": "MLST analysis results for species identification and typing"
        })
    sourmash_excluded_genomes: SourmashExcludedAssemblies | None = field(
        default=None,
        metadata={
            "title": "Sourmash excluded genomes",
            "description": "List of genomes excluded from sourmash database analysis"
        })
    plasmid_identification: PlasmidIdentification | None = field(
        default=None,
        metadata={
            "title": "Plasmid identification",
            "description": "MOB-suite plasmid identification results"
        })
    species_identification: SpeciesIdentification | None = field(
        default=None,
        metadata={
            "title": "Species identification",
            "description": "Taxonomic classification results from sourmash analysis"
        })
    antimicrobial_resistance: AntimicrobialResistance | None = field(
        default=None,
        metadata={
            "title": "Antimicrobial resistance",
            "description": "ResFinder analysis results for drug resistance prediction"
        })
    serotyping: Serotype | None = field(
        default=None,
        metadata={
            "title": "Serotyping",
            "description": "Salmonella serotype prediction results"
        })


@dataclass
class Sample(BaseSample):
    """A sample sheet entry and its corresponding checks and related results."""

    sample_pass: bool | None = field(
        default=None,
        metadata={
            "title": "Sample pass",
            "description": "If true the sample has passed workflow checks"})

    results: ResultsContents | None = field(
        default=None, metadata={
            "title": "Sample results",
            "description": """Results for the sample that are specific
            to the workflow in question."""})

    def has_assembly(self):
        """Check if assembly data is available."""
        return bool(
            self.results
            and self.results.assembly
            and self.results.assembly.flye_stats
        )

    def has_annotations(self):
        """Check if annotation data is available."""
        return bool(
            self.results
            and self.results.assembly
            and self.results.assembly.annotations
        )

    def has_mlst(self):
        """Check if MLST data is available."""
        if not (self.results and self.results.sequence_typing):
            return False
        species = self.results.sequence_typing.detected_species
        if not species:
            return False
        # Check it's not empty or a placeholder
        species_str = str(species).strip()
        return bool(species_str and species_str != "-" and species_str.lower() != "nan")

    def has_amr(self):
        """Check if AMR data is available."""
        return bool(
            self.results
            and self.results.antimicrobial_resistance
            and (
                self.results.antimicrobial_resistance.seq_regions
                or self.results.antimicrobial_resistance.seq_variations
            )
        )

    def has_serotype(self):
        """Check if serotype data is available."""
        return bool(
            self.results
            and self.results.serotyping
            and self.results.serotyping.predicted_serotype
        )

    def has_species_identification(self):
        """Check if species identification data is available."""
        return bool(
            self.results
            and self.results.species_identification
            and self.results.species_identification.detected_matches
        )

    # for summary report
    def get_summary_data(self, fill_na="-"):
        """Extract summary metrics from this sample for report display."""
        # Sample identifiers
        barcode = self.barcode or fill_na
        alias = self.alias or fill_na

        # Fastq metrics
        if self.results and self.results.fastq:
            n_seqs = self.results.fastq.n_seqs or 0
            median_len = self.results.fastq.median_length or 0
        else:
            n_seqs = 0
            median_len = 0

        # Assembly/contigs
        if self.results and self.results.assembly and self.results.assembly.contig:
            contig_count = len(self.results.assembly.contig)
        else:
            contig_count = 0

        # MLST
        if self.has_mlst():
            mlst = self.results.sequence_typing.detected_species or fill_na
        else:
            mlst = fill_na

        # Sourmash species identification
        if self.has_species_identification():
            sourmash_match = self.results.species_identification.detected_matches[0]
            species_call = sourmash_match.species_call or fill_na
            ani = sourmash_match.query_containment_ani
            f_query = sourmash_match.f_orig_query
        else:
            species_call = fill_na
            ani = None
            f_query = None

        return {
            'barcode': barcode,
            'alias': alias,
            'n_seqs': n_seqs,
            'median_len': median_len,
            'contig_count': contig_count,
            'mlst': mlst,
            'species_call': species_call,
            'ani': ani,
            'f_query': f_query
        }

    # for per_sample_report
    def get_run_summary(self, reference=None):
        """Get run summary statistics for sample."""
        run_dict = {}
        if self.results and self.results.fastq:
            total_yield = self.results.fastq.n_bases or 0
            median_read_length = self.results.fastq.median_length or 0
            median_read_quality = self.results.fastq.median_quality or 0
        else:
            total_yield = 0
            median_read_length = 0
            median_read_quality = 0
        if self.has_mlst():
            taxon = self.results.sequence_typing.detected_species
        else:
            taxon = "Unknown"
        if reference is None:
            run_dict["taxon"] = taxon.capitalize()
        else:
            run_dict["reference"] = os.path.basename(reference)
        run_dict["total_yield"] = convert_bp(total_yield)
        run_dict["median_read_length"] = convert_bp(median_read_length)
        run_dict["median_read_quality"] = f"{median_read_quality: .1f}"
        return run_dict

    def has_plasmids(self):
        """Check if plasmid identification data is available."""
        return bool(
            self.results
            and self.results.plasmid_identification
            and self.results.plasmid_identification.detected_contigs
            and any(
                contig.contig_type == ContigType.PLASMID
                for contig in self.results.plasmid_identification.detected_contigs
            )
        )

    def get_plasmid_summary(self):
        """Get plasmid identification summary statistics."""
        if not self.results or not self.results.plasmid_identification:
            return {
                'total_contigs': 0,
                'plasmid_count': 0,
                'chromosome_count': 0,
                'distinct_clusters': 0
            }
        plasmid_id = self.results.plasmid_identification
        # Count distinct clusters (only from plasmid contigs)
        plasmid_clusters = set()
        for contig in plasmid_id.detected_contigs:
            if (
                contig.contig_type == ContigType.PLASMID and
                contig.primary_cluster_id
            ):
                plasmid_clusters.add(contig.primary_cluster_id)
        return {
            'total_contigs': len(plasmid_id.detected_contigs),
            'plasmid_count': plasmid_id.n_plasmid_contigs,
            'chromosome_count': plasmid_id.n_chromosome_contigs,
            'distinct_clusters': len(plasmid_clusters)
        }

    def get_plasmid_table_data(self, fill_na="-"):
        """Get plasmid data formatted for table display."""
        if not self.has_plasmids():
            return pd.DataFrame()
        rows = []
        for contig in self.results.plasmid_identification.detected_contigs:
            # Only include plasmid contigs
            if contig.contig_type == ContigType.PLASMID:
                rows.append({
                    'Contig': contig.contig_id or fill_na,
                    'Length (bp)': contig.size or fill_na,
                    'MOB-cluster ID': contig.primary_cluster_id or fill_na,
                    'Replicon type': contig.rep_type or fill_na,
                    'Relaxase type': contig.relaxase_type or fill_na,
                    'MPF type': contig.mpf_type or fill_na,
                    'oriT type': contig.orit_type or fill_na,
                    'Mobility': contig.predicted_mobility or fill_na,
                })
        if not rows:
            return pd.DataFrame()
        df = pd.DataFrame(rows)
        # Placeholder fot 'null' cluster IDs for sorting
        df['cluster_sort'] = df['MOB-cluster ID'].fillna('zzz_no_cluster')
        # Sort by cluster ID (groups same plasmids), then by length descending
        sorted_df = df.sort_values(
            ['cluster_sort', 'Length (bp)'],
            ascending=[True, False]
        ).drop('cluster_sort', axis=1)
        return sorted_df

    def get_plasmid_clusters(self, fill_na="-"):
        """Get plasmids grouped by cluster with aggregated metadata."""
        if (
            not self.results or
            not self.results.plasmid_identification or
            not self.results.plasmid_identification.detected_contigs
        ):
            return {}
        plasmid_df = self.get_plasmid_table_data(fill_na)
        plasmid_objects = self.results.plasmid_identification.detected_contigs
        if plasmid_df.empty:
            return {}

        # Group plasmids by MOB-cluster ID
        plasmid_groups = {}
        for _, row in plasmid_df.iterrows():
            cluster_id = (
                row['MOB-cluster ID'] if row['MOB-cluster ID'] != fill_na
                else "Unassigned")
            if cluster_id not in plasmid_groups:
                plasmid_groups[cluster_id] = {
                    'rows': [],
                    'objects': [],
                    'metadata': {}
                }
            plasmid_groups[cluster_id]['rows'].append(row)
            contig_id = row['Contig']
            for obj in plasmid_objects:
                if obj.contig_id == contig_id and obj.contig_type == ContigType.PLASMID:
                    plasmid_groups[cluster_id]['objects'].append(obj)
                    break

        # Calculate metadata for each cluster
        for cluster_id, group_data in plasmid_groups.items():
            contigs = group_data['rows']
            contig_objects = group_data['objects']
            # Calculate totals
            num_contigs = len(contigs)
            total_length = sum(
                int(c['Length (bp)'])
                for c in contigs
                if c['Length (bp)'] != fill_na and str(c['Length (bp)']).isdigit()
            )
            # Summary information above table:
            # replicon type, mobility, mash neighbour, host range
            rep_types = {
                c['Replicon type'] for c in contigs if c['Replicon type'] != fill_na}
            mob_types = {
                c['Mobility'] for c in contigs if c['Mobility'] != fill_na}
            mash_neighbors = set()
            host_ranges = set()
            observed_host_ranges = set()
            for obj in contig_objects:
                if obj.mash_neighbor_identification:
                    neighbor_str = obj.mash_neighbor_identification
                    if obj.mash_neighbor_distance is not None:
                        neighbor_str += f" ({obj.mash_neighbor_distance: .4f})"
                    mash_neighbors.add(neighbor_str)
                if obj.predicted_host_range_rank and obj.predicted_host_range_name:
                    host_ranges.add(
                        f"{obj.predicted_host_range_rank} "
                        f"{obj.predicted_host_range_name}")
                if obj.observed_host_range_rank and obj.observed_host_range_name:
                    observed_host_ranges.add(
                        f"{obj.observed_host_range_rank} "
                        f"{obj.observed_host_range_name}")

            group_data['metadata'] = {
                'num_contigs': num_contigs,
                'total_length': total_length,
                'rep_types': rep_types,
                'mob_types': mob_types,
                'mash_neighbors': mash_neighbors,
                'host_ranges': host_ranges,
                'observed_host_ranges': observed_host_ranges
            }

        return plasmid_groups


@dataclass
class WorkflowResult(BaseWorkflowResult):
    """Definition for results that will be returned by this workflow."""

    samples: list[Sample] = field(
        metadata={
            "title": "Samples",
            "description": "Samples in this workflow instance"})

    def __post_init__(self):
        """Determine overall status for the workflow given the sample results."""
        self.workflow_pass = all(
            sample.sample_pass for sample in self.samples)
