"""Parsers for output of `bcftools stats` and `Bakta`."""

from collections import defaultdict
import os
import re

import pandas as pd


def format_gff3_string(text):
    """Decode URL-encoded characters for display."""
    text = text.replace("%2C", ",")
    text = text.replace("%3B", ";")
    text = text.replace("%3D", "=")
    text = text.replace("%26", "&")
    text = text.replace("%25", "%")
    text = text.replace("%09", " ")  # tab → space
    text = text.replace("%0A", " ")  # new line → space
    text = text.replace("%0D", " ")  # carriage return → space
    return text


def parse_bakta_gff3(filename):
    """Process a Bakta `.gff3` output file."""
    data = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                # skip empty lines
                continue
            if line.startswith("##"):
                if line == "##FASTA":
                    # only FASTA seqs from now on --> stop processing the file
                    break
                continue
            if line.startswith("#"):
                # skip comment lines in the file header
                continue

            contig, _, _, start, end, _, strand, _, attrs = line.split("\t")
            name, prod, ec_number, gene = "-", "-", "-", "-"
            for attr in attrs.split(";"):
                if attr.startswith("ID"):
                    name = attr.split("=")[1]
                elif attr.startswith("product"):
                    prod = format_gff3_string(attr.split("=")[1])
                elif attr.startswith("ec_number"):
                    ec_number = attr.split("=")[1]
                elif attr.startswith("gene"):
                    gene = attr.split("=")[1]
            data.append([contig, name, start, end, strand, gene, prod, ec_number])

    df = pd.DataFrame(
        data,
        columns=[
            "contig",
            "ID",
            "start",
            "end",
            "strand",
            "gene",
            "product",
            "EC number",
        ],
    )
    return df


def split_bcftools_stats_blocks(fname):
    """Split lines of a file into blocks based on comment lines."""
    comment = list()
    data = list()
    state = None
    with open(fname, "r") as fh:
        # TODO: handle empty tables
        for line in fh.readlines():
            if line.startswith("#"):
                if state != "comment" and state is not None:
                    yield comment, data
                    comment, data = list(), list()
                state = "comment"
                comment.append(line.strip("# ").rstrip())
            else:
                state = "data"
                data.append(line.rstrip())
        yield comment, data


def parse_bcftools_stats(fname):
    """Parse `bcftools stats` output.

    :param fname: file to parse.
    """
    tables = dict()
    for comment, data in split_bcftools_stats_blocks(fname):
        fields = [x.rstrip() for x in re.split(r"\[\d+\]", comment[-1])]
        section, fields = fields[0], fields[1:]
        rows = list()
        for d in data:
            items = d.split("\t")
            if items[0] != section:
                raise ValueError("first data field not equal to section key")
            rows.append(items[1:])
        tables[section] = pd.DataFrame(rows, columns=fields)
    # now some special handling
    sn_df = tables["SN"]
    sn_df["key"] = sn_df["key"].apply(lambda x: x.replace("number of ", "").rstrip(":"))
    tables["SN"] = pd.DataFrame(
        sn_df.pivot(index="id", columns="key", values="value").to_records()
    )
    return tables


def parse_bcftools_stats_multi(fnames, sample_names=None):
    """Parse multiple bcf stats outputs and combine.

    :param fnames: list of filenames output by `bcftools stats.`
    :param sample_names: list of names of each sample, else the basename
        of each input file is used.
    """
    if sample_names is not None:
        if len(sample_names) != len(fnames):
            raise TypeError("`fnames` and `sample_names` should be of equal length.")
    else:
        sample_names = [os.path.basename(f) for f in fnames]
    all_tables = defaultdict(list)
    for sample, fname in zip(sample_names, fnames):
        tables = parse_bcftools_stats(fname)
        # we don't use the 'id' field - it's for intersections
        # and unions of two VCFs given to bcftools stats
        for k, t in tables.items():
            t.drop(columns="id", inplace=True)
            t.insert(0, "sample", sample)
            all_tables[k].append(t)
    for k, t in all_tables.items():
        all_tables[k] = pd.concat(t)
    return dict(all_tables)  # don't want defaultdict


def parse_mlst(filename):
    """Parse mlst json to pandas df."""
    columns = ["id", "scheme", "sequence_type", "alleles"]
    df = pd.read_json(filename)
    df = df[columns]
    # Catches failed MLST run (i.e df["alleles"] is NaN)
    # Allele information is nested in second layer of json so needs to be extracted and
    # Added into dataframe
    if df["alleles"].isnull().values.any():
        return None
    else:
        alleles = pd.DataFrame(df["alleles"].tolist()).sort_index(axis=1)
        df = pd.concat([df.drop("alleles", axis=1), alleles], axis=1)
        df.columns = df.columns.str.replace("_", " ")
    return df
