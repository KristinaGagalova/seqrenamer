import sys
import argparse

import csv

from gffpal.gff import GFFRecord

from seqrenamer.exceptions import MapFileParseError
from seqrenamer.exceptions import MapFileKeyError
from seqrenamer.exceptions import XsvColumnNumberError
from seqrenamer.seq import Seq, Seqs
from seqrenamer.xsv import Xsv
from seqrenamer.scripts.encode import check_format, check_column, join_files
from seqrenamer.scripts.encode import FORMATS


def cli_decode(parser):

    parser.add_argument(
        "-f", "--format",
        choices=["auto"] + FORMATS,
        default="auto",
        help=(
            "The format of the infile. "
            "By default will attempt to determine format by file extension."
        )
    )

    parser.add_argument(
        "-c", "--column",
        default=None,
        help=(
            "Which column or field has the id that you want substituted. "
            "For tsv and csv, a 0-based column index can be used (Default 0). "
            "For fasta, either 'id' (Default) or 'description' may be used. "
            "For gff the columns 'seqid', 'id' (Default), 'name' may be used. "
            "If the gff 'id' is used, 'Parent' ids will be updated too."
        )
    )
    parser.add_argument(
        "-H", "--header",
        default=False,
        action="store_true",
        help=(
            "The first line is a header (skip it). "
            "This is only respected for 'csv' and 'tsv' formats."
        ),
    )

    parser.add_argument(
        "-C", "--comment",
        default="#",
        help=(
            "Comment signifier. If this string is encountered at the "
            "start of the line, the line is printed as is."
        ),
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the output to.",
    )

    parser.add_argument(
        "-m", "--map",
        type=argparse.FileType("r"),
        required=True,
        help="The id mapping file to substitute ids from."
    )

    parser.add_argument(
        "infiles",
        metavar="INFILE",
        nargs="+",
        type=argparse.FileType('r'),
        help="The input file to decode.",
    )

    return


def parse_map_file(infile):
    """ Parse a map file to use for renaming outputs. """

    from collections import defaultdict
    map_ = defaultdict(list)

    for i, line in enumerate(infile, 1):
        try:
            sline = line.strip().split("\t")
            map_[sline[0]].append(sline[1])
        except IndexError:
            raise MapFileParseError(
                f"Error parsing map file at line {i}. "
                f"The offending line was: {line}"
            )

    return map_


def get_from_map(map, key):
    try:
        return map[key]
    except KeyError:
        raise MapFileKeyError(
            f"Key {key} is not in map file. "
            "Have you selected the right column?"
        )


def get_from_map_single(map, key):
    try:
        return get_from_map(map, key)[0]
    except IndexError:
        raise MapFileKeyError(
            f"Key {key} is has multiple substitution candidates. "
            "For this format we expect only one-to-one mapping."
        )


def decode_seqs(infiles, outfile, map_, column):
    seqs = Seqs.parse_many(infiles)

    chunk = list()
    for i, seq in enumerate(seqs, 1):
        if column == "description":
            old_descs = get_from_map(map_, seq.desc)
            for old_desc in old_descs:
                chunk.append(str(Seq(seq.id, old_desc, seq.seq)))
        else:
            old_ids = get_from_map(map_, seq.id)
            for old_id in old_ids:
                chunk.append(str(Seq(old_id, seq.desc, seq.seq)))

        if i % 10000 == 0:
            outfile.write(''.join(chunk))
            chunk = []

    outfile.write(''.join(chunk))
    return


def decode_xsv(
    infiles,
    outfile,
    map_,
    column,
    comment,
    header,
    sep,
):
    xsv_writer = csv.writer(outfile, delimiter=sep, dialect='excel')

    inhandles = join_files(infiles, header)
    xsv_reader = Xsv(inhandles, comment, sep)

    first = True
    for row in xsv_reader:
        if header and first:
            xsv_writer.writerow(row)
            first = False
            continue

        try:
            old_val = get_from_map_single(map_, row[column])
            row[column] = old_val
        except IndexError:
            joined_line = sep.join(map(str, row))
            raise XsvColumnNumberError(
                f"Could not access column '{column}' in a line. "
                f"The offending line was: {joined_line}."
            )

        xsv_writer.writerow(row)
    return


def replace_gff_id(record, map_):
    old_id = record.attributes.id
    old_parents = record.attributes.parent

    if old_id is not None:
        new_id = get_from_map_single(map_, old_id)
        record.attributes.id = new_id

    new_parents = []
    for old_parent in old_parents:
        new_parent = get_from_map_single(map_, old_parent)
        new_parents.append(new_parent)

    record.attributes.parent = new_parents
    return record


def replace_gff_name(record, map_):
    old_name = record.attributes.name

    if old_name is not None:
        new_name = get_from_map_single(map_, old_name)
        record.attributes.name = new_name

    return record


def replace_gff_seqid(record, map_):
    old_seqid = record.seqid
    new_seqid = get_from_map_single(map_, old_seqid)
    record.seqid = new_seqid
    return record


def decode_gff(infiles, outfile, map_, column):
    inhandles = join_files(infiles, header=False)

    if column == "id":
        trans_function = replace_gff_id
    elif column == "name":
        trans_function = replace_gff_name
    elif column == "seqid":
        trans_function = replace_gff_seqid
    else:
        raise ValueError("This shouldn't ever happen")

    record_chunk = list()
    for i, line in enumerate(inhandles):
        if line.startswith("#"):
            record_chunk.append(line.strip())
            continue

        old_record = GFFRecord.parse(line)
        new_record = trans_function(old_record, map_)
        record_chunk.append(str(new_record))

        if i % 10000 == 0:
            outfile.write("\n".join(record_chunk))
            record_chunk = list()

    if len(record_chunk) > 0:
        outfile.write("\n".join(record_chunk))

    return


def decode(args):
    format = check_format(args.format, args.infiles)
    column = check_column(args.column, format)

    map_ = parse_map_file(args.map)

    if format == "fasta":
        decode_seqs(args.infiles, args.outfile, map_, column)
    elif format in ("csv", "tsv"):
        decode_xsv(
            args.infiles,
            args.outfile,
            map_,
            column,
            args.comment,
            args.header,
            ',' if format == "csv" else '\t',
        )
    elif format == "gff3":
        decode_gff(
            args.infiles,
            args.outfile,
            map_,
            column,
        )
    else:
        raise ValueError("This shouldn't ever happen")
    return
