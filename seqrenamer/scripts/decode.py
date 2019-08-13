import sys
import argparse

from seqrenamer.exceptions import MapFileParseError
from seqrenamer.exceptions import MapFileKeyError
from seqrenamer.seq import Seq, Seqs
from seqrenamer.scripts.encode import check_format, check_column
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


def decode_seqs(infiles, outfile, map_, column):
    seqs = Seqs.parse_many(infiles)

    chunk = list()
    for i, seq in enumerate(seqs, 1):
        if column == "description":
            old_desc = get_from_map(map_, seq.desc)
            chunk.append(str(Seq(seq.id, old_desc, seq.seq)))
        else:
            old_id = get_from_map(map_, seq.id)
            chunk.append(str(Seq(old_id, seq.desc, seq.seq)))

        if i % 10000:
            outfile.write(''.join(chunk))
            chunk = []

    outfile.write(''.join(chunk))
    return


def decode(args):
    format = check_format(args)
    column = check_column(args, format)
    map_ = parse_map_file(args.map)

    if format == "fasta":
        decode_seqs(args.infiles, args.outfile, map_, column)
    else:
        raise ValueError("This shouldn't ever happen")
    return
