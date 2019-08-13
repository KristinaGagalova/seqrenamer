import sys
import argparse

import csv
from os.path import splitext

from seqrenamer.seq import Seqs
from seqrenamer.id_generator import IdConverter
from seqrenamer.exceptions import InvalidArgumentError
from seqrenamer.xsv import Xsv


FORMATS = ["fasta", "tsv", "csv", "gff3"]
EXTENSIONS = {
    ".fasta": "fasta",
    ".faa": "fasta",
    ".fna": "fasta",
    ".fas": "fasta",
    ".fa": "fasta",
    ".tsv": "tsv",
    ".tab": "tsv",
    ".csv": "csv",
    ".gff3": "gff3",
    ".gff": "gff3",
}


def cli_encode(parser):

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
        "-l", "--length",
        default=5,
        type=int,
        help="The number of characters to use in the new id."
    )

    parser.add_argument(
        "-p", "--prefix",
        default="SR",
        type=str,
        help="The prefix to add to the beginning of the ids.",
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
        "-d", "--deduplicate",
        default=False,
        action="store_true",
        help=(
            "Remove duplicate sequences based on a checksum. "
            "This is only respected for 'fasta' formats."
        ),
    )

    parser.add_argument(
        "-s", "--strip",
        default=None,
        type=str,
        help=(
            "Strip these characters from the end of the sequence. "
            "This is only respected for 'fasta' formats."
        ),
    )

    parser.add_argument(
        "-U", "--upper",
        default=False,
        action="store_true",
        help=(
            "Convert sequences to uppercase. "
            "This is only used for 'fasta' format."
        ),
    )

    parser.add_argument(
        "--drop-desc",
        default=False,
        action="store_true",
        help=(
            "Remove descriptions from fasta headers. "
            "This is only used for 'fasta' format."
        ),
    )

    parser.add_argument(
        "-C", "--comment",
        default="#",
        help=(
            "Comment signifier. If this string is encountered at the "
            "start of the line, the line is excluded from output."
        ),
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="The ffindex .ffindex file.",
    )

    parser.add_argument(
        "-m", "--map",
        type=argparse.FileType("w"),
        required=True,
        help="Where to save the id mapping file."
    )

    parser.add_argument(
        "infiles",
        metavar="INFILE",
        nargs="+",
        type=argparse.FileType('r'),
        help="The input files to encode.",
    )

    return


def check_column(column, format):
    """ Gets the column to replace ids in.

    Examples:
    >>> check_column(0, "tsv")
    0
    >>> try:
    ...     check_column(-1, "csv")  # This should fail.
    ...     assert False  # We shouldn't reach this.
    ... except InvalidArgumentError:
    ...     pass
    >>> check_column(None, "csv")
    0
    >>> check_column(None, "fasta")
    'id'
    >>> check_column("id", "fasta")
    'id'
    >>> check_column("description", "fasta")
    'description'
    >>> try:
    ...     check_column("desc", "fasta")  # This should fail.
    ...     assert False  # We shouldn't reach this.
    ... except InvalidArgumentError:
    ...     pass
    >>> check_column(None, "gff3")
    'id'
    >>> check_column("id", "gff3")
    'id'
    >>> check_column("name", "gff3")
    'name'
    >>> check_column("seqid", "gff3")
    'seqid'
    >>> try:
    ...     check_column("type", "gff3")  # This should fail.
    ...     assert False  # We shouldn't reach this.
    ... except InvalidArgumentError:
    ...     pass
    """

    if format in ("tsv", "csv"):
        if column is None:
            return 0

        try:
            column = int(column)
            assert column >= 0
            return column
        except ValueError:
            raise InvalidArgumentError(
                "For tsv and csv formats a 0 based index must be used. "
                f"Could not parse {column} as an integer."
            )
        except AssertionError:
            raise InvalidArgumentError("Column indices must be >= 0.")

    elif format == "fasta":
        if column is None:
            return "id"

        if column in ("id", "description"):
            return column
        else:
            raise InvalidArgumentError(
                "For fasta format the column must "
                "be 'id' or 'description'."
            )

    elif format == "gff3":
        if column is None:
            return "id"

        if column in ("seqid", "id", "name"):
            return column
        else:
            raise InvalidArgumentError(
                "For gff3 format the column must be "
                "'seqid', 'id', or 'name'."
            )

    else:
        raise ValueError("This shouldn't ever happen.")


def check_format(format, infiles):
    """ Figures out what the format is if format is auto.

    Examples:
    >>> check_format("fasta", [None])
    'fasta'
    >>> class DummyFile:
    ...     def __init__(self, name):
    ...         self.name = name
    >>> check_format("auto", [DummyFile("test.fasta")])
    'fasta'
    >>> check_format("auto", [DummyFile("test.csv")])
    'csv'
    >>> check_format("auto", [DummyFile("test.csv"), DummyFile("test.tsv")])
    'csv'
    >>> try:
    ...     check_format("auto", [DummyFile("test")])  # No extension fail.
    ...     assert False  # We shouldn't reach this point.
    ... except InvalidArgumentError:
    ...     pass
    """

    if format != "auto":
        return format

    msg = (
        "Could not automatically determine the format. "
        "Please specify the --format parameter."
    )

    try:
        filename = infiles[0].name
    except AttributeError:
        raise InvalidArgumentError(msg)

    extension = splitext(filename)
    format = EXTENSIONS.get(extension[1], None)

    if format is None:
        raise InvalidArgumentError(msg)

    return format


def encode_seqs(
    infiles,
    outfile,
    mapfile,
    column,
    deduplicate,
    upper,
    strip,
    drop_desc,
    id_conv
):
    seqs = Seqs.parse_many(infiles)

    if strip is not None:
        seqs = seqs.map_seq(lambda s: s.rstrip(strip.encode()))

    if upper:
        seqs = seqs.map_seq(lambda s: s.upper())

    if deduplicate:
        seqs = seqs.deduplicated(lambda i: next(id_conv), column=column)
    else:
        seqs = seqs.replace_ids(lambda i: next(id_conv), column=column)

    chunk = list()
    for i, seq in enumerate(seqs, 1):
        if drop_desc:
            seq.desc = None

        chunk.append(str(seq))

        if i % 10000:
            outfile.write(''.join(chunk))
            chunk = []
            seqs.flush_ids(mapfile)

    outfile.write(''.join(chunk))
    seqs.flush_ids(mapfile)
    return


def join_files(infiles, header=False):
    first = True

    for f in infiles:
        if header and not first:
            # Drop the first line
            _ = next(f)

        first = False
        for line in f:
            yield line.rstrip("\r\n")
    return


def encode_xsv(
    infiles,
    outfile,
    mapfile,
    column,
    comment,
    header,
    sep,
    id_conv
):
    xsv_writer = csv.writer(outfile, delimiter=sep, dialect='excel')

    inhandles = join_files(infiles, header)
    xsv_reader = Xsv(inhandles, comment, sep)
    iterator = xsv_reader.replace_ids(lambda r: next(id_conv), column, header)

    for i, row in enumerate(iterator, 1):
        xsv_writer.writerow(row)

        if i % 10000:
            xsv_reader.flush_ids(mapfile)

    xsv_reader.flush_ids(mapfile)
    return


def encode(args):
    format = check_format(args.format, args.infiles)
    column = check_column(args.column, format)

    id_conv = IdConverter(prefix=args.prefix, length=args.length)

    if format == "fasta":
        encode_seqs(
            args.infiles,
            args.outfile,
            args.map,
            column,
            args.deduplicate,
            args.upper,
            args.strip,
            args.drop_desc,
            id_conv
        )
    elif format in ("csv", "tsv"):
        encode_xsv(
            args.infiles,
            args.outfile,
            args.map,
            column,
            args.comment,
            args.header,
            ',' if format == "csv" else '\t',
            id_conv,
        )
    else:
        raise ValueError("This shouldn't happen")

    return
