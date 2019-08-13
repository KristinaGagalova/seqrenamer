import sys
import argparse

from os.path import splitext

from seqrenamer.seq import Seqs
from seqrenamer.id_generator import IdConverter


FORMATS = ["fasta", "tsv", "csv", "gff3"]
EXTENSIONS = {
    "fasta": "fasta",
    "faa": "fasta",
    "fna": "fasta",
    "fas": "fasta",
    "fa": "fasta",
    "tsv": "tsv",
    "tab": "tsv",
    "csv": "csv",
    "gff3": "gff3",
    "gff": "gff3",
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
        type="str",
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
        help="The input file to encode.",
    )

    return


def check_column(args, format=None):  # noqa
    if format is None:
        format = args.format

    if format in ("tsv", "csv"):
        if args.column is None:
            return 0

        try:
            column = int(args.column)
            assert column >= 0
            return column
        except ValueError:
            args.error(
                "For tsv and csv formats a 0 based index must be used. "
                f"Could not parse {args.column} as an integer."
            )
        except AssertionError:
            args.error("Column indices must be >= 0.")

    elif format == "fasta":
        if args.column is None:
            return "id"

        if args.column not in ("id", "description"):
            args.error("For fasta format the column must "
                       "be 'id' or 'description'.")
        else:
            return args.column
    elif format == "gff3":
        if args.column is None:
            return "id"

        if args.column not in ("seqid", "id", "name"):
            args.error("For gff3 format the column must be "
                       "'seqid', 'id', or 'name'.")
        else:
            return args.column
    else:
        raise ValueError("This shouldn't ever happen.")


def check_format(args):
    """ Figures out what the format is if format is auto. """

    if args.format != "auto":
        return args.format

    msg = (
        "Could not automatically determine the format. "
        "Please specify the --format parameter."
    )

    try:
        filename = args.infiles[0].name
    except AttributeError:
        args.error(msg)

    extension = splitext(filename)
    format = EXTENSIONS.get(extension, None)

    if format is None:
        args.error(msg)

    return format


def encode_seqs(
    infiles,
    outfile,
    mapfile,
    column,
    deduplicate,
    upper,
    strip,
    id_conv
):
    seqs = Seqs.parse_many(infiles)

    if strip is not None:
        seqs = seqs.map_seq(lambda s: s.rstrip(strip))

    if upper:
        seqs = seqs.map_seq(lambda s: s.upper())

    if deduplicate:
        seqs = seqs.deduplicated(lambda i: next(id_conv), column=column)
    else:
        seqs = seqs.replace_ids(lambda i: next(id_conv), column=column)

    chunk = list()
    for i, seq in enumerate(seqs, 1):
        chunk.append(str(seq))

        if i % 10000:
            outfile.write(''.join(chunk))
            chunk = []
            seqs.flush_ids(mapfile)

    outfile.write(''.join(chunk))
    seqs.flush_ids(mapfile)
    return


def encode(args):
    format = check_format(args)
    column = check_column(args, format)

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
            id_conv
        )
    else:
        raise ValueError("This shouldn't happen")

    return
