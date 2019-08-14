# seqrenamer

Scripts for renaming sequence ids in fasta, gff3, tsv, and csv files.

The main use for these scripts is to deduplicate sequences,
giving them an unique id, and reduplicate the results after some
analysis has been done on the unique sequences.

Can also used to give things simple short ids for running software that panics
if it gets weird ids, or renaming ids given a tsv file mapping between
old and new ids.

Currently only a single type of id is available, which takes the form:
`{PREFIX}{ALPHANUM}` where `PREFIX` is user definable and alphanum is a
combination of integers and uppercase letters.

Eventually it'll support more arbitrary id substitution templates.


## Install

This is a python3 based package which should be installed with pip.
I recommend installing as user or inside a [virtual env](https://virtualenv.pypa.io/en/latest/).
 I haven't figured out how to get it working with conda environments yet.

For Linux, Mac, and WSL users using bash (or similar).

```bash
URL="https://github.com/darcyabjones/seqrenamer/archive/master.zip"

# To install as user.
python3 -m pip install --user "${URL}"

# Inside a virtualenv
python3 -m venv env
source env/bin/activate

env/bin/python3 -m pip install "${URL}"


# To install as root (not recommended)
sudo python3 -m pip install "${URL}"
```

Windows `CMD` and `powershell` are not officially supported (I have to way of testing it).
But something like this **might** work.

```cmd
Python -m pip install https://github.com/darcyabjones/seqrenamer/archive/master.zip
```

This will install the latest version of the code, and the master branch should generally be stable.
To install a specific version change the `master.zip` bit to whatever tag you like or clone the repo and call `pip install .` from inside that directory.


## Quick start

Seqrenamer is called using the command `sr` and has two sub-commands: `encode` and `decode`.
The `--help` parameter will print a list of commands and short descriptions for each subcommand.

`encode` replaces ids with a new one and can deduplicate sequences, outputting an "encoded" version of the input file and a tab-separated file mapping new ids to old ones.

```bash
# Rename fasta ids.
sr encode --map id_map.tsv my.fasta > encoded.fasta

# Remove duplicate sequences and strip descriptions for fasta headers.
sr encode --map id_map.tsv --outfile encoded.fasta --deduplicate --drop-desc my.fasta

# Take input from stdin. Strip "*" characters and change to uppercase (will affect deduplication).
  cat my.fasta \
| sr encode \
    --format fasta \
    --map id_map.tsv \
    --outfile encoded.fasta \
    --deduplicate \
    --drop-desc \
    --strip "*" \
    --upper \
    -
```

In each case, the file `encoded.fasta` will have the encoded fasta file,
and the file `id_map.tsv` will be a tab-separated file containing the columns: "old_id", "new_id", "description".
Where description is anything in the fasta header after the first space.
For `--deduplicate`-ed runs, an additional column containing the seguid checksum is included.


`decode` takes a file to be decoded and a map of old-ids to new ones, then outputs the results according to the id map.

```bash
sr decode \
  --map id_map.tsv \
  --outfile decoded_results.tsv \
  --header \
  encoded_results.tsv
```

The id map could be the output of a previous `encode` call, or it could be any tsv file with the encoded and decoded ids in the first two columns.
So you could create this in a spreadsheet, or manipulate tsvs using `cut`, `sed`, or `awk` to get arbitrary recoding.

Note that if an "encoded" id is present in the first column multiple times (e.g. from deduplicated sequences), the output will duplicate the result for each corresponding "decoded" id from the map file.
The exception to this reduplication is the `gff3` format with the `id` column, where only a 1-1 mapping is supported.
This is because both the `ID` and `Parent` attribute fields must be updated and a record may have multiple `Parent`s, so to reduplicate them the output would have to be combinatorial which doesn't make sense.
