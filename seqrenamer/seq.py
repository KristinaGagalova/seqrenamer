""" Simple fasta parser and utilities. """

from collections.abc import Iterator


class Seq(object):

    def __init__(self, id, desc, seq):
        """ Construct a new Seq object.

        Keyword arguments:
        id -- The sequence id <str>.
        desc -- A short description of the sequence <str>.
        seq -- The biological sequence <str>.

        Examples:

        >>> Seq("test", "description", "ATGCA")
        Seq(id='test', desc='description', seq='b'ATGCA'')
        """
        self.id = id
        self.desc = desc

        # Seq should be bytes
        if isinstance(seq, str):
            seq = seq.encode()

        self.seq = seq
        return

    def __str__(self):
        """ Returns a FASTA string from the object.

        Examples:

        >>> print(Seq("test", "description", "ATGCA".encode()))
        >test description
        ATGCA
        <BLANKLINE>
        """

        line_length = 60

        if self.desc is None:
            lines = [">{}".format(self.id)]
        else:
            lines = [">{} {}".format(self.id, self.desc)]

        for i in range(0, len(self), line_length):
            lines.append(self.seq[i:i+line_length].decode("utf-8"))

        return "\n".join(lines) + "\n"

    def __repr__(self):
        """ Returns a simple string representation of the object. """
        cls = self.__class__.__name__
        return "{}(id='{}', desc='{}', seq='{}')".format(cls, self.id,
                                                         self.desc, self.seq)

    def __getitem__(self, key):
        """ Allow us to access indices from the seq directly.

        Examples:
        >>> seq = Seq("test", "description", "ATGCA".encode())
        >>> seq[0]
        Seq(id='test', desc='description', seq='b'A'')
        """

        # Keep the output as a byteslice even when selecting single character.
        if isinstance(key, int):
            key = slice(key, key + 1)

        seq = self.seq[key]
        return self.__class__(self.id, self.desc, seq)

    def __eq__(self, other):
        """ Allows us to compare two Seq objects directly using '==' .
        NB. python internally implements != based on this too.

        Examples:
        >>> seq = Seq("test", "description", "ATGCA".encode())
        >>> assert seq == Seq("test2", "description", "ATGCA".encode())
        >>> assert seq == "ATGCA"
        """

        if isinstance(other, self.__class__):
            return self.seq == other.seq
        elif isinstance(other, bytes):
            return self.seq == other
        elif isinstance(other, str):
            return self.seq == other.encode()
        else:
            raise ValueError((
                "Equality comparisons not implemented between {} and {}."
                ).format(type(self), type(other)))
        return

    def __len__(self):
        """ The length of a Seq object should be the length of the seq.

        Examples:
        >>> seq = Seq("test", "description", "ATGCA")
        >>> len(seq)
        5
        """

        return len(self.seq)

    @classmethod
    def read(cls, handle):
        """ Read a single FASTA record.
        Parses a single FASTA record into a Seq object.
        Assumes that the first line will always contain the id line,
        and that there is a single FASTA sequence.

        Keyword arguments:
        handle -- An iterable containing lines (newline split) of the file.

        Returns:
        A Seq object.

        Examples:
        >>> fasta = [">test description", "ATGCA"]
        >>> Seq.read(fasta)
        Seq(id='test', desc='description', seq='b'ATGCA'')
        """
        if not isinstance(handle, Iterator):
            handle = iter(handle)

        try:
            id_, desc = cls._split_id_line(next(handle).strip())
        except ValueError as e:
            raise ValueError("Fasta parsing failed. " + str(e))

        # tuple comprehensions are generators so we're still doing lazy eval
        seq = "".join((l.strip() for l in handle))
        return Seq(id_, desc, seq.encode())

    @classmethod
    def parse(cls, handle, comment=";"):
        """ Parse multiple fasta records.
        Parses a multi-fasta formatted file-like object.

        Keyword arguments:
        handle -- A file-like object or any iterable over the fasta file lines.

        Yields:
        Seq objects.

        Examples:
        >>> fasta = [
        ...     ">test1 description",
        ...     "ATGCA",
        ...     ">test2 descr",
        ...     "TGACA",
        ... ]
        >>> seqs = Seq.parse(fasta)
        >>> next(seqs)
        Seq(id='test1', desc='description', seq='b'ATGCA'')
        >>> next(seqs)
        Seq(id='test2', desc='descr', seq='b'TGACA'')
        """

        # Store the initial state to avoid outputting empty record.
        first = True
        # Store lines for this block here.
        current_record = []

        for line in handle:
            if line.startswith(">"):
                if not first:
                    # Yield makes this function a generator.
                    # NB we reuse the read method to avoid repetition.
                    # It's also easier to test.
                    yield cls.read(current_record)
                else:
                    # Once we've passed the first sequence this passes.
                    first = False

                # Start a new block
                current_record = [line]
            elif line.startswith(comment):
                continue
            else:
                # For lines containing sequences we simply append the sequence.
                current_record.append(line)

        # The last sequence in the file won't have a ">" following it.
        # so we yield the last block too.
        yield cls.read(current_record)
        return

    @classmethod
    def parse_many(cls, handles, comment=";"):
        """ Parses many files yielding an iterator over all of them. """

        for handle in handles:
            for record in cls.parse(handle, comment=";"):
                yield record
        return

    @staticmethod
    def _split_id_line(line):
        """ Parse the FASTA header line into id and description components.
        NB expects the '>' character to be present at start of line.

        Keyword arguments:
        line -- A string containing the header.

        Returns:
        Tuple -- id and description strings.

        Examples:
        >>> Seq._split_id_line(">one two")
        ('one', 'two')
        >>> Seq._split_id_line(">one ")
        ('one', '')
        >>> Seq._split_id_line(">one")
        ('one', None)
        """

        if not line.startswith(">"):
            raise ValueError(("Encountered malformed fasta header. "
                              "Offending line is '{}'").format(line))
        # Strip the ">" character and split at most 1 time on spaces.
        sline = line[1:].split(" ", 1)

        if len(sline) == 1:
            return sline[0], None
        else:
            return sline[0], sline[1]

        # We should never reach this point.
        return

    def checksum(self):
        """ Returns the seguid checksum of a sequence. """

        from hashlib import sha1
        from base64 import b64encode
        hash_ = sha1(self.seq).digest()
        return b64encode(hash_).rstrip(b"=").decode("utf-8")

    def rstrip(self, chars):
        """ Strips some bytes from the end of the sequence.

        Examples:
        >>> seq = Seq("test", None, "MAGNIFIQUE*")
        >>> seq.rstrip(b"*")
        Seq(id='test', desc='None', seq='b'MAGNIFIQUE'')
        """

        return self.__class__(self.id, self.desc, self.seq.rstrip(chars))

    def upper(self):
        """ Converts the sequence to all uppercase.

        Examples:
        >>> seq = Seq('test', None, "AtGca")
        >>> seq.upper()
        Seq(id='test', desc='None', seq='b'ATGCA'')
        """

        return self.__class__(self.id, self.desc, self.seq.upper())


class Seqs(object):

    def __init__(self, seqs):
        self.seqs = seqs
        return

    @classmethod
    def parse(cls, handle, comment=";"):
        """ Parse multiple fasta records.
        Parses a multi-fasta formatted file-like object.

        Keyword arguments:
        handle -- A file-like object or any iterable over the fasta file lines.

        Returns:
        A Seqs object containing a generator of Seq objects.

        Examples:
        >>> fasta = [
        ...     ">test1 description",
        ...     "ATGCA",
        ...     ">test2 descr",
        ...     "TGACA",
        ... ]
        >>> seqs = iter(Seqs.parse(fasta))
        >>> next(seqs)
        Seq(id='test1', desc='description', seq='b'ATGCA'')
        >>> next(seqs)
        Seq(id='test2', desc='descr', seq='b'TGACA'')
        """
        return cls(Seq.parse(handle, comment=comment))

    @classmethod
    def parse_many(cls, handles, comment=";"):
        return cls(Seq.parse_many(handles, comment=comment))

    def filter(self, function):
        """ Filters out sequences if a function returns False.

        Examples:
        >>> fasta = [
        ...     ">test1 description",
        ...     "ATGCA",
        ...     ">test2 descr",
        ...     "TGACA",
        ... ]
        >>> seqs = Seqs.parse(fasta)
        >>> seqs = seqs.filter(lambda s: s.id == 'test2')
        >>> seqs = iter(seqs)
        >>> next(seqs)
        Seq(id='test2', desc='descr', seq='b'TGACA'')
        """

        return self.__class__(s for s in self.seqs if function(s))

    def map(self, function):
        """ Apply a function to all sequences.

        Examples:
        >>> fasta = [
        ...     ">test1 description",
        ...     "atgca",
        ...     ">test2 descr",
        ...     "TGACA",
        ... ]
        >>> seqs = Seqs.parse(fasta).map(lambda s: s.upper())
        >>> seqs = iter(seqs)
        >>> next(seqs)
        Seq(id='test1', desc='description', seq='b'ATGCA'')
        >>> next(seqs)
        Seq(id='test2', desc='descr', seq='b'TGACA'')
        """

        return self.__class__(function(s) for s in self.seqs)

    def map_id(self, function):
        """ Apply a function to all sequence ids.

        Examples:
        >>> fasta = [
        ...     ">test1 description",
        ...     "ATGCA",
        ...     ">test2 descr",
        ...     "TGACA",
        ... ]
        >>> seqs = Seqs.parse(fasta).map_id(lambda id: "new_id")
        >>> seqs = iter(seqs)
        >>> next(seqs)
        Seq(id='new_id', desc='description', seq='b'ATGCA'')
        >>> next(seqs)
        Seq(id='new_id', desc='descr', seq='b'TGACA'')
        """

        return self.__class__(
            Seq(function(s.id), s.desc, s.seq)
            for s
            in self.seqs
        )

    def map_desc(self, function):
        """ Apply a function to all sequence descriptions.

        Examples:
        >>> fasta = [
        ...     ">test1 description",
        ...     "ATGCA",
        ...     ">test2 descr",
        ...     "TGACA",
        ... ]
        >>> seqs = Seqs.parse(fasta).map_desc(lambda d: None)
        >>> seqs = iter(seqs)
        >>> next(seqs)
        Seq(id='test1', desc='None', seq='b'ATGCA'')
        >>> next(seqs)
        Seq(id='test2', desc='None', seq='b'TGACA'')
        """

        return self.__class__(
            Seq(s.id, function(s.desc), s.seq)
            for s
            in self.seqs
        )

    def map_seq(self, function):
        """ Apply a function to all sequence sequences.

        Examples:
        >>> fasta = [
        ...     ">test1 description",
        ...     "atgca",
        ...     ">test2 descr",
        ...     "TGACA",
        ... ]
        >>> seqs = Seqs.parse(fasta).map_seq(lambda s: s.lower())
        >>> seqs = iter(seqs)
        >>> next(seqs)
        Seq(id='test1', desc='description', seq='b'atgca'')
        >>> next(seqs)
        Seq(id='test2', desc='descr', seq='b'tgaca'')
        """

        return self.__class__(
            Seq(s.id, s.desc, function(s.seq))
            for s
            in self.seqs
        )

    def min_length(self, length):
        """ Filters out sequences that are shorter than length """

        return self.filter(lambda s: len(s) >= length)

    def max_length(self, length):
        """ Filters out sequences that are longer than length """

        return self.filter(lambda s: len(s) <= length)

    def deduplicated(self, id_conv, column="id"):
        """ Removes duplicates from a Seqs object and stores a mapping file.

        Examples:
        >>> inseqs = [
        ...     Seq('test1', None, "ATGCA"),
        ...     Seq('test2', None, "ATGCA"),
        ...     Seq('test3', None, "ACGTA"),
        ... ]
        >>> seqs = Seqs(inseqs).deduplicated(lambda id: id)
        >>> seq_iter = iter(seqs)
        >>> next(seq_iter)
        Seq(id='test1', desc='None', seq='b'ATGCA'')
        >>> next(seq_iter)
        Seq(id='test3', desc='None', seq='b'ACGTA'')
        >>> seqs.id_map[0][0], seqs.id_map[0][1]
        ('test1', 'test1')
        >>> seqs.id_map[1][0], seqs.id_map[1][1]
        ('test1', 'test2')
        >>> seqs.id_map[2][0], seqs.id_map[2][1]
        ('test3', 'test3')
        """

        return SeqDeduplicated(self.seqs, column=column, id_conv=id_conv)

    def replace_ids(self, id_conv, column="id"):
        """ Replaces all seq ids given a function and stores a mapping file.

        Examples:
        >>> inseqs = [
        ...     Seq('test1', None, "ATGCA"),
        ...     Seq('test2', None, "ATGCA"),
        ...     Seq('test3', None, "ACGTA"),
        ... ]
        >>> seqs = Seqs(inseqs).replace_ids(lambda id: "".join(reversed(id)))
        >>> seq_iter = iter(seqs)
        >>> next(seq_iter)
        Seq(id='1tset', desc='None', seq='b'ATGCA'')
        >>> next(seq_iter)
        Seq(id='2tset', desc='None', seq='b'ATGCA'')
        >>> next(seq_iter)
        Seq(id='3tset', desc='None', seq='b'ACGTA'')
        >>> seqs.id_map[0][0], seqs.id_map[0][1]
        ('1tset', 'test1')
        >>> seqs.id_map[1][0], seqs.id_map[1][1]
        ('2tset', 'test2')
        >>> seqs.id_map[2][0], seqs.id_map[2][1]
        ('3tset', 'test3')
        """
        return SeqReId(self.seqs, column=column, id_conv=id_conv)

    def __iter__(self):
        return iter(self.seqs)

    def __str__(self):
        """ Prints the sequences in fasta format.

        Examples:
        >>> inseqs = [
        ...     Seq('test1', None, "ATGCA"),
        ...     Seq('test2', None, "ATGCA"),
        ...     Seq('test3', None, "ACGTA"),
        ... ]
        >>> print(Seqs(inseqs))
        >test1
        ATGCA
        >test2
        ATGCA
        >test3
        ACGTA
        <BLANKLINE>
        """

        return "".join(map(str, self.seqs))


class SeqDeduplicated(Seqs):

    def __init__(self, seqs, column="id", id_conv=lambda x: x):
        self.id_map = list()
        if column == "desc":
            self.seqs = self._filter_desc(seqs, id_conv)
        else:
            self.seqs = self._filter_id(seqs, id_conv)
        return

    def _filter_id(self, seqs, id_conv=None):
        seen = dict()
        for record in seqs:
            checksum = record.checksum()

            if checksum in seen:
                new_id = seen[checksum]
                self.id_map.append((new_id, record.id, checksum, record.desc))
            else:
                new_id = id_conv(record.id)
                self.id_map.append((new_id, record.id, checksum, record.desc))
                seen[checksum] = new_id

                yield Seq(new_id, record.desc, record.seq)
        return

    def _filter_desc(self, seqs, id_conv=None):
        seen = dict()
        for record in seqs:
            checksum = record.checksum()

            if checksum in seen:
                new_desc = seen[checksum]
                self.id_map.append((
                    new_desc, record.desc, checksum, record.id
                ))
            else:
                new_desc = id_conv(record.desc)
                self.id_map.append((
                    new_desc, record.desc, checksum, record.id
                ))
                seen[checksum] = new_desc

                yield Seq(record.id, record.desc, record.seq)
        return

    def flush_ids(self, handle):
        for new_id, old_id, checksum, old_desc in self.id_map:
            if old_desc is None:
                old_desc = "."
            handle.write(f"{new_id}\t{old_id}\t{checksum}\t{old_desc}\n")

        self.id_map = list()
        return


class SeqReId(Seqs):

    def __init__(self, seqs, column="id", id_conv=lambda x: x):
        self.id_map = list()
        if column == "desc":
            self.seqs = self._filter_desc(seqs, id_conv)
        else:
            self.seqs = self._filter_id(seqs, id_conv)
        return

    def _filter_id(self, seqs, id_conv=None):
        seen = dict()
        for record in seqs:
            if record.id in seen:
                new_id = self.seen[record.id]
            else:
                new_id = id_conv(record.id)
                seen[record.id] = new_id

            self.id_map.append((new_id, record.id, record.desc))

            yield Seq(new_id, record.desc, record.seq)
        return

    def _filter_desc(self, seqs, id_conv=None):
        seen = dict()
        for record in seqs:
            if record.desc in seen:
                new_desc = seen[record.desc]
            else:
                new_desc = id_conv(record.desc)
                seen[record.desc] = new_desc

            self.id_map.append((new_desc, record.desc, record.id))

            yield Seq(record.id, new_desc, record.seq)
        return

    def flush_ids(self, handle):
        for new_id, old_id, old_desc in self.id_map:
            if old_desc is None:
                old_desc = "."
            handle.write(f"{new_id}\t{old_id}\t{old_desc}\n")

        self.id_map = list()
        return
