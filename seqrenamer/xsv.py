import csv
from seqrenamer.exceptions import XsvColumnNumberError


class Xsv(object):

    def __init__(self, handle, comment="#", sep=","):

        self.handle = self._filter_comments(handle, comment)
        self.id_map = list()
        self.sep = sep
        return

    def __iter__(self):
        csv_reader = csv.reader(
            self.handle,
            delimiter=self.sep,
            dialect='excel'
        )

        for row in csv_reader:
            yield row

        return

    def replace_ids(self, id_conv, column, header):
        csv_reader = iter(self)
        seen = dict()

        for row in csv_reader:
            if header:
                header = False
                yield row
                continue

            try:
                old_id = row[column]

                if old_id in seen:
                    new_id = seen[old_id]
                else:
                    new_id = id_conv(old_id)
                    seen[old_id] = new_id

                row[column] = new_id

            except IndexError:
                joined_line = self.sep.join(map(str, row))
                raise XsvColumnNumberError(
                    f"Could not access column '{column}' in a line. "
                    f"The offending line was: {joined_line}."
                )

            self.id_map.append((new_id, old_id))
            yield row

        return

    def flush_ids(self, handle):
        for new_id, old_id in self.id_map:
            handle.write(f"{new_id}\t{old_id}\n")

        self.id_map = list()
        return

    @staticmethod
    def _filter_comments(handle, comment):
        for line in handle:
            if not line.startswith(comment):
                yield line
        return
