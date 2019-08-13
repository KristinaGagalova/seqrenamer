class SRException(Exception):
    """ Subclasses of this should all have the msg attribute.

    Re-assign the ecode to whatever is appropriate.
    """

    ecode = 1

    def __init__(self, msg):
        self.msg = msg


class InvalidArgumentError(SRException):
    ecode = 1


class MapFileParseError(SRException):
    ecode = 2


class MapFileKeyError(SRException):
    ecode = 3


class XsvColumnNumberError(SRException):
    ecode = 4
