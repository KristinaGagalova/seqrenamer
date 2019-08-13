class SRException(Exception):
    """ Subclasses of this should all have the msg attribute.

    Re-assign the ecode to whatever is appropriate.
    """

    ecode = 1

    def __init__(self, msg):
        self.msg = msg


class InvalidArgumentException(SRException):
    ecode = 1


class MapFileParseError(SRException):
    ecode = 2


class MapFileKeyError(SRException):
    ecode = 3
