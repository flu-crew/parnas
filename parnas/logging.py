# -*- coding: utf-8 -*-
import logging
import sys


class ParnasLogger(logging.Logger):
    PLAIN = 55  # Add a new logging level.

    def plain(self, message, *args, **kws):
        self.log(self.PLAIN, message, *args, **kws)
logging.addLevelName(ParnasLogger.PLAIN, 'PLAIN')


class LogFormatter(logging.Formatter):
    """
    Adapted from https://stackoverflow.com/a/62488520.
    """

    def __init__(self):
        super().__init__(datefmt="%H:%M:%S")

    def format(self, record):
        if record.levelno == ParnasLogger.PLAIN:
            self._style._fmt = "%(message)s"
        else:
            if record.levelno == logging.INFO:
                self._style._fmt = "[parnas %(asctime)s] %(message)s"
            else:
                color = {
                    logging.WARNING: 33,
                    logging.ERROR: 31,
                    logging.FATAL: 31,
                    logging.DEBUG: 36
                }.get(record.levelno, 0)
                self._style._fmt = f"[\033[{color}m%(levelname)s\033[0m] %(message)s"
        return super().format(record)


logging.setLoggerClass(ParnasLogger)
parnas_logger = logging.getLogger('Parnas logger')
formatter = LogFormatter()
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(formatter)

parnas_logger.addHandler(handler)
parnas_logger.setLevel(logging.INFO)
