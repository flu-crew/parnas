# -*- coding: utf-8 -*-
import logging
import sys


class LogFormatter(logging.Formatter):
    """
    Adapted from https://stackoverflow.com/a/62488520.
    """
    def format(self, record):
        if record.levelno == logging.INFO:
            self._style._fmt = "   [parnas] %(message)s"
        else:
            color = {
                logging.WARNING: 33,
                logging.ERROR: 31,
                logging.FATAL: 31,
                logging.DEBUG: 36
            }.get(record.levelno, 0)
            self._style._fmt = f"\033[{color}m%(levelname)s\033[0m: %(message)s"
        return super().format(record)


parnas_logger = logging.getLogger('Parnas logger')
formatter = LogFormatter()
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(formatter)

parnas_logger.addHandler(handler)
parnas_logger.setLevel(logging.DEBUG)
