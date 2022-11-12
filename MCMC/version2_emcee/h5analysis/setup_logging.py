import logging
import os

def setup_logging(parsed_args):

    pid = os.getpid()

    logging_fname = f'{pid}.log'

    logging.basicConfig(filename=logging_fname,
                        level=logging.DEBUG)
