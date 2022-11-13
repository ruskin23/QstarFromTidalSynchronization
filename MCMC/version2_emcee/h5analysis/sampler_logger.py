import logging
import os

def setup_logging():

    pid = os.getpid()

    std_out_err_fname = f'{pid}.outerr'

    io_destination = os.open(
        std_out_err_fname,
        os.O_WRONLY | os.O_TRUNC | os.O_CREAT | os.O_DSYNC,
        mode=0o666
    )
    os.dup2(io_destination, 1)
    os.dup2(io_destination, 2)


    logging_fname = f'{pid}.log'

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        handler.close()

    logging.basicConfig(filename=logging_fname,
                        level=logging.DEBUG)
