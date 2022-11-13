import logging
import os


def setup_logging(parse_args):

    pid = os.getpid()

    std_out_err_fname = parse_args.std_out_err_path + f'/{pid}.err'

    io_destination = os.open(
        std_out_err_fname,
        os.O_WRONLY | os.O_TRUNC | os.O_CREAT | os.O_DSYNC,
        mode=0o666
    )
    os.dup2(io_destination, 1)
    os.dup2(io_destination, 2)


    logging_fname = parse_args.logging_path + f'/{pid}.log'

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        handler.close()

    logging.basicConfig(filename=logging_fname,
                        level=logging.DEBUG)
