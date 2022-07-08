
import os
import os.path
from datetime import datetime
import logging



def setup_process(config):
    """Logging and I/O setup for the current processes."""

    def ensure_directory(fname):
        """Make sure the directory containing the given name exists."""

        dirname = os.path.dirname(fname)
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    fname_substitutions = dict(
        now=datetime.now().strftime(config.fname_datetime_format),
        system=config.system,
        pid=os.getpid()
    )

    std_out_err_fname = config.std_out_err_fname % fname_substitutions
    ensure_directory(std_out_err_fname)

    io_destination = os.open(
        std_out_err_fname,
        os.O_WRONLY | os.O_TRUNC | os.O_CREAT | os.O_DSYNC,
        mode=0o666
    )
    os.dup2(io_destination, 1)
    os.dup2(io_destination, 2)

    logging_fname = config.logging_fname % fname_substitutions
    ensure_directory(logging_fname)

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        handler.close()
    logging_config = dict(
        filename=logging_fname,
        level=logging.DEBUG,#getattr(logging, config.logging_verbosity.upper()),
        format=config.logging_message_format,
    )
    if config.logging_datetime_format is not None:
        logging_config['datefmt'] = config.logging_datetime_format
    logging.basicConfig(**logging_config)
