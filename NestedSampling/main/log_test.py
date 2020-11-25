# logging_example.py

import logging

def add_handler(logger):

    # Create handlers
    f_handler = logging.FileHandler('file.log',mode='a')

    # Create formatters and add it to handlers
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(f_handler)

    logger.warning('This is a warning')
    logger.error('This is an error')
    logger.info('this is info rofl')


def remove_handler(logger):

    logger.removeHandler(logging.getLogger().handlers)

# Create a custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# f_handler.setLevel(logging.DEBUG)

add_handler(logger)
remove_handler(logger)
add_handler(logger)

