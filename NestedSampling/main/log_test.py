# logging_example.py

import logging

# Create a custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Create handlers
f_handler = logging.FileHandler('file.log',mode='w')
# f_handler.setLevel(logging.DEBUG)

# Create formatters and add it to handlers
f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
f_handler.setFormatter(f_format)

# Add handlers to the logger
logger.addHandler(f_handler)

logger.warning('This is a warning')
logger.error('This is an error')
logger.info('this is info rofl')
