import logging
from datetime import datetime
import os

# Create a custom logger
logger = logging.getLogger(__name__)

# Set the default logging level
logger.setLevel(logging.DEBUG)

# Create handlers
c_handler = logging.StreamHandler()

# Generate a timestamp for the log file name
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
log_file_name = f'file_{timestamp}.log'

f_handler = logging.FileHandler(os.path.join("logs",log_file_name))

# Set logging levels for handlers
c_handler.setLevel(logging.WARNING)
f_handler.setLevel(logging.DEBUG)

# Create formatters and add them to handlers
c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
c_handler.setFormatter(c_format)
f_handler.setFormatter(f_format)

# Add handlers to the logger
logger.addHandler(c_handler)
logger.addHandler(f_handler)
