from datetime import datetime
from projection.logging_config import logger

def printDatetime(message = "", date = datetime.now()):
    message += str(date)
    print(message)
    logger.info(message)

