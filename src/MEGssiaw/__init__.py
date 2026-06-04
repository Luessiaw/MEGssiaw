__version__ = 26.0508
__author__ = "Luessiaw"

from .mainSolver import *
from .version import get_git_commit
from datetime import datetime
import logging

meta = {
    "time": str(datetime.now()),
    "LuessiawMNE_commit": get_git_commit(),
}

logger = logging.getLogger("LuessiawMNE")
logger.addHandler(logging.NullHandler())

logging.basicConfig(
    format='%(asctime)s [%(levelname)s] %(message)s', 
    level=logging.INFO,
    datefmt='%H:%M:%S'
)

