import logging
logger = logging.getLogger(__name__)

class MEGssiawError(Exception):
    """
    Base exception for LuessiawMNE.
    """

    def __init__(self, message: str):
        super().__init__(message)
        self.message = message

    def __str__(self):
        return f"[LuessiawMNE] {self.message}"