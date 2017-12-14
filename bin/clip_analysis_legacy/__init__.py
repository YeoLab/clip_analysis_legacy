import os

def data_dir():
    """
    Returns the data directory that contains files for data and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'data')