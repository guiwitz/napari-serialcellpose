
try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"


from ._reader import get_serial_reader


from .serial_widget import SerialWidget
