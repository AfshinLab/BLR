import os
import tempfile
from contextlib import contextmanager


@contextmanager
def tempinput(string: bytes):
    # Temporarly generate file object
    # https://stackoverflow.com/questions/11892623/stringio-and-compatibility-with-with-statement-context-manager
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.write(string)
    temp.close()
    try:
        yield temp.name
    finally:
        os.unlink(temp.name)
