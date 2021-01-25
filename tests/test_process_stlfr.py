import os
import tempfile
from contextlib import contextmanager

from blr.cli.process_stlfr import BarcodeGenerator


@contextmanager
def tempinput(data):
    # Temporarly generate file object
    # https://stackoverflow.com/questions/11892623/stringio-and-compatibility-with-with-statement-context-manager
    temp = tempfile.NamedTemporaryFile(delete=False)
    temp.write(data)
    temp.close()
    try:
        yield temp.name
    finally:
        os.unlink(temp.name)


def test_generate_barcodes_no_reference():
    index_to_barcode = BarcodeGenerator(barcodes_file=None)
    indexes = ["1_1_1", "1_1_2", "1_1_3", "1_1_1"]
    barcodes = [index_to_barcode.get(i) for i in indexes]
    assert barcodes[0] != barcodes[1] != barcodes[2]
    assert barcodes[0] == barcodes[-1]


def test_generate_barcodes_with_reference():
    with tempinput(b"AA      1\nTT      2\nCC      3") as barcode_file:
        index_to_barcode = BarcodeGenerator(barcodes_file=barcode_file)
    indexes = ["1_1_1", "1_1_2", "1_1_3", "1_1_1"]
    barcodes = [index_to_barcode.get(i) for i in indexes]
    correct_barcodes = ["AAAAAA", "AAAATT", "AAAACC", "AAAAAA"]
    assert all(b == b_ref for b, b_ref in zip(barcodes, correct_barcodes))
