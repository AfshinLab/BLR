from blr.cli.process_stlfr import BarcodeGenerator


def test_generate_barcodes_no_reference():
    index_to_barcode = BarcodeGenerator()
    indexes = ["1_1_1", "1_1_2", "1_1_3", "1_1_1"]
    barcodes = [index_to_barcode.get(i) for i in indexes]
    assert barcodes[0] != barcodes[1] != barcodes[2]
    assert barcodes[0] == barcodes[-1]
