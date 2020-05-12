Plattfroms
==========
This text contains information about the different linked-read plattforms supported by the pipeline and which requirements for input and possible preprocessing they have.

- BLR_ (supported)
- `10x Genomics Chromium`_ (supported)
- stLFR_ (comming soon)
- CPT-seq_ (comming soon)
- TELL-seq_ (comming soon)


.. _BLR::
.. _`10x Genomics Chromium`::
.. _stLFR::
.. _CPT-seq::
.. _TELL-seq::


BLR
---
Barcode Linked Reads (BLR) is based on the technology described in `Redin el al. 2019 <https://doi.org/10.1038/s41598-019-54446-x>`_. 


**Input files**

- Read1 FASTQ
- Read2 FASTQ


10x Genomics Chromium
---------------------
**Input files**

- Read1 FASTQ
- Read2 FASTQ
- Barcode whitelist


stLFR
-----
**Input files**

- Read1 FASTQ
- Read2 FASTQ
- Barcode whitelist

**External preprocessing**

The barcode demultiplexing has to be performed unsing the `stLFR_read_demux <https://github.com/stLFR/stLFR_read_demux>`_ tool.  


CPT-seq
-------
**Input files**

- Read1 FASTQ
- Read2 FASTQ


TELL-seq
--------
**Input files**

- Read1 FASTQ
- Read2 FASTQ
- Barcode whitelist.
