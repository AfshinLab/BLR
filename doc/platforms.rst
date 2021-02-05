Platforms
==========
This text contains information about the different linked-read platforms supported by the pipeline and which requirements for input and possible preprocessing they have.

- BLR_ (fully supported)
- `10x Genomics`_ (supported)
- stLFR_ (supported)
- TELL-seq_ (supported)
- CPT-seq_ (comming soon)
- Haplotagging_ (comming soon)



BLR
---
Barcode Linked Reads (BLR) is based on the technology described in `Redin el al. 2019 <https://doi.org/10.1038/s41598-019-54446-x>`__. Long DNA fragments are subjected to tagmentation using Tn5-covered beads to cut, tag and wrap the fragment around the beads. The DNA-wrapped beads are then used in emmulsion PCR along with barcoded oligo. Within each emmulsion droplet the barcode and tagged fragments are amplified independently and then linked using overlap-extension. Barcode-linked fragments are recovered and indexed for Illumina sequencing.

:Barcodes:
    Semi-degerate sequence of 20 bases with about 3.5 billion possible sequences.
:Preprocessing:
    Barcode sequences are extracted from read 1 and clustered using `starcode <https://github.com/gui11aume/starcode>`__ to error-correct barcodes. Read 1 and 2 is then trimmed using `cutadapt <https://github.com/marcelm/cutadapt>`__ and the corrected barcode attached to the header.
:Reguired input files:
    Read1 FASTQ,
    Read2 FASTQ

10x Genomics
------------
10x Genomics linked-read technology comes in two versions; the older `GemCode (v1) <https://doi.org/10.1038/nbt.3432>`__ and more recent `Chromium Genome (v2) <https://doi.org/10.1101/gr.234443.118>`__. Long DNA fragments are combined in droplets with barcode-containing gel-beads to create GEMs ((Gel Bead-In EMulsions). The fragments are amplified and barcoded using a combination of free random hexamers and barcode-linked random hexamers from the gel beads. Following this barcoded fragments are recovered and fragments before ligation of 3' sequencing adaptor. Libraries are sequenced using Illumina Sequencing. The commercial version of the technology is `currently discontinued <https://www.10xgenomics.com/products/linked-reads>`__.

:Barcodes:
    10x Genomics uses a barcode library of 16 base sequences. GemCode libraries have maximum of 737 thousand sequences while the Chromium Genome has about 4 million.
:Preprocessing:
    Barcode need to be extracted from read 1 and error-corrected using a whitelist (Download from `here <https://github.com/10XGenomics/supernova/tree/master/tenkit/lib/python/tenkit/barcodes>`__). Reads are trimmed and the corrected barcode appended to the header. Preprocessing uses the `ema <https://github.com/arshajii/ema>`__ ``count`` and ``preproc`` tools.
:Required input files:
    Read1 FASTQ,
    Read2 FASTQ,
    Barcode whitelist


stLFR
-----
stLFR (single-tube long fragment read) is based on the technology described in `Wang et al. 2019 <https://doi.org/10.1101/gr.245126.118>`__ and is commercially available from `MGI <https://en.mgi-tech.com/products/reagents_info/18/>`__. The technology uses tagmentation to individually cut-and-hold long DNA fragments in solution. The tagmentase-DNA complex is then hybridized and individual wrapped around barcoded beads through the adaptor introduced by the tagmentation. The barcode is then ligated to each subfragment before recovery and final library prepration. Sequencing is preformed on the DNBSEQ platfroms.

:Barcodes:
    The barcodes are genereated using a combinatorial split-and-pool approach. The barcode is a combination of three barcodes from a 1,536 barcode library in which each barcode is 10 bp. The total number of possible barcodes is about 3.6 billion.
:Preprocessing:
    **The barcode needs to be extracted and error-corrected externaly using the** `stLFR_read_demux <https://github.com/stLFR/stLFR_read_demux>`__ **tool before beign inputed to the pipeline**. Following this the reads are trimmed using `cutadapt <https://github.com/marcelm/cutadapt>`__. The stLFR_read_demux inserts a index based barcode to the read header (e.g. ``#1024_323_231``) based on which three barcodes were detected. This is not compatible with some aligners such as `ema <https://github.com/arshajii/ema>`__ due to not cosisting of IUPAC base symbols. Therefore the index barcodes are replaced with either (A) a concatemer of the three detect barcodes or (B) a uniquely generated 16 base sequence (recommended). For option B a whitelist of barcodes from the stLFR_read_demux (Download from `here <https://github.com/stLFR/stLFR_read_demux/blob/master/scripts/barcode.list>`__ tools needs to be provided.
:Required input files:
    Read1 FASTQ,
    Read2 FASTQ,
    Barcode whitelist (Optional)


TELL-seq
--------
TELL-seq is based on the technology from `Chen et al. 2020 <https://doi.org/10.1101/gr.260380.119>`_ and is commercially available from the company `Universal Sequencing <https://www.universalsequencing.com/>`__. The method uses clonaly barcode beads with attacted tagmentases to cut and barcode individual long DNA fragments in solution. A second tagmentation is also preformed in solution to introduce a second adaptor. The library is sequenced using Illumina sequencing with special setup to sequence the barcode as index 1.

:Barcodes:
    Semi-degenerated sequence of 18 bases with about 2.4 billion possible sequences.
:Preprocessing:
    Barcodes are either (A) clustered using `starcode <https://github.com/gui11aume/starcode>`__ as for BLR_ or (B) single count barcodes are corrected to any barcode within a hamming distance of 1 or discarded. Option B follows the method used in `Chen et al. 2020`_. Reads are subsequently tagged in the header with the corrected barcode.
:Required input files:
    Read1 FASTQ,
    Read2 FASTQ,
    Index1 FASTQ (contains barcodes)


CPT-seq
-------
Technology based `Amini et al. 2014 <https://doi.org/10.1038/ng.3119>`_ and the follow-up CPTv2-seq from `Zhang et al. 2017 <https://doi.org/10.1038/nbt.3897>`_. These technologies were developed by Illumin but are not commercially available.

Haplotagging
------------
Haplotagging is based on the technology presented in `Meier et al. 2020 <https://doi.org/10.1101/2020.05.25.113688>`_. The technology uses barcoded beads covered with Tn5 tagmentase to cut and barcode individual long DNA fragments in solution. The beads are coated in a combination of two barcodes AB and CB that become inserted at the 5' and 3' of each cut fragment. Barcodes are combinatorialy generated with about 85 million possible combinations in total.
