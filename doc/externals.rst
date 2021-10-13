External Tools
==============
This document contains information about external software and fixes done in order to run it with the BLR pipeline.

- NAIBR_
- EMA_


NAIBR
-----
`NAIBR <https://github.com/raphael-group/NAIBR>`__ is used to call structural variants.

- `mpmath` and `future` is not in original git environment description but needed to run.
- Needs to be run with work directory = NAIBR folder (Otherwise won't find its own modules).
- Input BAM file needs to be indexed.

Citation
  Rebecca Elyanow, Hsin-Ta Wu, Benjamin J Raphael, Identifying structural variants using linked-read sequencing data, Bioinformatics, Volume 34, Issue 2, 15 January 2018, Pages 353–360, https://doi.org/10.1093/bioinformatics/btx712


EMA
---
`EMA <https://github.com/arshajii/ema>`__ or Emerald is used for preprocessing of 10x reads and alignment for all platforms. 

Citation
  Ariya Shajii, Ibrahim Numanagić, Christopher Whelan, Bonnie Berger, Statistical Binning for Barcoded Reads Improves Downstream Analyses, Cell Systems, Volume 7, Issue 2, 2018, p219-226.e5, ISSN 2405-4712, https://doi.org/10.1016/j.cels.2018.07.005.
