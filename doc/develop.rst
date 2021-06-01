Development
===========

Conda environment files
-----------------------

There are two types of files that describe Conda environments.

- The file ``environment.yml`` contains abstract dependencies such as ``pysam`` or
  ``bowtie2``. This file is managed manually and needs to be
  updated whenever there are new dependencies or when the required version for a
  dependency changes.

- The ``environment.linux.lock.yml`` and ``environment.osx.lock.yml`` files
  (lock files) contain a fully specified description of the entire environment,
  with locked-down versions.  These files are used to create the test
  environment.

Use the script ``misc/condalock.sh`` to update the lock files whenever you make
changes to ``environment.yml``.


SAM Tags
--------
Specifications on SAM-tags used for holding information during data processing and which argparse
option flags to use when specifying them in python scripts. The `10x Genomics barcoded BAM format
<https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam>`_ is followed
where that information is defined.

..  csv-table::
    :header: "SAM-tag", "Option flag", "Description"
    :widths: 20, 20, 40

    "BX", "-b, --barcode-tag", "SAM tag for storing the error corrected barcode."
    "MI", "-m, --molecule-tag", "SAM tag for storing molecule index specifying a identified molecule
    for each barcode. Note that the index is only unique within the particular chunk."
    "RX", "-s, --sequence-tag", "SAM tag for storing the uncorrected barcode sequence."

Profiling
---------

To run profiling on a particular subcommand you can use the ``--profile`` argument. For example with 
the subcommand ``tagbam`` the command is:

..  code-block:: bash

    blr --profile tagbam input.bam -o output.bam

This command will generate a file called ``blr_tagbam.prof`` with all the profiling information. This 
can then be used with Python's standard library module 
`pstat <https://docs.python.org/3/library/profile.html#pstats.Stats>`_ 
or for example `Snakeviz <https://jiffyclub.github.io/snakeviz/>`_ which allows interaction through the browser. 
