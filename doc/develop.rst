Development
===========

- `Installing BLR`_
- Testing_
- Linting_
- `Conda environment files`_
- Profiling_
- `Snakemake Pipeline`_

  + `Chunk handlig`_


Installing BLR
--------------

To setup BLR for development it is benefitial to install it in editable mode. Use the command below:

.. code-block::

    pip install -e .[dev]

This also install development tools such as ``pytest`` and ``flake8``.


Testing
-------

Testdata needs to be downloaded and setup for running the tests. Copy and run the block the block below.

.. code-block::

    testdata_version=0.6
    wget -nv https://export.uppmax.uu.se/uppstore2018173/blr-testdata-${testdata_version}.tar.gz
    tar xf blr-testdata-${testdata_version}.tar.gz
    ln -s blr-testdata-${testdata_version} blr-testdata

Now that everything is setup tests can be run using the `tests/run.sh` script.

.. code-block::

    bash tests/run.sh


Linting
-------

``flake8`` is used for linting. Either run the command below before commiting or setup pre-commit to run this automatically.

.. code-block::

    flake8 src tests


Conda environment files
-----------------------

There are two types of files in this repository that describe conda environments.

- The file ``environment.yml`` contains abstract dependencies such as ``pysam`` or
  ``bowtie2``. This file is managed manually and needs to be
  updated whenever there are new dependencies or when the required version for a
  dependency changes.

- The ``environment.linux-64.lock`` and ``environment.osx-64.lock`` files
  (lock files) contain explicitly-defined environments which are reproducible and platform
  dependant.  These files are used to create the test environments.

Whenever the ``environment.yml`` file is updated, you need to run:

.. code-block::

    conda-lock -f environment.yml -p linux-64 -p osx-64 --filename-template "environment.{platform}.lock"`

to generate the ``environment.{platform}.lock`` files.

Install conda-lock using pip

.. code-block::

    pip install conda-lock

or conda

.. code-block::

    conda install conda-lock -c conda-forge


SAM Tags
--------

Specifications on SAM-tags used for holding information during data processing and which argparse
option flags to use when specifying them in python scripts. The `10x Genomics barcoded BAM format
<https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bam>`_ is followed
where that information is defined.

======= ================== ==============================================================
SAM-tag Option flag        Description
======= ================== ==============================================================
BX      -b, --barcode-tag  String for the error-corrected barcode                        
MI      -m, --molecule-tag Integer index for an identified molecule for each barcode [*]_ 
RX      -s, --sequence-tag String for the uncorrected barcode sequence                   
HP                         Integer (1 or 2) for the read haplotype assigned
PS                         Integer for the phase set (phaseblock) that the read is part of
PC                         Integer for the quality of the phase set (phaseblock)
======= ================== ==============================================================

.. [*] Note that the index is only unique within the particular chunk.


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


Snakemake Pipeline
------------------

Chunk handlig
^^^^^^^^^^^^^^
Chunks are the separate portions of the mapped bam that go through postprocessing. Each chunk might contain one or 
more contigs. Chunks are handled in the Snakemake workflow through a dictionary called ``chunks`` which define contigs 
as different sets. Each set contains a list of chunks which inturn contain lists of contigs composing each chunk. Chunks are referred to by 
the first contig name. The three primary sets are accessed by the following keys:

  ``'all'`` = handles every contig in reference

  ``'primary'`` = handle every contig in reference that should go through certain post-processing steps (see below). Is a subset of 'all'.

  ``'phased'`` = handles every contig in reference that is diploid i.e. can be phased. Is a subset of 'primary'.

Several subsets of these are also defined for convinence.

  ``'not_phased' = 'all' - 'phased'``

  ``'not_primary' = 'all' - 'primary'``

  ``'primary_not_phased' = 'primary' - 'phased'`` 

These sets are used to control which contigs go through which processing steps. Which contigs are included are defined 
through the ``phasing_contigs`` (for ``'phased'``) and ``contigs_skipped`` (for ``'not_primary'``) parameters in the 
config file ``blr.yml``. 

Processing steps run by ``'primary'`` contigs but not ``'all'``:

- find_clusterdups
- get_barcode_merges
- concat_molecule_stats
- get_barcodes_to_filter
- call_variants
- lsv_calling

Processing steps run by ``'phased'`` contigs but not ``'primary'``:

- hapcut2_extracthairs
- hapcut2_linkfragments
- hapcut2_phasing
- build_config
