Installing Lariat
=================

To use lariat_ for alignment you need to manually install it within your
environment. See instructions below for how to install it on MacOS and on Uppmax.

.. _lariat: https://github.com/10XGenomics/lariat

MacOS install instructions
--------------------------

First create a new environment with whish to build lariat from source

.. code-block::

 conda create -n lariat-build
 conda activate lariat-build
 conda install go
 condat install clangxx_osx-64

Clone and build lariat.

.. code-block::

 git clone https://github.com/10XGenomics/lariat
 cd lariat
 git submodule update --init --recursive
 cd go
 make

Test install by running:

.. code-block::

 bin/lariat -h


Uppmax install instructions
---------------------------

Run the following within a the directory you which to install lariat.
The is also a bash script ``/misc/lariat_uppmax_install.sh`` for this purpose.

First clone the lariat repository

.. code-block::

 git clone --recursive  https://github.com/10XGenomics/lariat.git


Edit the file 'bwa_bridge.c' to comment out the final line. The final line should be:
 ``//char * bwa_pg = "10X Genomics";``.

.. code-block::

 vi lariat/go/src/gobwa/bwa_bridge.c

Loading the required packages for building lariat:

.. code-block::

 module load gcc zlib go


Build lariat

.. code-block::
  
  cd go 
  make 
  bin/lariat -h 

Add lariat to environment
-------------------------
For the pipeline to recognize the lariat install it must be available within
the blr environment. To add lariat to your blr environment, symlink it to the ``bin``
folder of your env.

.. code-block::

    ln -s /path/to/lariat/go/bin/lariat /path/to/miniconda/envs/my-blr-env/bin/.

