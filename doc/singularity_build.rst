Instructions to build a singularity file for BLR
=================================================
This is practical in many scenarios to avoid installations steps for BLR and all the dependencies (in offline systems for instance).
Thus, to be able to run BLR directly from a single file.


Steps used to build the file:
-----------------------------

- Make sure you have singualrity installed in your system, and you have a root access in it.
- Make a new directory.

- Clone Blr to it.

- Add the following text in a "build.txt" file:

.. code-block::

      Bootstrap: docker

      From: continuumio/miniconda3

      %files
          BLR/ /opt
      %environment
        export LC_ALL=C
        export PATH=/opt/conda/envs/blr/bin:$PATH


      %post
          /opt/conda/bin/conda env create --name blr -f /opt/BLR/environment.yml
          . /opt/conda/etc/profile.d/conda.sh
          conda activate blr
          python3 -m pip install --no-cache-dir /opt/BLR/
          conda clean --all --yes

      %runscript
          exec /opt/conda/envs/blr/bin/"$@"


- Run the command:

.. code-block::

     sudo singularity build blr_image build.txt

- Test BLR as follows:

.. code-block::

     ./blr_image blr -h
