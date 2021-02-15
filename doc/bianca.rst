Installing on Bianca
====================
Instruction for how to install the full BLR pipe for running on the UPPMAX Bianca cluster.

Pack enviroment on Rackham 
------------------------------------
- Login to rackham

- Install conda-pack (see install instructions on https://conda.github.io/conda-pack/) or just run

.. code-block::
  
  pip install conda-pack
  
- Clone the BLR repo and move inside

.. code-block::

  git clone https://github.com/NBISweden/BLR.git
  cd BLR

- If you want to use a different version than the latest checkout the version using ``git checkout``. You can also check out any tagged release listed from running ``git tag  -l``.
- Create a new conda environment. Note that the file ``environment.linux.lock.yml`` is used!

.. code-block::

  conda env create -n my-env -f environment.linux.lock.yml 

- Activate the created environment


.. code-block::

  conda activate my-env


- Install BLR in non-editable mode using

.. code-block::

  pip install .

- Run conda-pack to create an archive of the environment.

.. code-block::

  conda pack -n my-env  

- This creates a archive file called ``my-env.tar.gz``. Move this file to the bianca cluster, see instructions in https://www.uppmax.uu.se/support/user-guides/bianca-user-guide/. 

Unpack archive on Bianca
------------------------
- Create a new directory and unpack the archive therein.

.. code-block::

  mkdir -p path/to/my-env
  tar -xzf my-env.tar.gz -C path/to/my-env
  
- Activate the environment using the following command

.. code-block::
  
  source path/to/my-env/bin/activate

- Run the following command to unpack the environment.

.. code-block::
  
  conda-unpack
  
- Test that BLR is successfully installed by doing a test run
- To deactivate the enviroment run

.. code-block::
  
  source path/to/my-env/bin/deactivate


