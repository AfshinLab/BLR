To get Lariat correctly installed:
----------------------------------
- git clone --recursive  https://github.com/10XGenomics/lariat.git


Edit the file 'bwa_bridge.c':
---------------------------

- vi lariat/go/src/gobwa/bwa_bridge.c

comment the last line: 

- '// char * bwa_pg = "10X Genomics";' 


Loading required packages on Uppmax:
-------------------------
- module load gcc zlib go


Make it!
--------
.. code-block::
  
  cd go 
  make 
  bin/lariat -h 
