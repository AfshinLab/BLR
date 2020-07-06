# MultiQC BLR testdata

Testdata is divided into two folders. 

1. `data`: This contains the files that tests are runned on and are the 
files that the respective modules are expected to handle.
2. `reference`: This contains the expected output generated from running 
MultiQC on the files in data. The file are generated in the `multiqc_data` 
folder that is normaly output when running MultiQC. These files are compared 
in the test line-by-line to the output from running on the files in `data`.      