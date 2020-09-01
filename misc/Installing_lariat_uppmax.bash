#/bin/bash

git clone --recursive https://github.com/10XGenomics/lariat.git

sed -i 's/char \* bwa_pg/\/\/char \* bwa_pg/g' lariat/go/src/gobwa/bwa_bridge.c

module load gcc zlib go

cd lariat/go
make

echo -e '\n\n========\nTest run\n========'
bin/lariat -h
