#!/bin/bash

# This script installs the lariat aligner from source in the current folder on Uppmax.
# For more help see: doc/lariat_install.rst

set -euo pipefail

git clone --recursive https://github.com/10XGenomics/lariat.git

sed -i 's/char \* bwa_pg/\/\/char \* bwa_pg/g' lariat/go/src/gobwa/bwa_bridge.c

module load gcc zlib go

cd lariat/go
make

echo -e '\n\n========\nTest run\n========'
bin/lariat -h

echo -e 'Lariat installed successfully!'
echo -e 'Symlink the current installation to your blr enviroment to make it available for running the pipeline'
