#! /bin/bash

v=$1
sed -i -e "s/version:.*$/version: $v/g" wlogdate/meta.yaml

conda-build wlogdate/
conda convert --platform all //anaconda3/conda-bld/osx-64/wlogdate-$v*.tar.bz2 -o //anaconda3/conda-bld/
anaconda upload --force //anaconda3/conda-bld/*/wlogdate-$v*.tar.bz2
