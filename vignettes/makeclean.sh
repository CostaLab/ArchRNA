#!/bin/bash
lst=($(ls  -I ifnb_integration.ipynb -I pbmc_h5.ipynb -I brain_10xmatrix.ipynb -I .gitignore -I makeclean.sh))


for files in "${lst[@]}"; do
    echo 'remove '  $files
    rm -r "$files"
done

