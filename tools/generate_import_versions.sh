#!/bin/bash

set -ex

mv import_versions.tsv{,.bak}

(
    echo -e "build\ttable_group\tversion"
    for i in $(find GRCh37 GRCh38 -maxdepth 2 -mindepth 2 -type d)
    do
        echo $i | sed 's/\//\t/g'
    done
) > import_versions.tsv

