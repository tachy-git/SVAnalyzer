#!/bin/bash

for dir in */ ; do
    [ -d "$dir" ] || continue

    shopt -s nullglob
    files=("$dir"/histo_*.root)
    shopt -u nullglob

    if [ ${#files[@]} -gt 0 ]; then
        echo "Merging in $dir"
        ( cd "$dir" && hadd -fk histo.root histo_*.root )
    else
        echo "Skipping $dir (no histo_*.root)"
    fi
done

shopt -s nullglob
all_histo_roots=(*/histo.root)
shopt -u nullglob

if [ ${#all_histo_roots[@]} -gt 0 ]; then
    echo "Merging all subdir histo.root into ./histo.root"
    hadd -fk histo.root "${all_histo_roots[@]}"
else
    echo "No subdir histo.root found"
fi
