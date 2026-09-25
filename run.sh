#!/bin/sh

# write using /instances as base dir and do not add extension
FILE="/random/dna/DNA_5_20_2"

cmake --build build

if [ $? -eq 0 ]; then
    time build/csp_with_levenshtein_distance "instances/$FILE.in" "results/$FILE.out"
fi