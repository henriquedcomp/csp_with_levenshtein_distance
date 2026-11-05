#!/bin/sh

MODELS=(
    "1"
    "2"
    "3"
    "4"
    "5"
    "6"
    "7"
    "8"
    "9"
    "10"
    "11"
    "12"
    "13"
    "14"
    "15"
    "16"
    "17"
    "18"
    "19"
    "20"
)

# write using /instances as base dir and do not add extension
FILE="/random/protein/PROTEIN_5_40_0"

cmake --build build

if [ $? -eq 0 ]; then
    for model in "${MODELS[@]}"; do
        build/csp_with_levenshtein_distance "instances/$FILE.in" "results/$FILE.out" "$model"
    done
fi