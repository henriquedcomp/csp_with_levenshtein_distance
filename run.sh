#!/bin/bash

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
FILE="/random/binary/BINARY_20_40_0"

cmake --build build

if [ $? -eq 0 ]; then
    path=$(realpath --relative-to="." "instances/$FILE")
    dir=$(dirname "$path")

    output_dir=$(echo "$dir" | sed "s|^[^/]*|results|")

    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
    fi

    for model in "${MODELS[@]}"; do
        build/csp_with_levenshtein_distance "instances/$FILE.in" "results/$FILE.out" "$model"
    done
fi