#!/bin/bash

cmake --build build

if [ $? -eq 0 ]; then
    find "./instances" -type f -name "*.in" -print0 | while IFS= read -r -d '' file; do
        path=$(realpath --relative-to="." "$file")
        dir=$(dirname "$path")

        output_dir=$(echo "$dir" | sed "s|^[^/]*|results|")

        if [ ! -d "$output_dir" ]; then
            mkdir -p "$output_dir"
        fi

        filename=$(basename "$file" ".in")

        time build/csp_with_levenshtein_distance "$path" "$output_dir/$filename.out"
    done
fi