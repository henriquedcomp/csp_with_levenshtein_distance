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

RESET='\e[0m'
BLUE='\e[34m'

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

        echo -e "${BLUE}[INFO] Processing instance $path${RESET}"

        for model in "${MODELS[@]}"; do
            build/csp_with_levenshtein_distance "$path" "$output_dir/$filename.out" "$model"

            exit_code=$?

            if (( exit_code != 0 )); then
                {
                    if (( exit_code == 143 )); then
                        echo "FAILURE: Terminated by SIGTERM (possible memory exhaustion)"
                    elif (( exit_code == 137 )); then
                        echo "FAILURE: Killed by SIGKILL (possible OOM)"
                    else
                        echo "FAILURE: Exit code $exit_code"
                    fi

                    echo ""
                } >> "$output_dir/$filename.out"
            fi
        done
    done
fi