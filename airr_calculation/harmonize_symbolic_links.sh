#!/bin/bash

set -euo pipefail

# This is where ALL symbolic links will be created
DEST="/projects/rpci/songyao/pnfioric/arc_project/airr_calculation/fastq_files"

# Create the destination directory if it does not already exist
mkdir -p "${DEST}"


# Array containing all source directories
# Each of these directories contains subfolders we want linked
SOURCE_DIRS=(
"/projects/rpci/songliu/liyan/LiYan/RQ-026109_1to4/RQ026109-Yao"
"/projects/rpci/songliu/liyan/LiYan/RQ-026109_5to8/RQ026109_Yao/"
"/projects/rpci/songliu/liyan/LiYan/RQ-026110-Yao"
"/projects/rpci/songliu/liyan/LiYan/RQ-026110_13to16/RQ026110-Yao"
"/projects/rpci/songliu/liyan/LiYan/RQ026111/RQ026111-Yao"
)

for SRC in "${SOURCE_DIRS[@]}"; do
    echo "Processing: ${SRC}"

    for SUBDIR in "${SRC}"/*; do

        # Skip if not a directory
        [ -d "${SUBDIR}" ] || continue

        BASENAME=$(basename "${SUBDIR}")
        LINK_PATH="${DEST}/${BASENAME}"

        if [ -e "${LINK_PATH}" ]; then
            echo "Skipping existing: ${LINK_PATH}"
        else
            ln -s "${SUBDIR}" "${LINK_PATH}"
            echo "Linked:"
            echo "  ${LINK_PATH} -> ${SUBDIR}"
        fi
    done
done

echo "Done."