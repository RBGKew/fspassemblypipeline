#!/bin/bash

# Directory containing the fastq.gz files
INPUT_DIR="InputDir"

echo "Starting job on $HOSTNAME"

# Create sample list
cd $INPUT_DIR
touch sample.list
for file in *.gz; do
    first_element="${file%%_R*}"
    if ! grep -q "^$first_element$" sample.list; then
        echo "$first_element" >> sample.list
		mkdir $first_element
		mv $first_element\_*gz $first_element
    fi
done

echo "Job finished"
