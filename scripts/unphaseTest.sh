#!/bin/bash

# Region to run script on
region=$1
options=$2

# Record commands
set -o xtrace

# Create temp dir for output results
tempDir=${1}_unphaseTest
mkdir ${tempDir}
cd ${tempDir}

# Run margin
echo Running Margin
time ../margin ../../tests/data/diploidTestExamples/${region}/*.bam ../../tests/data/diploidTestExamples/${region}/HG002.shasta.*.fasta ../../params/allParams.np.json ${options} --logLevel DEBUG

# Calculate identity
echo Comparing predicated sequence to true haplotype1
time python3 ../../scripts/dirty_assembly_compare.py ../../tests/data/diploidTestExamples/${region}/HG002_h1*.fa output.fa

echo Comparing predicated sequence to true haplotype2
time python3 ../../scripts/dirty_assembly_compare.py ../../tests/data/diploidTestExamples/${region}/HG002_h2*.fa output.fa

# Move back
cd ..
