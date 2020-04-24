#!/bin/bash

# Region to run script on
pathToData=$1
region=$2
company=$3
seq=$4
params=$5
coverage=${6}

# Record commands
set -o xtrace

# Create temp dir for output results
tempDir=${region}_${company}_${seq}_unphaseTest
mkdir ${tempDir}
cd ${tempDir}

# Run margin
echo Running Margin
time ../margin ${pathToData}/haploidTestExamples/${company}/${seq}/${region}/*.${coverage}x.bam ${pathToData}/haploidTestExamples/${company}/${seq}/${region}/*shasta.fasta ../../params/${company}/${seq}/${params} --logLevel DEBUG

# Calculate identity
echo Comparing predicated sequence to true sequence
time python3 ../../scripts/dirty_assembly_compare.py ${pathToData}/haploidTestExamples/${company}/${seq}/${region}/*truth.fasta output.fa

# Move back
cd ..
