#!/bin/bash

# Region to run script on
region=$1
options=$2

# Record commands
set -o xtrace

# Create temp dir for output results
tempDir=${region}_phaseTest
mkdir ${tempDir}
cd ${tempDir}

# Run margin
echo Running Margin
time ../margin ../../tests/data/diploidTestExamples/${region}/*.bam ../../tests/data/diploidTestExamples/${region}/HG002.shasta.*.fasta ../../params/allParams.np.json ${options} --logLevel DEBUG --diploid

# Calculate identity
echo Comparing predicated haplotype1 to true haplotype1
time python3 ../../scripts/dirty_assembly_compare.py ../../tests/data/diploidTestExamples/${region}/HG002_h1*.fa output.fa.hap1

echo Comparing predicated haplotype2 to true haplotype1
time python3 ../../scripts/dirty_assembly_compare.py ../../tests/data/diploidTestExamples/${region}/HG002_h1*.fa output.fa.hap2

echo Comparing predicated haplotype1 to true haplotype2
time python3 ../../scripts/dirty_assembly_compare.py ../../tests/data/diploidTestExamples/${region}/HG002_h2*.fa output.fa.hap1

echo Comparing predicated haplotype2 to true haplotype2
time python3 ../../scripts/dirty_assembly_compare.py ../../tests/data/diploidTestExamples/${region}/HG002_h2*.fa output.fa.hap2

# Build collection of differences
time python3 ../../scripts/dirty_assembly_compare.py output.fa.hap1 output.fa.hap2 verbose >predictedMismatches.txt
time python3 ../../scripts/dirty_assembly_compare.py ../../tests/data/diploidTestExamples/${region}/HG002_h1*.fa ../../tests/data/diploidTestExamples/${region}/HG002_h2*.fa verbose >trueMismatches.txt

# Compare differences
echo Comparing predicted hets
time python3 ../../scripts/compareHets.py trueMismatches.txt predictedMismatches.txt > validatedMismatches.txt
cat validatedMismatches.txt

# Make POA visualizations of each
echo Making POA graph visualizations of predicted hets
mkdir visualizations
for i in `cat validatedMismatches.txt| grep 'Predicted mismatch:' | cut -f5 -d' '`
do
  j=`echo ${i} - 5 | bc`
  python3 ../../scripts/phasedPoaToDot.py output_poa.csv.hap1 --start=${j} --length=10 --output=visualizations/hap1_poa_${j}
done

# Exit
cd ..