#!/bin/bash

# Region to run script on
pathToData=$1
region=$2
company=$3
seq=$4
params=$5
makePlots=$6

# Record commands
set -o xtrace

# Create temp dir for output results
tempDir=${region}_${company}_${seq}_phaseTest
mkdir ${tempDir}
cd ${tempDir}

# Run margin
echo Running Margin
time ../margin ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/*.bam ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/*shasta*.fasta ../../params/${company}/${seq}/${params} --logLevel DEBUG --diploid

# Calculate identity
echo Comparing predicated haplotype1 to true haplotype1
time python3 ../../scripts/dirty_assembly_compare.py ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/HG002_h1*.fa output.fa.hap1 verbose >th1_ph1_mismatches.txt

echo Comparing predicated haplotype2 to true haplotype1
time python3 ../../scripts/dirty_assembly_compare.py ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/HG002_h1*.fa output.fa.hap2 verbose >th1_ph2_mismatches.txt

echo Comparing predicated haplotype1 to true haplotype2
time python3 ../../scripts/dirty_assembly_compare.py ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/HG002_h2*.fa output.fa.hap1 verbose >th2_ph1_mismatches.txt

echo Comparing predicated haplotype2 to true haplotype2
time python3 ../../scripts/dirty_assembly_compare.py ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/HG002_h2*.fa output.fa.hap2 verbose >th2_ph2_mismatches.txt

# Build collection of differences
time python3 ../../scripts/dirty_assembly_compare.py output.fa.hap1 output.fa.hap2 verbose >predictedMismatches.txt
time python3 ../../scripts/dirty_assembly_compare.py ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/HG002_h1*.fa ${pathToData}/diploidTestExamples/${company}/${seq}/${region}/HG002_h2*.fa verbose >trueMismatches.txt

# Compare differences
echo Comparing predicted hets
time python3 ../../scripts/compareHets.py trueMismatches.txt predictedMismatches.txt >validatedMismatches.txt
cat validatedMismatches.txt

if [ ${makePlots} = TRUE ]; then
  # Make POA visualizations of each difference
  echo Making POA graph visualizations of difference
  mkdir visualizations
  for i in $(cat ${1} | grep 'Mismatch:' | cut -f4 -d' '); do
    python3 ../../scripts/phasedPoaToDot.py ${2} --start=${i} --length=10 --output=visualizations/diff_${i}
  done

  # Make POA visualizations of each predicted het
  echo Making POA graph visualizations of predicted hets
  mkdir visualizations
  for i in $(cat validatedMismatches.txt | grep 'Predicted mismatch:' | cut -f5 -d' '); do
    j=$(echo ${i} - 5 | bc)
    python3 ../../scripts/phasedPoaToDot.py output_poa.csv.hap1 --start=${j} --length=10 --output=visualizations/hap1_poa_${j}
  done
fi

# Exit
cd ..
