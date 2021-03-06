# Margin #

**Margin** is a suite of tools used for analysis of long-read data using Hidden Markov Models.  It has two submodules: `polish` and `phase`.

**MarginPhase** is a utility for read haplotyping and variant phasing.
It takes a BAM, a VCF, and a reference FASTA as input, and it produces a BAM with haplotagged reads and a VCF with phased variants.
Margin's phase submodule is used in the haplotyping and variant phasing steps of the pipeline [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper/).

**MarginPolish** is a diploid-aware assembly polisher. 
It takes as input a FASTA assembly and an indexed BAM (ONT reads aligned to the assembly), and it produces a haploid or diploid polished FASTA assembly. 
While Margin's polish submodule serves as a standalone assembly polisher, it is also part of an assembly pipeline which includes an ultrafast nanopore assembler [Shasta](https://github.com/chanzuckerberg/shasta) and a multi-task RNN polisher [HELEN](https://github.com/kishwarshafin/helen). 
HELEN operates on images generated by Margin to generate a polished haploid sequence.

## Installation ##

### Dependencies ###

If compiling on Ubuntu, this will install all required packages:
```
apt-get -y install git make gcc g++ autoconf zlib1g-dev libcurl4-openssl-dev libbz2-dev libhdf5-dev
```

Note that libhdf5-dev is required for HELEN image generation with MarginPolish.  Both submodules will work without this package, but will not include image generation functionality for HELEN.

Margin is compiled with cmake.  We recommend using the latest cmake version, but 3.7 and higher are supported:
```
wget https://github.com/Kitware/CMake/releases/download/v3.14.4/cmake-3.14.4-Linux-x86_64.sh && sudo mkdir /opt/cmake && sudo sh cmake-3.14.4-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && sudo ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
cmake --version
```

### Compilation ###

```
# Check out the repository and submodules:
git clone https://github.com/UCSC-nanopore-cgl/margin.git
cd margin
git submodule update --init

# Make build directory:
mkdir build
cd build

# Generate Makefile and run:
cmake ..
make
./margin
 ```

### Tests ###

Unit tests are contained in the 'allTests' executable, which can be run in your build directory. 
This runs every test and can take up to an hour. 
To test individual components, you can comment out ones you don't want to run in tests/allTests.c.

## Running Margin ##

### Data Formats ###

Margin requires that the input BAM is indexed.  Both submodules can accept read information both in FASTA and FASTQ formats (aligned).  If quality scores are present, the base likelihood is factored into alignment weight estimation.  CRAM format is currently unsupported.

Margin requires an unzipped FASTA reference.  It uses an index to access regions of the FASTA, but can generate the index without requiring the user to if the process has write access to the data directory.

Margin supports both gzipped and uncompressed VCFs.  No index is needed.

### Configuration ###

Margin uses a JSON configuration file which contains model and runtime parameterization for:
- Pairwise Alignment
- Graph Alignment
- Run Length Estimation
- Chunking
- Read/Depth Management

Both tools have specific parameter files for read data and output goals, all of which extend from a base parameters file containing shared configuration.

For phasing, there are two parameter files for whether your primary objective is to haplotag reads or to phase variants, found in the `misc` folder below.

For polishing, the run-length estimation correctness is tied to the basecaller version that was used to generate the reads. 
We have trained models for the ONT r9.4 and r10.3 pores with multiple versions of the Guppy basecaller for human and microbial samples, provided in the `ont` folder. 
Additionally we have provided a model trained on PacBio HiFi reads at `pacbio/hifi/allParams.hifi.json`.

```
params/
├── base_params.json
├── misc
│   ├── allParams.ccs_haplotag.json
│   ├── allParams.no_rle.json
│   ├── allParams.ont_haplotag.json
│   └── allParams.phase_vcf.json
├── ont
│   ├── r10.3
│   │   ├── allParams.np.human.r103-g3210.json
│   │   └── allParams.np.microbial.r103g324.json
│   └── r9.4
│       ├── allParams.np.human.r94-g235.json
│       ├── allParams.np.human.r94-g305.json
│       ├── allParams.np.human.r94-g344.json
│       ├── allParams.np.human.r94-g360.json
│       ├── allParams.np.microbial.r94-g305.json
│       └── allParams.np.microbial.r94-g344.json
└── pacbio
    └── hifi
        └── allParams.hifi.json
```
### Chunking ###

Margin first determines where there are read alignments on the initial assembly segments.  It uses these positions to determine chunk boundaries, including a configurable overlap between chunks. 
Each chunk is handled separately (by a single thread) and the final result (phasing information or consensus sequence) is stored in memory. 
After all chunks are complete, Margin stitches them together by comparing read haplotype assignments (phase), or by aligning the overlap between chunks and finding a position to stitch at (polish).  


## MarginPhase ##

Margin's phase submodule is a tool developed primarily to assist in the long-read variant calling pipeline PEPPER-Margin-DeepVariant.
Margin serves two uses in this capacity: haplotagging reads to assist with genotyping, and phasing variants in the final genotyped VCF.
For both use cases, MarginPhase requires a BAM, reference, and VCF.

For each chunk, Margin extracts read substrings aligned a configurable distance up- and downstream from each variant position.
Reads are downsampled to a configurable depth prioritizing longer reads.
Alignment likelihoods are generated between all read substrings and the alleles (generated by mutating the reference sequence). 
These likelihoods are fed into the phasing algorithm to haplotype alleles. 
Any filtered read substrings (during initial extraction from the BAM or during downsampling) are aligned to the two haplotypes, and are assigned based on which haplotype they align best to.
Stitching is performed before final results are output.

During phaseset determination, adjacent variants are phased together unless there are no reads spanning adjacent variants,
if the ratio of the number of reads in trans to the number of reads in cis is above a threshold, 
or if the likelihood of the read partitioning between haplotypes at the variant site is below a threshold.
The reason for a phase break is documented in the `OUTPUT.phaseset.bed` file.

### Variant Calling Workflow ###

The PEPPER-Margin-DeepVariant pipeline is documented [here](https://github.com/kishwarshafin/pepper) .

### Running MarginPhase ###

```
$ ./margin phase
usage: margin phase <ALIGN_BAM> <REFERENCE_FASTA> <VARIANT_VCF> <PARAMS> [options]
Version: 2.0

Tags reads in ALIGN_BAM using variants in VARIANT_VCF.

Required arguments:
    ALIGN_BAM is the alignment of reads to the reference.
    REFERENCE_FASTA is the reference sequence BAM file in fasta format.
    VARIANT_VCF is the set of variants to use for phasing.
    PARAMS is the file with margin parameters.

Default options:
    -h --help                : Print this help screen
    -a --logLevel            : Set the log level [default = info]
    -t --threads             : Set number of concurrent threads [default = 1]
    -o --outputBase          : Name to use for output files [default = 'output']
    -r --region              : If set, will only compute for given chromosomal region
                                 Format: chr:start_pos-end_pos (chr3:2000-3000)
    -p --depth               : Will override the downsampling depth set in PARAMS
    -k --tempFilesToDisk     : Write temporary files to disk (for --diploid or supplementary output)

Output options:
    -M --skipHaplotypeBAM    : Do not write out phased BAM
    -V --skipPhasedVCF       : Do not write out phased VCF
```

### Sample Execution ###
```
# haplotagging mode
./margin phase \
        /path/to/margin/tests/data/realData/HG002.r94g360.chr20_59M_100k.bam \
        /path/to/margin/tests/data/realData/hg38.chr20_59M_100k.fa \
        /path/to/margin/tests/data/realData/HG002.r94g360.chr20_59M_100k.vcf \
        /path/to/margin/tests/data/realData/params/misc/allParams.ont_haplotag.json \
        --skipPhasedVCF
                
# variant phasing mode
./margin phase \
        /path/to/margin/tests/data/realData/HG002.r94g360.chr20_59M_100k.bam \
        /path/to/margin/tests/data/realData/hg38.chr20_59M_100k.fa \
        /path/to/margin/tests/data/realData/HG002.r94g360.chr20_59M_100k.vcf \
        /path/to/margin/tests/data/realData/params/misc/allParams.ont_phase_vcf.json \
        --skipHaplotypeBam
```

### Resource Requirements ###

Margin took between 19m (35x PacBio-HiFi) and 80m (75x ONT) for a single whole genome with 64 threads and peak memory usage of 35GB during phasing.




## MarginPolish ##

Margin's polish submodule takes reads and alignments from the BAM and generates an initial graph describing matches, inserts, and deletes observed in the alignments.  It aligns each read to this graph and determines the probabilities for multiple likely alignments.  After all reads are aligned, it generates weighted alignment scores for each node in the graph and uses these to determine a most-likely path through the graph.  This path becomes the inital graph for the next iteration of the process.  Iteration continues until the total weighted likelihood of alignments decreases between iteration steps, or until a configured maximum iteration count is reached.

All of these alignments are done in [run-length space](https://en.wikipedia.org/wiki/Run-length_encoding), which simplifies and cleans the alignments.  This reduces effects from errors in homopolymer runs, which are the primary source of error in nanopore reads.  After determining a most-likely run-length-encoded sequence, it expands the run-lengths to generate a final consensus sequence.  This expansion is done using a Bayesian model which predicts the most likely true run-length from all run-length observations in the reads aligned to the node.

In diploid polishing mode, Margin identifies candidate sites where variants may be present (alternatively, a VCF can be supplied and the tool will use these sites instead of its own detection strategy). 
It finds read substrings aligned around these sites, uses all unique read substrings as candidate alleles, and aligns all read substrings to all candidate alleles.
These alleles and the read alignments to them are used to construct a graph containing read linkage information between sites. 
The phasing HMM is applied to this graph, and we use this to identify the most likely alleles at each site and the haplotype to which they belong. 
After all chunks are complete, Margin stitches haplotypes together using a set similarity metric between reads assigned to the adjacent chunks' haplotypes, and stitches sequences together using the haploid stitching method.

### Assembly Workflow ###

For a detailed description of the end-to-end assembly workflow, see the documentation provided in the [HELEN](https://github.com/kishwarshafin/helen) repository.


### Running MarginPolish ###

``` 
$ ./margin polish
usage: margin polish <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options]
Version: 2.0

Polishes the ASSEMBLY_FASTA using alignments in BAM_FILE.

Required arguments:
    BAM_FILE is the alignment of reads to the assembly (or reference).
    ASSEMBLY_FASTA is the reference sequence BAM file in fasta format.
    PARAMS is the file with marginPolish parameters.

Default options:
    -h --help                : Print this help screen
    -a --logLevel            : Set the log level [default = info]
    -t --threads             : Set number of concurrent threads [default = 1]
    -o --outputBase          : Name to use for output files [default = 'output']
    -r --region              : If set, will only compute for given chromosomal region
                                 Format: chr:start_pos-end_pos (chr3:2000-3000)
    -p --depth               : Will override the downsampling depth set in PARAMS
    -k --tempFilesToDisk     : Write temporary files to disk (for --diploid or supplementary output)

Diploid options:
    -2 --diploid             : Will perform diploid phasing.
    -v --vcf                 : VCF with sites for phasing (will not perform variant detection if set)
    -S --skipFilteredReads   : Will NOT attempt to haplotype filtered reads (--diploid only)
    -R --skipRealignment     : Skip realignment (for haplotyping only)
    -A --onlyVcfAlleles      : Only use alleles specified in the VCF. Requires NO RLE and 
                                 Requires NO RLE and --skipOutputFasta

HELEN feature generation options:
    -f --produceFeatures     : output splitRleWeight or diploidRleWeight (based on -2 flag) features for HELEN
    -F --featureType         : output specific feature type for HELEN (overwrites -f).  Valid types:
                                 splitRleWeight:   [default] run lengths split into chunks
                                 channelRleWeight: run lengths split into per-nucleotide channels
                                 simpleWeight:     weighted likelihood from POA nodes (non-RLE)
                                 diploidRleWeight: [default] produces diploid features 
    -L --splitRleWeightMaxRL : max run length (for RLE feature types) 
                                 [split default = 10, channel default = 10]
    -u --trueReferenceBam    : true reference aligned to ASSEMBLY_FASTA, for HELEN
                               features.  Setting this parameter will include labels
                               in output.  If -2/--diploid is set, this parameter must
                               contain two comma-separated values

Miscellaneous supplementary output options:
    -c --supplementaryChunks : Write supplementary files for each chunk (in additon to writing
                               whole genome information)
    -d --outputPoaDot        : Write out the poa as DOT file (only done per chunk)
    -i --outputRepeatCounts  : Write out the repeat counts as CSV file
    -j --outputPoaCsv        : Write out the poa as CSV file
    -n --outputHaplotypeReads: Write out phased reads and likelihoods as CSV file (--diploid only)
    -s --outputPhasingState  : Write out phasing likelihoods as JSON file (--diploid only)
    -M --skipHaplotypeBAM    : Do not write out phased BAMs (--diploid only, default is to write)
    -T --skipOutputFasta     : Do not write out phased fasta (--diploid only, default is to write)
```

#### Sample Execution ####

```
# haploid mode
./margin polish \
        /path/to/margin/tests/data/realData/HG002.r94g360.chr20_59M_100k.bam \
        /path/to/margin/tests/data/realData/hg38.chr20_59M_100k.fa \
        /path/to/margin/tests/data/realData/params/ont/r9.4/allParams.np.human.r94-g360.json
        
# diploid mode
./margin polish \
        /path/to/margin/tests/data/realData/HG002.r94g360.chr20_59M_100k.bam \
        /path/to/margin/tests/data/realData/hg38.chr20_59M_100k.fa \
        /path/to/margin/params/ont/r9.4/allParams.np.human.r94-g360.json \
        --diploid
        -v /path/to/margin/tests/data/realData/HG002.r94g360.chr20_59M_100k.vcf
```


### HELEN Image Generation ###

[HELEN](https://github.com/kishwarshafin/helen) is a multi-task RNN polisher which operates on images produced by MarginPolish.  The images summarize the state of the nodes in the alignment graph before run-length expansion.  They include the weights associated with read observations aligned at each node.

If MarginPolish is configured to generate images (the -f option), it will output a single .h5 file for each thread. 

MarginPolish produces different image types (used during development) which can be configured with the -F flag, but users must use the default type 'splitRleWeight' for the trained models HELEN provides.

To produce HELEN training images, run MarginPolish with the -u flag. 
This takes an argument of an indexed BAM alignment of the truth sequence to the assembly. 
MarginPolish will extract the alignments from this BAM for each analyzed chunk. 
If there is a single alignment at this location with a sequence that approximately matches the chunk's size, it is used to label the images for both nucleotide and run-length. 

### Resource Requirements ###

While comprehensive resource usage profiling has not been done yet, we find that memory usage scales linearly with thread count, read depth, and chunk size. 
For this reason, our default parameters downsample read depth to 64 and restrict chunk size to 100000 bases.

We found that 2GB of memory per thread is sufficient to run Margin's polish submodule in haploid mode on genome-scale assemblies and alignment.

Across 13 whole-genome runs, we averaged roughly 350 CPU hours per gigabase of assembled sequence.



© 2019 by Benedict Paten (benedictpaten@gmail.com), Trevor Pesout (tpesout@ucsc.edu)
