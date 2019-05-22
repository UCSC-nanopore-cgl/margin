# MarginPolish #

MarginPolish is a graph-based assembly polisher that takes advantage of multiple probable alignment paths in run-length space.  It takes as input a FASTA assembly and BAM, and it produces a polished FASTA assembly.  

While MarginPolish serves as a standalone assembly polisher, it is also  part of an assembly pipeline which includes an ultrafast nanopore assembler [Shasta](https://github.com/chanzuckerberg/shasta) and a multi-task RNN polisher [HELEN](https://github.com/kishwarshafin/helen).  HELEN consumes images generated by MarginPolish.

## Overview ##

MarginPolish takes reads and alignments from the BAM and generates an initial graph describing matches, inserts, and deletes observed in the alignments.  It aligns each read to this graph and determines the probabilities for multiple likely alignments.  After all reads are aligned, it generates weighted alignment scores for each node in the graph and uses these to determine a most-likely path through the graph.  This path becomes the inital graph for the next iteration of the process.  Iteration continues until the total weighted likelihood of alignments decreases between iteration steps, or until a configured maximum iteration count is reached.

All of these alignments are done in [run-length space](https://en.wikipedia.org/wiki/Run-length_encoding), which simplifies and cleans the alignments.  This reduces effects from errors in homopolymer runs, which are the primary source of error in nanopore reads.  After determining a most-likely run-length-encoded sequence, it decompresses the run-lengths to generate a final consensus sequence.  This expansion is done using a Bayesian model which predicts the most likely true run-length from all run-length observations in the reads aligned to the node.


### Execution Flow ###

MarginPolish first determines where on the initial assembly segments there are read alignments.  It uses these positions to determine chunk boundaries, including a configurable overlap between chunks.  Each chunk is polished separately and the final consensus sequence is stored in memory.  After all chunks are complete, MarginPolish stitches them together by aligning the overlap between chunks and finding a position to stitch together.  

### Assembly Workflow ###

For a detailed description of the end-to-end assembly workflow, see the documentation provided in the [HELEN](https://github.com/kishwarshafin/helen) repository.

## Installation ##

### Dependencies ###
cmake version 3.7 (or higher):
```
wget https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.sh && mkdir /opt/cmake && sh cmake-3.7.2-Linux-x86_64.sh --prefix=/opt/cmake --skip-license && ln -s /opt/cmake/bin/cmake /usr/local/bin/cmake
```

Many of these dependencies are due to integration with htslib.  If compiling on Ubuntu, this will install all required packages:
```
apt-get -y install make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libhdf5-dev
```

### Compilation ###

- Check out the repository and submodules:
```
git clone https://github.com/UCSC-nanopore-cgl/marginPolish.git
cd marginPolish
git submodule update --init
```

- Make build directory:
```
mkdir build
cd build
```

- Generate Makefile and run:
```
cmake ..
make
./marginPolish
 ```

## Running MarginPolish ##


- to run marginPolish:
``` marginPolish <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options] ```

- program OPTIONS:
```
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
    -r --region              : If set, will only compute for given chromosomal region.
                                 Format: chr:start_pos-end_pos (chr3:2000-3000).

HELEN feature generation options:
    -f --produceFeatures     : output features for HELEN.
    -F --featureType         : output features of chunks for HELEN.  Valid types:
                                 splitRleWeight:  [default] run lengths split into chunks
                                 nuclAndRlWeight: split into nucleotide and run length (RL across nucleotides)
                                 rleWeight:       weighted likelihood from POA nodes (RLE)
                                 simpleWeight:    weighted likelihood from POA nodes (non-RLE)
    -L --splitRleWeightMaxRL : max run length (for 'splitRleWeight' type only) [default = 10]
    -u --trueReferenceBam    : true reference aligned to ASSEMBLY_FASTA, for HELEN
                               features.  Setting this parameter will include labels
                               in output.

Miscellaneous supplementary output options:
    -i --outputRepeatCounts  : Output base to write out the repeat counts [default = NULL]
    -j --outputPoaTsv        : Output base to write out the poa as TSV file [default = NULL]
```


#### Sample Execution

```./marginPolish ../tests/NA12878.np.chr3.5kb.bam ../tests/hg19.chr3.9mb.fa ../params/allParams.np.human.json -o example_out```


### Configuration ###

MarginPolish uses a configuration file which contains parameters for:
- Alignment
- Chunking
- Read Filtering
- Run Length Management

We find that the run-length estimation correctness is strongly tied to the basecaller version that was used


### HELEN Image Generation ###




MarginPolish produces different image types (used during development) which can be configured with the -F flag, but we encourage users to only use the default type 'splitRleWeight' for use with HELEN.


### Miscellaneous Output ###

The -i flag will output a TSV file for each chunk describing the observed run lengths at each node in the final alignment.  This can be used to train a Bayesian model which can predict run lengths.  The -j flag will output a representation of the POA for each chunk.

### MarginPhase ###

MarginPhase is a program for simultaneous haplotyping and genotyping with long reads.  It relies on much of the same infrastructure as marginPolish.  We refer you to the marginPhase repo (from which this is forked) at:

https://github.com/benedictpaten/marginPhase

There is infrastructure in this codebase which is used in support of MarginPhase and is not directly tied to MarginPolish.  It is our intention to eventually integrate these two applications to produce polished haplotypes instead of a single consensus. 


### Tests ###

After building, run the 'allTests' executable in your build directory.  This runs every test. You can comment out ones you don't want to run in tests/allTests.c
