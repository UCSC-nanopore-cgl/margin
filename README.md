# Margin #

**Margin** is a suite of tools used for analysis of long-read data using Hidden Markov Models.  It has two submodules: `polish` and `phase`.

**MarginPhase** is a utility for read haplotyping and variant phasing.  

**MarginPolish** is a diploid-aware assembly polisher.  Documentation is available [here](docs/MarginPolish.md).

## Quickstart ##

The PEPPER-Margin-DeepVariant pipeline is documented [here](https://github.com/kishwarshafin/pepper).

Margin as a standalone tool is most easily used with a Docker image available here:
```
docker pull kishwars/pepper_deepvariant:r0.4
docker run kishwars/pepper_deepvariant:r0.4 margin phase -h
```

Margin requires an indexed BAM, a reference FASTA, and a VCF (supports .vcf and .vcf.gz)

To haplotag ONT reads:
```
docker run \
    -v `pwd`:/data \
    kishwars/pepper_deepvariant:r0.4 \
    margin phase \
    /data/$YOUR_ALIGNMENT_HERE.bam \
    /data/$YOUR_REFERENCE_HERE.fasta \
    /data/$YOUR_VARIANTS_HERE.vcf \
    /opt/margin_dir/params/misc/allParams.ont_haplotag.json \
    -t $THREAD_COUNT \
    -o /data/$OUTPUT_PREFIX \
    --skipPhasedVCF
```

To haplotag PacBio-HiFi reads:
```
docker run \
    -v `pwd`:/data \
    kishwars/pepper_deepvariant:r0.4 \
    margin phase \
    /data/$YOUR_ALIGNMENT_HERE.bam \
    /data/$YOUR_REFERENCE_HERE.fasta \
    /data/$YOUR_VARIANTS_HERE.vcf \
    /opt/margin_dir/params/misc/allParams.ccs_haplotag.json \
    -t $THREAD_COUNT \
    -o /data/$OUTPUT_PREFIX \
    --skipPhasedVCF
```

To phase a VCF:
```
docker run \
    -v `pwd`:/data \
    kishwars/pepper_deepvariant:r0.4 \
    margin phase \
    /data/$YOUR_ALIGNMENT_HERE.bam \
    /data/$YOUR_REFERENCE_HERE.fasta \
    /data/$YOUR_VARIANTS_HERE.vcf \
    /opt/margin_dir/params/misc/allParams.phase_vcf.json \
    -t $THREAD_COUNT \
    -o /data/$OUTPUT_PREFIX \
    --skipHaplotypeBAM 
```

## Installation ##

### Dependencies ###

If compiling on Ubuntu, this will install all required packages:
```
sudo apt-get install git make gcc g++ autoconf zlib1g-dev libcurl4-openssl-dev libbz2-dev libhdf5-dev
```

Note that libhdf5-dev is required for HELEN image generation with MarginPolish but is not required. 

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
 
### Verification ###
```
# haplotagging mode
./margin phase \
    ../tests/data/realData/HG002.r94g360.chr20_59M_100k.bam \
    ../tests/data/realData/hg38.chr20_59M_100k.fa \
    ../tests/data/realData/HG002.r94g360.chr20_59M_100k.vcf \
    ../params/misc/allParams.ont_haplotag.json \
    --skipPhasedVCF
# verification: expect 145, 137
samtools view output.haplotagged.bam | grep "HP:i:1" | wc -l
samtools view output.haplotagged.bam | grep "HP:i:2" | wc -l
                
# variant phasing mode
./margin phase \
    ../tests/data/realData/HG002.r94g360.chr20_59M_100k.bam \
    ../tests/data/realData/hg38.chr20_59M_100k.fa \
    ../tests/data/realData/HG002.r94g360.chr20_59M_100k.vcf \
    ../params/misc/allParams.phase_vcf.json \
    --skipHaplotypeBAM
# verification: expect 105
cat output.phased.vcf | grep -e "1|0\|0|1" | wc -l
```

### Resource Requirements ###

Margin took between 19m (35x PacBio-HiFi) and 80m (75x ONT) for a single whole genome with 64 threads and peak memory usage of 35GB during phasing.

© 2019 by Benedict Paten (benedictpaten@gmail.com), Trevor Pesout (tpesout@ucsc.edu)
