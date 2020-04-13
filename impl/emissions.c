/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

/*
 * Allele alphabet and substitutions
 */

uint16_t *stSite_getSubstitutionProb(stSite *site, int64_t from, int64_t to) {
    /*
     * Gets the (log) substitution probability of getting the derived (to) allele
     * given the source (from/haplotype) allele.
     */
    return &site->substitutionLogProbs[from * site->alleleNumber + to];
}

/*
 * Emission probabilities with optimization to make fast
 */

/*
 * Following implement Hamming weight for uint64_t ints, taken from
 * https://en.wikipedia.org/wiki/Hamming_weight
 */

//types and constants used in the functions below
//uint64_t is an unsigned 64-bit integer variable type (defined in C99 version of C language)
const uint64_t m1 = 0x5555555555555555; //binary: 0101...
const uint64_t m2 = 0x3333333333333333; //binary: 00110011..
const uint64_t m4 = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...

inline int popcount64NoBuiltIn(uint64_t x) {
    /*
     * Returns Hamming weight of input unsigned integer.
     */
    //This uses fewer arithmetic operations than any other known
    //implementation on machines with fast multiplication.
    //This algorithm uses 12 arithmetic operations, one of which is a multiply.
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits
    return (x * h01) >> 56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}

inline int popcount64(uint64_t x) {
    /*
     * Returns Hamming weight of input unsigned integer.
     */
    //return popcount64NoBuiltIn(x); // Use this line if the builtin is unavailable
    return __builtin_popcountll(x);
}

inline int popcount128(uint128_t n) {
    /*
     *  128 bit popcount. Currently naive. TODO: planning to allow partitions to be up to depth 128
     */
    uint64_t parts[2];
    memcpy(parts, &n, sizeof(n));
    return popcount64(parts[0]) + popcount64(parts[1]);
}

static inline uint64_t *retrieveBitCountVector(uint64_t *bitCountVector,
                                               uint64_t siteOffset, uint64_t allele, uint64_t bit) {
    /*
     * Returns a pointer to a bit count vector for a given site, allele and bit.
     */
    return &bitCountVector[siteOffset * ALLELE_LOG_PROB_BITS
                           + allele * ALLELE_LOG_PROB_BITS
                           + bit];
}

uint64_t calculateBitCountVector(uint8_t **seqs, uint64_t depth,
                                 uint64_t siteOffset, uint64_t allele, uint64_t bit) {
    /*
     * Calculates the bit count vector for a given site, allele and bit.
     */
    uint64_t bitCountVector = 0;
    for (uint64_t i = 0; i < depth; i++) {
        uint8_t *p = &(seqs[i][siteOffset]);
        bitCountVector |= ((((uint64_t) p[allele] >> bit) & 1) << i);
    }

    return bitCountVector;
}

uint64_t *calculateCountBitVectors(uint8_t **seqs, stReference *ref,
                                   uint64_t firstSite, uint64_t length, uint64_t depth) {
    /*
     * Calculates the bit count vector for every site, allele and bit in the given range of sites.
     */

    if (ref->length == 0) {
        return st_malloc(0);
    }

    // Index of first (inclusive) and last (exclusive) allele in the column
    uint64_t firstAllele = ref->sites[firstSite].alleleOffset;
    uint64_t lastAllele =
            firstSite + length < ref->length ? ref->sites[firstSite + length].alleleOffset : ref->totalAlleles;

    // Array of bit vectors, for each site, for each allele and for each bit in uint8_t
    uint64_t *bitCountVectors = st_malloc((lastAllele - firstAllele) * ALLELE_LOG_PROB_BITS * sizeof(uint64_t));

    // For each site
    for (uint64_t i = firstSite; i < firstSite + length; i++) {
        // For each allele
        uint64_t siteOffset = ref->sites[i].alleleOffset - firstAllele;
        for (uint64_t j = 0; j < ref->sites[i].alleleNumber; j++) {
            // For each bit
            for (uint64_t k = 0; k < ALLELE_LOG_PROB_BITS; k++) {
                *retrieveBitCountVector(bitCountVectors, siteOffset, j, k) =
                        calculateBitCountVector(seqs, depth, siteOffset, j, k);
            }
        }
    }

    return bitCountVectors;
}

uint64_t getLogProbOfAllele(uint64_t *bitCountVectors, uint64_t depth, uint64_t partition,
                            uint64_t siteOffset, uint64_t allele) {
    /*
     * Returns the -log prob of the reads in a given partition being generated by a given allele.
     */
    uint64_t *j = retrieveBitCountVector(bitCountVectors, siteOffset, allele, 0);
    uint64_t negLogProb = popcount64(j[0] & partition);

    for (uint64_t i = 1; i < ALLELE_LOG_PROB_BITS; i++) {
        negLogProb += (popcount64(j[i] & partition) << i);
    }

    return negLogProb;
}

static inline uint64_t minP(uint64_t a, uint64_t b) {
    return a < b ? a : b;
}

static inline void alleleLogHapProbabilities(stRPColumn *column, stSite *site, uint64_t siteOffset,
                                             uint64_t partition, uint64_t *bitCountVectors,
                                             uint64_t *alleleLogProbs) {
    /*
     * For each allele calculate the -log probability of the
     * sub-partition.
     */
    for (uint64_t i = 0; i < site->alleleNumber; i++) {
        alleleLogProbs[i] = getLogProbOfAllele(bitCountVectors, column->depth, partition, siteOffset, i);
    }
}

static inline void ancestorHapProbabilities(stSite *site, uint64_t *alleleLogProbs,
                                            uint64_t *ancestorAlleleProbs) {
    /*
     * Get the -log probabilities of the ancestor alleles for a given set of allele log probabilities.
     */

    // Calculate the probability of haplotype alleles and read alleles for each ancestor allele
    for (uint64_t i = 0; i < site->alleleNumber; i++) {
        uint16_t *j = stSite_getSubstitutionProb(site, i, 0);
        ancestorAlleleProbs[i] = alleleLogProbs[0] + j[0];
        for (uint64_t k = 1; k < site->alleleNumber; k++) {
            ancestorAlleleProbs[i] =
                    minP(ancestorAlleleProbs[i],
                         alleleLogProbs[k] + j[k]);
        }
    }
}

static inline uint64_t getMaxAlleleLogProb(stSite *site, uint64_t *alleleLogProbs) {
    /*
     * Get the min -log probability of an allele.
     */

    // Calculate the probability of haplotype alleles and read alleles for each ancestor allele
    uint64_t p = alleleLogProbs[0];
    for (uint64_t i = 1; i < site->alleleNumber; i++) {
        p = minP(p, alleleLogProbs[i]);
    }
    return p;
}

static inline uint64_t genotypeLogProbability(stRPColumn *column, stSite *site, uint64_t siteOffset,
                                              uint64_t partition, uint64_t *bitCountVectors,
                                              bool includeAncestorSubProb) {
    /*
     * Get the -log probability of the alleles in a given position within a column for a given partition.
     */
    // Get the sum of log probabilities of the derived alleles over the possible source alleles

    // For each allele calculate the log probability of the
    // partition and store counts in an array
    uint64_t alleleLogProbsHap1[site->alleleNumber];
    alleleLogHapProbabilities(column, site, siteOffset, partition, bitCountVectors, alleleLogProbsHap1);
    uint64_t alleleLogProbsHap2[site->alleleNumber];
    alleleLogHapProbabilities(column, site, siteOffset, ~partition, bitCountVectors, alleleLogProbsHap2);

    if (!includeAncestorSubProb) {
        return getMaxAlleleLogProb(site, alleleLogProbsHap1) + getMaxAlleleLogProb(site, alleleLogProbsHap2);
    }

    uint64_t ancestorAlleleProbsHap1[site->alleleNumber];
    ancestorHapProbabilities(site, alleleLogProbsHap1, ancestorAlleleProbsHap1);
    uint64_t ancestorAlleleProbsHap2[site->alleleNumber];
    ancestorHapProbabilities(site, alleleLogProbsHap2, ancestorAlleleProbsHap2);

    // Combine the probabilities to calculate the overall probability of the genotype
    uint64_t logGenotypeProb = ancestorAlleleProbsHap1[0] + ancestorAlleleProbsHap2[0] + site->allelePriorLogProbs[0];
    for (uint64_t i = 1; i < site->alleleNumber; i++) {
        logGenotypeProb = minP(logGenotypeProb,
                               ancestorAlleleProbsHap1[i] + ancestorAlleleProbsHap2[i] + site->allelePriorLogProbs[i]);
    }

    return logGenotypeProb;
}

double emissionLogProbability(stRPColumn *column,
                              stRPCell *cell, uint64_t *bitCountVectors, stReference *ref,
                              stRPHmmParameters *params) {
    /*
     * Get the log probability of a set of reads for a given column.
     */
    assert(column->length > 0);
    uint64_t logPartitionProb = 0;
    uint64_t firstAllele = ref->sites[column->refStart].alleleOffset;
    for (uint64_t i = column->refStart; i < column->refStart + column->length; i++) {
        stSite *site = &(ref->sites[i]);
        uint64_t siteOffset = site->alleleOffset - firstAllele;

        // Get the reference prior probabilities
        logPartitionProb += genotypeLogProbability(column, site, siteOffset, cell->partition, bitCountVectors,
                                                   params->includeAncestorSubProb);
    }

    return -((double) logPartitionProb);
}

/*
 * Functions for calculating genotypes/haplotypes
 */

static uint64_t getMLAllele(stSite *site, uint64_t *alleleLogProbs, uint64_t maxProbAncestorAllele) {
    /*
     * Return the allele with maximum probability for a given ancestor allele.
     */
    uint64_t maxAllele = 0;
    uint64_t maxProb = alleleLogProbs[0] + *stSite_getSubstitutionProb(site, maxProbAncestorAllele, 0);
    for (uint64_t i = 1; i < site->alleleNumber; i++) {
        uint64_t hapProb = alleleLogProbs[i] + *stSite_getSubstitutionProb(site, maxProbAncestorAllele, i);
        if (hapProb < maxProb) {
            maxProb = hapProb;
            maxAllele = i;
        }
    }
    return maxAllele;
}

static void fillInPredictedGenomePosition(stGenomeFragment *gF, uint64_t siteIndex, uint64_t partition,
                                          stRPColumn *column, uint64_t *bitCountVectors) {
    /*
     * Computes the most probable haplotype alleles / genotype and associated posterior
     * probabilities for a given position within a cell/column.
     */

    // Coordinates
    stSite *site = &(gF->reference->sites[siteIndex]);
    uint64_t firstAllele = gF->reference->sites[column->refStart].alleleOffset;
    uint64_t siteOffset = site->alleleOffset - firstAllele;

    // Get the log-prob of ancestral alleles

    // For each allele calculate the log probability of the
    // partition and store counts in an array
    uint64_t alleleLogProbsHap1[site->alleleNumber];
    alleleLogHapProbabilities(column, site, siteOffset, partition, bitCountVectors, alleleLogProbsHap1);
    uint64_t alleleLogProbsHap2[site->alleleNumber];
    alleleLogHapProbabilities(column, site, siteOffset, ~partition, bitCountVectors, alleleLogProbsHap2);

    uint64_t ancestorAlleleProbsHap1[site->alleleNumber];
    ancestorHapProbabilities(site, alleleLogProbsHap1, ancestorAlleleProbsHap1);
    uint64_t ancestorAlleleProbsHap2[site->alleleNumber];
    ancestorHapProbabilities(site, alleleLogProbsHap2, ancestorAlleleProbsHap2);

    // Combine the probabilities to calculate the log prob of the ml ancestor allele
    uint64_t maxLogColumnProb = ancestorAlleleProbsHap1[0] + ancestorAlleleProbsHap2[0] + site->allelePriorLogProbs[0];
    uint64_t ancestorAllele = 0;
    for (uint64_t i = 1; i < site->alleleNumber; i++) {
        uint64_t j = ancestorAlleleProbsHap1[i] + ancestorAlleleProbsHap2[i] + site->allelePriorLogProbs[i];
        if (j < maxLogColumnProb) {
            maxLogColumnProb = j;
            ancestorAllele = i;
        }
    }

    // Given ml ancestor allele, figure prob of haplotype alleles
    uint64_t hapAllele1 = getMLAllele(site, alleleLogProbsHap1, ancestorAllele);
    uint64_t hapAllele2 = getMLAllele(site, alleleLogProbsHap2, ancestorAllele);

    uint64_t k = siteIndex - gF->refStart;
    // Fill in the genome fragment
    gF->ancestorString[k] = ancestorAllele;
    gF->haplotypeString1[k] = hapAllele1;
    gF->haplotypeString2[k] = hapAllele2;
    // Get combined genotype
    gF->genotypeString[k] = hapAllele1 < hapAllele2 ? hapAllele1 * site->alleleNumber + hapAllele2 :
                            hapAllele2 * site->alleleNumber + hapAllele1;
    gF->genotypeProbs[k] = -((float) maxLogColumnProb);
    gF->haplotypeProbs1[k] = -(float) alleleLogProbsHap1[hapAllele1];
    gF->haplotypeProbs2[k] = -(float) alleleLogProbsHap2[hapAllele2];

    // Fill in read coverages
    assert(column->depth >= popcount64(partition));
    gF->readsSupportingHaplotype1[k] = popcount64(partition);
    gF->readsSupportingHaplotype2[k] = column->depth - popcount64(partition);
    //st_uglyf(" Hello %i %i %i %f %f \n", (int)popcount64(partition), (int)column->depth, (int)column->depth - popcount64(partition),
    //		(float)binomialPValue(column->depth, popcount64(partition)), (float)binomialPValue(column->depth, column->depth-popcount64(partition)));
}

void fillInPredictedGenome(stGenomeFragment *gF, uint64_t partition,
                           stRPColumn *column, stRPHmmParameters *params) {
    /*
     * Computes the most probable haplotype alleles / genotypes and associated posterior
     * probabilities for a given interval defined by a cell/column. Fills in these values in the
     * genome fragment argument.
     */

    // Calculate the bit vectors
    uint64_t *bitCountVectors = calculateCountBitVectors(column->seqs, gF->reference, column->refStart,
                                                         column->length, column->depth);

    assert(column->length > 0);
    for (uint64_t i = 0; i < column->length; i++) {
        fillInPredictedGenomePosition(gF, i + column->refStart, partition, column,
                                      bitCountVectors);
    }

    // Cleanup
    free(bitCountVectors);
}
