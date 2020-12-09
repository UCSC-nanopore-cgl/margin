/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

//#include <getopt.h>
//#include <stdio.h>
//#include <ctype.h>
//#include <memory.h>
//#include <hashTableC.h>
//#include <unistd.h>
//#include <time.h>
#include <math.h>
#include "marginVersion.h"

#include "margin.h"
#include "localPhasingCorrectness.h"
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>


PhasedVariant *phasedVariant_construct(const char *refSeqName, int64_t refPos, double quality, stList *alleles, int64_t gt1, int64_t gt2, char * phaseSet) {
    PhasedVariant *pv = st_calloc(1, sizeof(PhasedVariant));
    pv->refSeqName = stString_copy(refSeqName);
    pv->refPos = refPos;
    pv->quality = quality;
    pv->alleles = alleles;
    pv->gt1 = gt1;
    pv->gt2 = gt2;
    pv->phaseSet = stString_copy(phaseSet);
    return pv;
}

void phasedVariant_destruct(PhasedVariant *pv) {
    free(pv->refSeqName);
    stList_destruct(pv->alleles);
    free(pv->phaseSet);
    free(pv);
}
int phasedVariant_positionCmp(const void *a, const void *b) {
    PhasedVariant *A = (PhasedVariant*) a;
    PhasedVariant *B = (PhasedVariant*) b;
    if (A->refPos == B->refPos) {
        st_logCritical("Encountered two variants at same position: %s:%"PRId64"\n", A->refSeqName, A->refPos);
        return 0;
    }
    return A->refPos < B->refPos ? -1 : 1;
}

stHash *getPhasedVariants(const char *vcfFile) {
    // what we're saving into
    stHash *entries = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, (void(*)(void*))stList_destruct);
    time_t start = time(NULL);

    // open file
    htsFile *fp = hts_open(vcfFile,"rb");
    if (fp == NULL) {
        st_errAbort("Could not open VCF %s\n", vcfFile);
    }

    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    int nsmpl = bcf_hdr_nsamples(hdr);
    if (nsmpl > 1) {
        st_logCritical("Got %d samples reading %s, will only take VCF records for the first\n", nsmpl, vcfFile);
    }

    // find type of phaseSet
    bool phaseSetIsInt = FALSE;
    int psId = bcf_hdr_id2int(hdr, BCF_DT_ID, "PS");
    if (psId < 0) {
        st_errAbort("PS tag not present in VCF header for %s", vcfFile);
    }
    int psType = bcf_hdr_id2type(hdr, BCF_HL_FMT, psId);
    if (psType == BCF_HT_INT) {
        phaseSetIsInt = TRUE;
    } else if (psType == BCF_HT_STR) {
        phaseSetIsInt = FALSE;
    } else {
        st_errAbort("Unknown PS type in VCF header for %s", vcfFile);
    }


    // tracking
    int64_t totalEntries = 0;
    int64_t skippedForNotPass = 0;
    int64_t skippedForHomozygous = 0;
    int64_t skippedForNoPhaseset = 0;
    int64_t totalSaved = 0;

    // iterate over records
    bcf1_t *rec = bcf_init();
    while ( bcf_read(fp, hdr, rec) >= 0 )
    {
        //unpack for read REF,ALT,INFO,etc
        bcf_unpack(rec, BCF_UN_ALL);
        totalEntries++;

        // pass variant
        if (!bcf_has_filter(hdr, rec, "PASS")) {
            skippedForNotPass++;
            continue;
        }

        // genotype
        int gt1 = -1;
        int gt2 = -1;
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if (ngt>0 && !bcf_gt_is_missing(gt_arr[0])  && gt_arr[1] != bcf_int32_vector_end) {
            gt1 = bcf_gt_allele(gt_arr[0]);
            gt2 = bcf_gt_allele(gt_arr[1]);
        }
        free(gt_arr);
        if (gt1 == gt2) {
            skippedForHomozygous++;
            continue;
        }

        // phase set
        char *phaseset = NULL;
        if (phaseSetIsInt) {
            int mPSs = 0, nPSs;
            int32_t **PSs = NULL;
            nPSs = bcf_get_format_int32(hdr, rec, "PS", &PSs, &mPSs);
            if (nPSs <= 0 || PSs[0] == 0) {
                skippedForNoPhaseset++;
                continue;
            }
            phaseset = stString_print("%"PRId32, PSs[0]);
        } else {
            int mPSs = 0, nPSs;
            char **PSs = NULL;
            nPSs = bcf_get_format_string(hdr, rec, "PS", &PSs, &mPSs);
            if (nPSs <= 0 || stString_eq(PSs[0], ".")) {
                skippedForNoPhaseset++;
                continue;
            }
            phaseset = stString_copy(PSs[0]);
        }


        // location data
        const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
        int64_t pos = rec->pos;

        // qual
        double quality = rec->qual;

        // get alleles
        stList *alleles = stList_construct3(0, (void (*)(void*)) free);
        for (int i=0; i<rec->n_allele; ++i) {
            stList_append(alleles, stString_copy(rec->d.allele[i]));
        }

        // save it
        PhasedVariant *pv = phasedVariant_construct(chrom, pos, quality, alleles, gt1, gt2, phaseset);
        stList *contigList = stHash_search(entries, pv->refSeqName);
        if (contigList == NULL) {
            contigList = stList_construct3(0, (void(*)(void*))phasedVariant_destruct);
            stHash_insert(entries, stString_copy(pv->refSeqName), contigList);
        }
        stList_append(contigList, pv);
        totalSaved++;
    }

    // cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) ) {
        st_logCritical("> Failed to close VCF %s with code %d\n", vcfFile, ret);
    }

    // loggit
    st_logInfo("Read %"PRId64" variants from %s over %"PRId64" contigs in %"PRId64"s, keeping %"PRId64" phased variants"
               " and discarding %"PRId64" for not PASS, %"PRId64" for HOM, %"PRId64" for not phased.\n",
            totalEntries, vcfFile, stHash_size(entries), time(NULL) - start, totalSaved, skippedForNotPass,
            skippedForHomozygous, skippedForNoPhaseset);

    // ensure sorted
    stHashIterator *itor = stHash_getIterator(entries);
    char *contigName = NULL;
    while ((contigName = stHash_getNext(itor)) != NULL) {
        stList *contigEntries = stHash_search(entries, contigName);
        stList_sort(contigEntries, phasedVariant_positionCmp);
        assert(((PhasedVariant*) stList_get(contigEntries, 0))->refPos <= ((PhasedVariant*) stList_get(contigEntries, stList_length(contigEntries) - 1))->refPos);
    }
    stHash_destructIterator(itor);

    return entries;
}

stList *getSharedContigs(stHash *entry1, stHash *entry2) {

    // get contigs
    stSet *contigs1 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stSet *contigs2 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stList *keys = stHash_getKeys(entry1);
    for (int64_t i = 0; i < stList_length(keys); i++) {
        stSet_insert(contigs1, stString_copy(stList_get(keys, i)));
    }
    stList_destruct(keys);
    keys = stHash_getKeys(entry2);
    for (int64_t i = 0; i < stList_length(keys); i++) {
        stSet_insert(contigs2, stString_copy(stList_get(keys, i)));
    }
    stList_destruct(keys);

    // save intersection
    stList *sharedContigs = stList_construct3(0, free);
    stSet *intersection = stSet_getIntersection(contigs1, contigs2);
    stSetIterator *itor = stSet_getIterator(intersection);
    char *key = NULL;
    while ((key = stSet_getNext(itor)) != NULL) {
        stList_append(sharedContigs, stString_copy(key));
    }
    stSet_destructIterator(itor);

    // sort
    stList_sort(sharedContigs, (int (*)(const void *, const void *)) strcmp);

    // cleanup
    stSet_destruct(intersection);
    stSet_destruct(contigs1);
    stSet_destruct(contigs2);

    return sharedContigs;
}

PartialPhaseSums* partialPhaseSums_construct(const char *queryPhaseSet, const char *truthPhaseSet) {
    PartialPhaseSums* pps = (PartialPhaseSums*) malloc(sizeof(PartialPhaseSums));
    pps->queryPhaseSet = stString_copy(queryPhaseSet);
    pps->truthPhaseSet = stString_copy(truthPhaseSet);
    pps->unphasedSum = 0.0;
    pps->phaseSum1 = 0.0;
    pps->phaseSum2 = 0.0;
    return pps;
}

void partialPhaseSums_destruct(PartialPhaseSums *pps) {
    free(pps->queryPhaseSet);
    free(pps->truthPhaseSet);
    free(pps);
}

stHash *phaseSetIntervals(stList *phasedVariants) {
    
    stHash *intervals = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    int64_t prevPos = -1;
    for (int64_t i = 0; i < stList_length(phasedVariants); ++i) {
        PhasedVariant *pv = stList_get(phasedVariants, i);
        if (prevPos > pv->refPos) {
            st_errAbort("error: phased variant at position %lli on sequence %s is out of order with position %lli\n", pv->refPos, pv->refSeqName, prevPos);
        }
        prevPos = pv->refPos;
        int64_t *interval = stHash_search(intervals, pv->phaseSet);
        if (interval == NULL) {
            interval = (int64_t*) malloc(2 * sizeof(int64_t));
            interval[0] = i;
            stHash_insert(intervals, stString_copy(pv->phaseSet), interval);
            
        }
        interval[1] = i;
    }
    return intervals;
}

double *phasingCorrectnessInternal(stList *queryPhasedVariants, stList *truthPhasedVariants, double decay,
                                   stHash *queryPhaseSetIntervals, stHash *truthPhaseSetIntervals, bool forward) {
    
    // TODO: what's the most sensible way to handle het variants that only occur in one VCF?
    // for now, just skipping them
    
    // holds the partial sums for each active pair of phase sets
    stList *phaseSetPartialSums = stList_construct3(0, (void (*)(void *)) partialPhaseSums_destruct);
    
    // accumulator for the sum
    double totalSum = 0.0;
    
    // accumulator for the max value of the sum
    double partitionSum = 0.0;
    double partitionTotalSum = 0.0;
    
    // accumulator for the unphased partial sums of phase set pairs that have fallen out of scope
    double outOfScopeSum = 0.0;
    
    // which direction are we iterating down the list of variants
    int64_t i, j, incr;
    if (forward) {
        i = 0;
        j = 0;
        incr = 1;
    }
    else {
        i = stList_length(queryPhasedVariants) - 1;
        j = stList_length(truthPhasedVariants) - 1;
        incr = -1;
    }
    
    while (i >= 0 && i < stList_length(queryPhasedVariants) && j >= 0 && j < stList_length(truthPhasedVariants)) {
        PhasedVariant *qpv = stList_get(queryPhasedVariants, i);
        PhasedVariant *tpv = stList_get(truthPhasedVariants, j);
        
        if ((qpv->refPos < tpv->refPos && forward) ||
            (qpv->refPos > tpv->refPos && !forward)) {
            // variant only in query
            i += incr;
        }
        else if ((tpv->refPos < qpv->refPos && forward) ||
                 (tpv->refPos > qpv->refPos && !forward)) {
            // variant only in truth
            j += incr;
        }
        else {
            
            // match up the alleles
            bool match11 = (strcmp(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt1)) == 0);
            bool match12 = (strcmp(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt2)) == 0);
            bool match21 = (strcmp(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt1)) == 0);
            bool match22 = (strcmp(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt2)) == 0);
            
            if (!(match11 || match12) || !(match21 || match22)) {
                // the site is shared, but the alleles are not, just skip this variant
                // TODO: is this the best way to handle this case?
                continue;
            }
            
            if ((int) match11 + (int) match12 + (int) match21 + (int) match22 > 2) {
                // at least one allele must be duplicated in the list of alts
                st_errAbort("error: duplicate alleles detected at position %lli on sequence %s\n",
                        qpv->refPos, qpv->refSeqName);
            }
            
            // do we find a phase set pair that matches this variant's phase set pair?
            bool foundCophasedSum = false;
            
            // add the contribution to the total of the appropriate partial sum,
            // and register a correctly phased variant pair wherever necessary
            for (int64_t k = 0; k < stList_length(phaseSetPartialSums); ++k) {
                
                PartialPhaseSums *sums = stList_get(phaseSetPartialSums, k);

                if (strcmp(qpv->phaseSet, sums->queryPhaseSet) == 0 &&
                    strcmp(tpv->phaseSet, sums->truthPhaseSet) == 0) {
                    // the current pair of variants are co-phased with the variants
                    // that make up this partial sum
                    
                    foundCophasedSum = true;
                    
                    // because we've filtered down to 1) only het sites, and  2) sites
                    // where the alleles match, the only two combinations of matching
                    // that are allowed are 1-1/2-2 or 1-2/2-1
                    if (match11) {
                        totalSum += sums->phaseSum1;
                        sums->phaseSum1 += 1.0;
                    }
                    else {
                        totalSum += sums->phaseSum2;
                        sums->phaseSum2 += 1.0;
                    }
                }
                else {
                    // this is a different phase set
                    totalSum += sums->unphasedSum;
                }
                // the unphased sum acts as if always correctly phased
                sums->unphasedSum += 1.0;
            }
            totalSum += outOfScopeSum;
            
            // partition function is the max value, always counts pairs as phased
            partitionTotalSum += partitionSum;
            partitionSum += 1.0;
            
            if (!foundCophasedSum) {
                // this is the first time we've found this phase set pair, we need
                // to initialize a new partial sum for it
                stList_append(phaseSetPartialSums,
                              partialPhaseSums_construct(qpv->phaseSet, tpv->phaseSet));
            }
            
            // decay all of the partial sums to prepare for the next iteration
            for (int64_t k = 0; k < stList_length(phaseSetPartialSums); ++k) {
                PartialPhaseSums *sums = stList_get(phaseSetPartialSums, k);
                sums->unphasedSum *= decay;
                sums->phaseSum1 *= decay;
                sums->phaseSum2 *= decay;
            }
            partitionSum *= decay;
            outOfScopeSum *= decay;
            
            i += incr;
            j += incr;
        }
        
        // check if any phase set pairs have fallen out of scope
        for (int64_t k = 0; k < stList_length(phaseSetPartialSums);) {
            
            PartialPhaseSums *sums = stList_get(phaseSetPartialSums, k);
            int64_t *queryInterval = stHash_search(queryPhaseSetIntervals, sums->queryPhaseSet);
            int64_t *truthInterval = stHash_search(truthPhaseSetIntervals, sums->truthPhaseSet);
            
            if (i < queryInterval[0] || i > queryInterval[1]
                || j < truthInterval[0] || j > truthInterval[1]) {
                // one of the phase sets has fallen out of scope, the unphased summands in this phase set
                // pair can now be accumulated in the out of scope accumulator
                
                outOfScopeSum += sums->unphasedSum;
                
                // remove this sum from the list of partial sums
                partialPhaseSums_destruct(sums);
                stList_set(phaseSetPartialSums, k, stList_get(phaseSetPartialSums, stList_length(phaseSetPartialSums) - 1));
                stList_pop(phaseSetPartialSums);
            }
            else {
                ++k;
            }
        }
    }
    
    stList_destruct(phaseSetPartialSums);
    
    // package the two values and return them
    double *returnVal = malloc(2 * sizeof(double));
    returnVal[0] = totalSum;
    returnVal[1] = partitionTotalSum;
    return returnVal;
}

double switchCorrectness(stList *queryPhasedVariants, stList *truthPhasedVariants) {
     
    // TODO: finish this
    
    int64_t prev_i = -1;
    int64_t prev_j = -1;
    for (int64_t i = 0, j = 0; i < stList_length(queryPhasedVariants) && j < stList_length(truthPhasedVariants);) {
        PhasedVariant *qpv = stList_get(queryPhasedVariants, i);
        PhasedVariant *tpv = stList_get(truthPhasedVariants, j);
        if (qpv->refPos < tpv->refPos) {
            // variant only in query
            ++i;
        }
        else if (tpv->refPos < qpv->refPos) {
            // variant only in truth
            ++j;
        }
        else {
            
            ++i;
            ++j;
        }
    }
    
    return 0.0;
}

double phasingCorrectness(stList *queryPhasedVariants, stList *truthPhasedVariants, double decay) {
    
    if (decay < 0.0 || decay > 1.0) {
        st_errAbort("error: decay factor is %d, must be between 0.0 and 1.0\n", decay);
    }
    
    if (decay == 0.0) {
        // this has to be handled as a special case, because it's actually a limit rather
        // than direct evaluation. if computed directly, leads to division by 0
        
        return switchCorrectness(queryPhasedVariants, truthPhasedVariants);
    }
    
    // the interval of variant indexes that each phase set is contained in
    stHash *queryPhaseSetIntervals = phaseSetIntervals(queryPhasedVariants);
    stHash *truthPhaseSetIntervals = phaseSetIntervals(truthPhasedVariants);
    
    double *forwardSums = phasingCorrectnessInternal(queryPhasedVariants, truthPhasedVariants, decay,
                                                     queryPhaseSetIntervals, truthPhaseSetIntervals, true);
    double *reverseSums = phasingCorrectnessInternal(queryPhasedVariants, truthPhasedVariants, decay,
                                                     queryPhaseSetIntervals, truthPhaseSetIntervals, false);
    
    
    stHash_destruct(queryPhaseSetIntervals);
    stHash_destruct(truthPhaseSetIntervals);
    
    double correctness = (forwardSums[0] + reverseSums[0]) / (forwardSums[1] + reverseSums[1]);
    free(forwardSums);
    free(reverseSums);
    
    return correctness;
}

void usage() {
    fprintf(stderr, "usage: localPhasingCorrectness <TRUTH_VCF> <QUERY_VCF> \n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Generate LPC data for phase sets in both VCFs.\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {

    if (argc < 3) {
        usage();
        return 0;
    }

    char *truthVcfFile = stString_copy(argv[1]);
    char *queryVcfFile = stString_copy(argv[2]);

    // for logging
    st_setLogLevel(info);

    // sanity check (verify files are accessible)
    if (access(truthVcfFile, R_OK) != 0) {
        st_errAbort("Could not read from truth vcf file: %s\n", truthVcfFile);
    }
    if (access(queryVcfFile, R_OK) != 0) {
        st_errAbort("Could not read from query vcf file: %s\n", queryVcfFile);
    }

    stHash *truthVariants = getPhasedVariants(truthVcfFile);
    stHash *queryVariants = getPhasedVariants(queryVcfFile);
    
    free(truthVcfFile);
    free(queryVcfFile);

    stList *sharedContigs = getSharedContigs(truthVariants, queryVariants);
    st_logCritical("Got %"PRId64" shared contigs (truth %"PRId64", query %"PRId64")\n", stList_length(sharedContigs),
            stHash_size(truthVariants), stHash_size(queryVariants));
    
    // TODO: magic number, make configurable
    int64_t numLengthScales = 20;
    
    
    assert(numLengthScales >= 2);
    double *decayValues = (double*) malloc(numLengthScales * sizeof(double));
    double *lengthScales = (double*) malloc(numLengthScales * sizeof(double));
    
//    // measure the longest contig, so we have an idea of what scale we need to represent
//    double longestContigLength = 0.0;
//    for (int64_t i = 0; i < stList_length(sharedContigs); ++i) {
//        stList *contigTruthVariants = stHash_search(truthVariants, stList_get(sharedContigs, i));
//        stList *contigQueryVariants = stHash_search(queryVariants, stList_get(sharedContigs, i));
//        longestContigLength = fmax(longestContigLength, stList_length(contigTruthVariants));
//        longestContigLength = fmax(longestContigLength, stList_length(contigQueryVariants));
//    }
//
//    double logLongest = log(longestContigLength);
//    for (int64_t i = 1; i + 1 < numLengthScales; ++i) {
//        double x = logLongest - i * logLongest / (numLengthScales - 1);
//        double a = exp(x);
//        decayValues[i] = 1.0 - exp(-log(2.0) / a);
//
//    }
    
    // an even grid for now, just to have something up and running
    double denom = numLengthScales - 1;
    for (int64_t i = 0; i < numLengthScales; ++i) {
        decayValues[i] = i / denom;
        lengthScales[i] = 1.0 - log(2.0) / log(decayValues[i]);
    }
    
    printf("decay\tlength_scale");
    for (int64_t i = 0; i < stList_length(sharedContigs); ++i) {
        printf("\t%s", (char*) stList_get(sharedContigs, i));
    }
    printf("\n");
    for (int64_t i = 0; i < numLengthScales; ++i) {
        printf("%f\t%f", decayValues[i], lengthScales[i]);
        for (int64_t j = 0; j < stList_length(sharedContigs); ++j) {
            stList *contigTruthVariants = stHash_search(truthVariants, stList_get(sharedContigs, j));
            stList *contigQueryVariants = stHash_search(queryVariants, stList_get(sharedContigs, j));
            
            double correctness = phasingCorrectness(contigTruthVariants, contigQueryVariants, decayValues[i]);
            
            printf("\t%f", correctness);
        }
        printf("\n");
    }
        
    stList_destruct(sharedContigs);
    stHash_destruct(truthVariants);
    stHash_destruct(queryVariants);
    
    return 0;
}

