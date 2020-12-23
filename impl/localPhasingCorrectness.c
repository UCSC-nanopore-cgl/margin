/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "localPhasingCorrectness.h"

PhasedVariant *phasedVariant_construct(const char *refSeqName, int64_t refPos, double quality, stList *alleles, int64_t gt1, int64_t gt2, const char * phaseSet) {
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
        st_errAbort("error: Could not open VCF %s\n", vcfFile);
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
        st_errAbort("error: PS tag not present in VCF header for %s", vcfFile);
    }
    int psType = bcf_hdr_id2type(hdr, BCF_HL_FMT, psId);
    if (psType == BCF_HT_INT) {
        phaseSetIsInt = TRUE;
    } else if (psType == BCF_HT_STR) {
        phaseSetIsInt = FALSE;
    } else {
        st_errAbort("error: Unknown PS type in VCF header for %s", vcfFile);
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
    st_logCritical("Read %"PRId64" variants from %s over %"PRId64" contigs in %"PRId64"s, keeping %"PRId64" phased variants"
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

VariantCorrectness *variantCorrectness_construct(int64_t refPos, double correctness, double maxCorrectness) {
    VariantCorrectness* vc = malloc(sizeof(VariantCorrectness));
    vc->refPos = refPos;
    vc->correctness = correctness;
    vc->maxCorrectness = maxCorrectness;
    return vc;
}

void variantCorrectness_destruct(VariantCorrectness* vc) {
    free(vc);
}

double meanVariantDist(stHash *query, stHash *truth, stList *sharedContigs) {
    
    int64_t distSum = 0;
    int64_t numPairs = 0;
    for (int64_t i = 0; i < stList_length(sharedContigs); ++i) {
        char *contig = stList_get(sharedContigs, i);
        stList *queryPhasedVariants = stHash_search(query, contig);
        stList *truthPhasedVariants = stHash_search(truth, contig);
        int64_t prevPos = -1;
        for (int64_t j = 0, k = 0; j < stList_length(queryPhasedVariants) && k < stList_length(truthPhasedVariants); ) {
            
            PhasedVariant *qpv = stList_get(queryPhasedVariants, j);
            PhasedVariant *tpv = stList_get(truthPhasedVariants, k);
            
            if (qpv->refPos < tpv->refPos) {
                // variant only in query
                ++j;
            }
            else if (tpv->refPos < qpv->refPos) {
                // variant only in truth
                ++k;
            }
            else {
                
                // TODO: duplicative with localCorrectnessInternal
                
                // match up the alleles
                bool match11 = stString_eq(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt1));
                bool match12 = stString_eq(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt2));
                bool match21 = stString_eq(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt1));
                bool match22 = stString_eq(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt2));
                
                ++j;
                ++k;
                
                if (!(match11 || match12) || !(match21 || match22)) {
                    // the site is shared, but the alleles are not, just skip this variant
                    // TODO: is this the best way to handle this case?
                    continue;
                }
                
                if ((int) match11 + (int) match12 + (int) match21 + (int) match22 > 2) {
                    // at least one allele must be duplicated in the list of alts
                    st_logCritical("error: duplicate alleles detected at position %"PRId64" on sequence %s\n",
                                   qpv->refPos, qpv->refSeqName);
                    continue;
                }
                
                if (prevPos != -1) {
                    distSum += (qpv->refPos - prevPos);
                    ++numPairs;
                }
                
                prevPos = qpv->refPos;
            }
        }
    }
    
    return ((double) distSum) / numPairs;
}

PartialPhaseSums* partialPhaseSums_construct(const char *queryPhaseSet, const char *truthPhaseSet) {
    PartialPhaseSums* pps = (PartialPhaseSums*) malloc(sizeof(PartialPhaseSums));
    pps->queryPhaseSet = stString_copy(queryPhaseSet);
    pps->truthPhaseSet = stString_copy(truthPhaseSet);
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
            st_errAbort("error: Phased variant at position %"PRId64" on sequence %s is out of order with position %"PRId64"\n", pv->refPos, pv->refSeqName, prevPos);
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

double *phasingCorrectnessInternal(stList *queryPhasedVariants, stList *truthPhasedVariants,
                                   double decay, bool bySeqDist, bool crossBlockCorrect,
                                   stHash *queryPhaseSetIntervals, stHash *truthPhaseSetIntervals, bool forward,
                                   stList *variantCorrectnessOut) {
    
    // TODO: what's the most sensible way to handle het variants that only occur in one VCF?
    // for now, just skipping them
    
    // holds the partial sums for each active pair of phase sets
    stList *phaseSetPartialSums = stList_construct3(0, (void (*)(void *)) partialPhaseSums_destruct);
    
    // accumulator for the sum
    double totalSum = 0.0;
    
    // accumulator for the max value of the sum
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
    
    int64_t prevPosition = -1;
    
    st_logDebug("beginning %s sum\n", forward ? "forward" : "backward");
    
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
            bool match11 = stString_eq(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt1));
            bool match12 = stString_eq(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt2));
            bool match21 = stString_eq(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt1));
            bool match22 = stString_eq(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt2));
            
            i += incr;
            j += incr;
            
            if (!(match11 || match12) || !(match21 || match22)) {
                // the site is shared, but the alleles are not, just skip this variant
                // TODO: is this the best way to handle this case?
                continue;
            }
            
            if ((int) match11 + (int) match12 + (int) match21 + (int) match22 > 2) {
                // at least one allele must be duplicated in the list of alts
                continue;
            }
            
            double decayValue;
            if (bySeqDist) {
                decayValue = pow(decay, fabs((double) (qpv->refPos - prevPosition)));
            }
            else {
                decayValue = decay;
            }
            
            // decay all of the previous partial sums
            for (int64_t k = 0; k < stList_length(phaseSetPartialSums); ++k) {
                PartialPhaseSums *sums = stList_get(phaseSetPartialSums, k);
                sums->phaseSum1 *= decayValue;
                sums->phaseSum2 *= decayValue;
            }
            outOfScopeSum *= decayValue;
            
            st_logDebug("going into iteration query %"PRId64", truth %"PRId64":\n\tref pos %"PRId64"\n\titer decay %f\n\ttotal %f\n\tpartition total %f\n\tout of scope sum %f\n", i - incr, j - incr, qpv->refPos, decayValue, totalSum, partitionTotalSum, outOfScopeSum);
            
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
                    
                    // partition functions acts as if correctly phased with everything
                    partitionTotalSum += sums->phaseSum1 + sums->phaseSum2;
                    
                    // because we've filtered down to 1) only het sites, and  2) sites
                    // where the alleles match, the only two combinations of matching
                    // that are allowed are 1-1/2-2 or 1-2/2-1
                    if (match11) {
                        totalSum += sums->phaseSum1;
                        sums->phaseSum1 += 1.0;
                        
                        if (variantCorrectnessOut) {
                            stList_append(variantCorrectnessOut, variantCorrectness_construct(qpv->refPos, sums->phaseSum1,
                                                                                              sums->phaseSum1 + sums->phaseSum2));
                        }
                    }
                    else {
                        totalSum += sums->phaseSum2;
                        sums->phaseSum2 += 1.0;
                        
                        if (variantCorrectnessOut) {
                            stList_append(variantCorrectnessOut, variantCorrectness_construct(qpv->refPos, sums->phaseSum2,
                                                                                              sums->phaseSum1 + sums->phaseSum2));
                        }
                    }
                }
                else if (crossBlockCorrect) {
                    // this is a different phase set, but we're counting those summands
                    // as correct
                    totalSum += sums->phaseSum1 + sums->phaseSum2;
                    partitionTotalSum += sums->phaseSum1 + sums->phaseSum2;
                    
                    if (variantCorrectnessOut) {
                        stList_append(variantCorrectnessOut, variantCorrectness_construct(qpv->refPos, sums->phaseSum1 + sums->phaseSum2,
                                                                                          sums->phaseSum1 + sums->phaseSum2));
                    }
                }
            }
            
            // add any summands from out of scope phase set pairs
            totalSum += outOfScopeSum;
            partitionTotalSum += outOfScopeSum;
            
            if (!foundCophasedSum) {
                // this is the first time we've found this phase set pair, we need
                // to initialize a new partial sum for it
                PartialPhaseSums *sums = partialPhaseSums_construct(qpv->phaseSet, tpv->phaseSet);
                if (match11) {
                    sums->phaseSum1 = 1.0;
                }
                else {
                    sums->phaseSum2 = 1.0;
                }
                stList_append(phaseSetPartialSums, sums);
                
                if (variantCorrectnessOut) {
                    stList_append(variantCorrectnessOut, variantCorrectness_construct(qpv->refPos, 0.0, 0.0));
                }
            }
            
            if (variantCorrectnessOut) {
                VariantCorrectness *vc = stList_get(variantCorrectnessOut, stList_length(variantCorrectnessOut) - 1);
                vc->correctness += outOfScopeSum;
                vc->maxCorrectness += outOfScopeSum;
            }
            
            prevPosition = qpv->refPos;
        }
        
        // check if any phase set pairs have fallen out of scope
        for (int64_t k = 0; k < stList_length(phaseSetPartialSums);) {
            
            PartialPhaseSums *sums = stList_get(phaseSetPartialSums, k);
            int64_t *queryInterval = stHash_search(queryPhaseSetIntervals, sums->queryPhaseSet);
            int64_t *truthInterval = stHash_search(truthPhaseSetIntervals, sums->truthPhaseSet);
            
            st_logDebug("end of iter partial sum %"PRId64":\n\tquery phase set: %s\n\ttruth phase set: %s\n\tphased sum 1: %f\n\tphased sum 2: %f\n", k, sums->queryPhaseSet, sums->truthPhaseSet, sums->phaseSum1, sums->phaseSum2);
            
            if (i < queryInterval[0] || i > queryInterval[1]
                || j < truthInterval[0] || j > truthInterval[1]) {
                // one of the phase sets has fallen out of scope, the unphased summands in this phase set
                // pair can now be accumulated in the out of scope accumulator
                
                st_logDebug("\t\tthis sum falls out of scope at this iteration\n");
                
                if (crossBlockCorrect) {
                    outOfScopeSum += sums->phaseSum1 + sums->phaseSum2;
                }
                
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

double switchCorrectness(stList *queryPhasedVariants, stList *truthPhasedVariants, bool bySeqDist,
                         bool crossBlockCorrect, double *phasedPairsOut, stList *variantCorrectnessOut) {
    
    char *prevQueryPhaseSet = NULL;
    char *prevTruthPhaseSet = NULL;
    bool prevVarInPhase = false;
    int64_t prevPosition = -1;
    int64_t minAdjacentDist = LONG_MAX;
    int64_t numCorrectlyPhasedPairs = 0;
    int64_t numPossiblyPhasedPairs = 0;
    int64_t minCounted = 0;
    bool prevPairCounted = false;
    bool prevPairCorrect = false;
    bool pairCounted = false;
    bool pairCorrect = false;
    
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
            
            // TODO: duplicative with phasingCorrectnessInternal
            
            // match up the alleles
            bool match11 = stString_eq(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt1));
            bool match12 = stString_eq(stList_get(qpv->alleles, qpv->gt1), stList_get(tpv->alleles, tpv->gt2));
            bool match21 = stString_eq(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt1));
            bool match22 = stString_eq(stList_get(qpv->alleles, qpv->gt2), stList_get(tpv->alleles, tpv->gt2));
            
            ++i;
            ++j;
            
            if (!(match11 || match12) || !(match21 || match22)) {
                // the site is shared, but the alleles are not, just skip this variant
                // TODO: is this the best way to handle this case?
                continue;
            }
            
            if ((int) match11 + (int) match12 + (int) match21 + (int) match22 > 2) {
                // at least one allele must be duplicated in the list of alts
                continue;
            }
            
            pairCounted = false;
            pairCorrect = false;
            
            if (prevQueryPhaseSet != NULL && prevTruthPhaseSet != NULL) {
                
                int64_t phasePairDist = qpv->refPos - prevPosition;
                
                st_logDebug("checking phasing between var at position %"PRId64" in phase %d and var at position %"PRId64" in phase %d\n", prevPosition, prevVarInPhase, qpv->refPos, match11);
                
                // TODO: this could break if the minimum distance is achieved at pairs that aren't
                
                bool phaseSetPairMatch = (stString_eq(qpv->phaseSet, prevQueryPhaseSet)
                                          && stString_eq(tpv->phaseSet, prevTruthPhaseSet));
                
                if (phasePairDist < minAdjacentDist && bySeqDist && (phaseSetPairMatch || crossBlockCorrect)) {
                    // only the nearest pairs together count, so we throw away all previously
                    // recorded pairs that had a larger distance between
                    // note: the min distance will not be maintained unless using by seq dist
                    
                    st_logDebug("new minimum distance of %"PRId64" from %"PRId64", discarding previous pair phasing data\n", phasePairDist, minAdjacentDist);
                    
                    numPossiblyPhasedPairs = 0;
                    numCorrectlyPhasedPairs = 0;
                    prevPairCounted = false;
                    minAdjacentDist = phasePairDist;
                    
                    if (variantCorrectnessOut) {
                        minCounted = stList_length(variantCorrectnessOut);
                    }
                }
                
                if (phasePairDist == minAdjacentDist || !bySeqDist) {
                    
                    pairCounted = (phaseSetPairMatch || crossBlockCorrect);
                    pairCorrect = ((phaseSetPairMatch && match11 == prevVarInPhase) ||
                                   (!phaseSetPairMatch && crossBlockCorrect));
                    
                    if (pairCounted) {
                        ++numPossiblyPhasedPairs;
                    }
                    if (pairCorrect) {
                        ++numCorrectlyPhasedPairs;
                    }
                    
                    st_logDebug("num phased incremented to %"PRId64", possibly phased incremented to %"PRId64"\n", numCorrectlyPhasedPairs, numPossiblyPhasedPairs);
                }
            }
            
            if (variantCorrectnessOut) {
                stList_append(variantCorrectnessOut, variantCorrectness_construct(qpv->refPos, 0.0, 0.0));
                if (stList_length(variantCorrectnessOut) > 1) {
                    VariantCorrectness* pvc = stList_get(variantCorrectnessOut,
                                                         stList_length(variantCorrectnessOut) - 2);
                    
                    pvc->correctness = ((int) (prevPairCorrect && prevPairCounted)) + ((int) (pairCorrect && pairCounted));
                    pvc->maxCorrectness = ((int) prevPairCounted + (int) pairCounted);
                }
            }
            
            prevVarInPhase = match11;
            prevQueryPhaseSet = qpv->phaseSet;
            prevTruthPhaseSet = tpv->phaseSet;
            prevPosition = qpv->refPos;
            prevPairCorrect = pairCorrect;
            prevPairCounted = pairCounted;
        }
    }
    
    if (variantCorrectnessOut && stList_length(variantCorrectnessOut) != 0) {
        VariantCorrectness* vc = stList_get(variantCorrectnessOut,
                                            stList_length(variantCorrectnessOut) - 1);
        vc->correctness = pairCorrect && pairCounted;
        vc->maxCorrectness = pairCounted;
        
        // reset any variants that were counted before finding the min distance
        // (if doing sequence-based distance)
        for (int64_t i = 0; i < minCounted; ++i) {
            VariantCorrectness* uncounted = stList_get(variantCorrectnessOut, i);
            uncounted->correctness = 0.0;
            uncounted->maxCorrectness = 0.0;
        }
    }
    
    if (phasedPairsOut != NULL) {
        *phasedPairsOut = numPossiblyPhasedPairs;
    }
    
    return ((double) numCorrectlyPhasedPairs) / numPossiblyPhasedPairs;
}

double phasingCorrectness(stList *queryPhasedVariants, stList *truthPhasedVariants, double decay,
                          bool bySeqDist, bool crossBlockCorrect, double *effectivePairCountOut,
                          stList* variantCorrectnessOut) {
    
    if (decay < 0.0 || decay > 1.0) {
        st_errAbort("error: Decay factor is %d, must be between 0.0 and 1.0\n", decay);
    }
    
    st_logDebug("calculating correctness with decay %f\n", decay);
    
    if (decay == 0.0) {
        // this has to be handled as a special case, because it's actually a limit rather
        // than direct evaluation. if computed directly, leads to division by 0
        
        double correctness = switchCorrectness(queryPhasedVariants, truthPhasedVariants, bySeqDist, crossBlockCorrect,
                                               effectivePairCountOut, variantCorrectnessOut);
        
        return correctness;
    }
    
    // the interval of variant indexes that each phase set is contained in
    stHash *queryPhaseSetIntervals = phaseSetIntervals(queryPhasedVariants);
    stHash *truthPhaseSetIntervals = phaseSetIntervals(truthPhasedVariants);
    
    stList *revVariantCorrectness = NULL;
    if (variantCorrectnessOut) {
        revVariantCorrectness = stList_construct3(0, (void (*)(void*)) variantCorrectness_destruct);
    }
    
    double *forwardSums = phasingCorrectnessInternal(queryPhasedVariants, truthPhasedVariants, decay,
                                                     bySeqDist, crossBlockCorrect,
                                                     queryPhaseSetIntervals, truthPhaseSetIntervals,
                                                     true, variantCorrectnessOut);
    double *reverseSums = phasingCorrectnessInternal(queryPhasedVariants, truthPhasedVariants, decay,
                                                     bySeqDist, crossBlockCorrect,
                                                     queryPhaseSetIntervals, truthPhaseSetIntervals,
                                                     false, revVariantCorrectness);
    
    // add the backward half of the variants' sum into the forward half and normalize
    if (variantCorrectnessOut) {
        for (int64_t i = 0; i < stList_length(variantCorrectnessOut); ++i) {
            VariantCorrectness *fvc = stList_get(variantCorrectnessOut, i);
            VariantCorrectness *rvc = stList_get(revVariantCorrectness,
                                                 stList_length(revVariantCorrectness) - i - 1);
            fvc->correctness += rvc->correctness;
            fvc->maxCorrectness += rvc->maxCorrectness;
        }
        stList_destruct(revVariantCorrectness);
    }
    
    stHash_destruct(queryPhaseSetIntervals);
    stHash_destruct(truthPhaseSetIntervals);
    
    double correctness = (forwardSums[0] + reverseSums[0]) / (forwardSums[1] + reverseSums[1]);
    if (effectivePairCountOut != NULL) {
        *effectivePairCountOut = forwardSums[1] + reverseSums[1];
    }
    
    st_logDebug("fwd numer %f, bwd numer %f, fwd denom %f, bwd denom %f, final answer %f\n",
                forwardSums[0], reverseSums[0], forwardSums[1], reverseSums[1], correctness);
    
    free(forwardSums);
    free(reverseSums);
    
    return correctness;
}
