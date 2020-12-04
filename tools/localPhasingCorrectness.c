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

stHash *getPhasedVariants(char *vcfFile) {
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

    stList *sharedContigs = getSharedContigs(truthVariants, queryVariants);
    st_logCritical("Got %"PRId64" shared contigs (truth %"PRId64", query %"PRId64")\n", stList_length(sharedContigs),
            stHash_size(truthVariants), stHash_size(queryVariants));


    return 0;
}

