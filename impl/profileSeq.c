/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"

/*
 * Functions for profile sequence
 */

stProfileSeq *stProfileSeq_constructEmptyProfile(stReference *ref, char *readId,
                                                 int64_t referenceStart, int64_t length) {
    /*
     * Creates an empty profile sequence, with all the profile probabilities set to 0.
     */

    stProfileSeq *seq = st_calloc(1, sizeof(stProfileSeq));
    seq->ref = ref;
    seq->readId = stString_copy(readId);
    seq->refStart = referenceStart;
    seq->length = length;
    seq->alleleOffset = ref->sites[referenceStart].alleleOffset;
    uint64_t lastAllele = referenceStart + length < ref->length ? ref->sites[referenceStart + length].alleleOffset
                                                                : ref->totalAlleles;
    seq->profileProbs = st_calloc(lastAllele - seq->alleleOffset, sizeof(uint8_t));
    return seq;
}

void stProfileSeq_destruct(stProfileSeq *seq) {
    /*
     * Cleans up memory for profile sequence.
     */

    free(seq->profileProbs);
    free(seq->readId);
    free(seq);
}

uint8_t *stProfileSeq_getProb(stProfileSeq *seq, uint64_t site, uint64_t allele) {
    /*
     * Gets probability of a given character as a float.
     */
    return &(seq->profileProbs[seq->ref->sites[site].alleleOffset - seq->alleleOffset +
                               allele]); //((float)p[characterIndex])/255;
}

void stProfileSeq_print(stProfileSeq *seq, FILE *fileHandle) {
    /*
     * Prints a debug representation of a profile sequence.
     */

    char profileString[seq->length + 1];
    profileString[seq->length] = '\0';
    for (int64_t i = 0; i < seq->length; i++) {
        stSite *site = &(seq->ref->sites[i + seq->refStart]);
        uint8_t maxProb = seq->profileProbs[site->alleleOffset - seq->alleleOffset];
        int64_t maxAllele = 0;
        for (int64_t j = 1; j < site->alleleNumber; j++) {
            uint8_t prob = seq->profileProbs[site->alleleOffset - seq->alleleOffset + j];
            if (prob < maxProb) {
                maxProb = prob;
                maxAllele = j;
            }
        }
        //TODO this prints weird boxes 'cause this isn't a char
        profileString[i] = maxAllele;
    }

    fprintf(fileHandle, "\tSEQUENCE REF_NAME: %s REF_START %"
                        PRIi64 " REF_LENGTH: %" PRIi64 " ML_STRING: %s\n",
            seq->ref->referenceName, seq->refStart, seq->length, profileString);
}

static int cmpint64(int64_t i, int64_t j) {
    return i > j ? 1 : i < j ? -1 : 0;
}
