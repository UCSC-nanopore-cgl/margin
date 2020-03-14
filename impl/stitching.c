//
// Created by Benedict Paten on 3/14/20.
//
// Code for stitching together "chunks" of inferred sequence
//

#include "margin.h"

int64_t removeOverlap(char *prefixString, char *suffixString, int64_t approxOverlap, PolishParams *polishParams,
                      int64_t *prefixStringCropEnd, int64_t *suffixStringCropStart) {

    // Align the overlapping suffix of the prefixString and the prefix of the suffix string
    int64_t prefixStringLength = strlen(prefixString);
    int64_t suffixStringLength = strlen(suffixString);

    // Get coordinates of substrings to be aligned
    int64_t i = (prefixStringLength - approxOverlap) < 0 ? 0 : prefixStringLength - approxOverlap;
    int64_t j = approxOverlap < suffixStringLength ? approxOverlap : suffixStringLength;

    // Crop suffix
    char c = suffixString[j];
    suffixString[j] = '\0';

    // Symbol strings
    SymbolString sX = symbolString_construct(&(prefixString[i]), 0, strlen(&(prefixString[i])), polishParams->alphabet);
    SymbolString sY = symbolString_construct(suffixString, 0, strlen(suffixString), polishParams->alphabet);


    // Use default state machine for alignment
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);

    // Run the alignment
    stList *alignedPairs = getAlignedPairs(sM, sX, sY, polishParams->p, 1, 1);

    //
    symbolString_destruct(sX);
    symbolString_destruct(sY);
    stateMachine_destruct(sM);

    /*for(uint64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *aPair = stList_get(alignedPairs, i);
        st_uglyf("Boo %i %i %i\n", (int)stIntTuple_get(aPair, 0), (int)stIntTuple_get(aPair, 1),
                (int)stIntTuple_get(aPair, 2));
    }*/

    if(stList_length(alignedPairs) == 0 && st_getLogLevel() >= info) {
        st_logInfo("    Failed to find good overlap. Suffix-string: %s\n", &(prefixString[i]));
        st_logInfo("    Failed to find good overlap. Prefix-string: %s\n", suffixString);
    }

    // Remove the suffix crop
    suffixString[j] = c;

    // Pick the median point
    stIntTuple *maxPair = NULL;
    for(int64_t k=0; k<stList_length(alignedPairs); k++) {
        stIntTuple *aPair = stList_get(alignedPairs, k);
        if(maxPair == NULL || stIntTuple_get(aPair, 0) > stIntTuple_get(maxPair, 0)) {
            maxPair = aPair;
        }
    }
    if(maxPair == NULL) {
        st_logCritical("    Failed to find any aligned pairs between overlapping strings, not "
                       "doing any trimming (approx overlap: %i, len x: %i, len y: %i)\n", approxOverlap, prefixStringLength, suffixStringLength);
        *prefixStringCropEnd = prefixStringLength;
        *suffixStringCropStart = 0;
    }
    else {
        *prefixStringCropEnd = stIntTuple_get(maxPair, 1) + i; // Exclusive
        *suffixStringCropStart = stIntTuple_get(maxPair, 2);  // Inclusive
    }

    int64_t overlapWeight = maxPair == NULL ? 0 : stIntTuple_get(maxPair, 0);

    stList_destruct(alignedPairs);

    return overlapWeight;
}