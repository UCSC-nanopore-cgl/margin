//
// Created by Benedict Paten on 3/14/20.
//

#include "margin.h"

/*
 * Functions for modeling repeat counts
 */

inline double *
repeatSubMatrix_setLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, int64_t observedRepeatCount,
                           int64_t underlyingRepeatCount) {
    // santiy check
    // TODO fix this!! filter reads before this point? rely on a prior (GC AT N)?
    if (base >= repeatSubMatrix->alphabet->alphabetSize - 1) {
        char *logIdentifier = getLogIdentifier();
        if (base == repeatSubMatrix->alphabet->alphabetSize - 1) {
            st_logInfo(" %s [repeatSubMatrix_setLogProb] base 'Nn' (%d) not supported for repeat estimation! "
                       "Setting to 'A' (0)\n", logIdentifier, base);
        } else {
            st_errAbort(" %s [repeatSubMatrix_setLogProb] base > 'Nn' (%d) not supported for repeat estimation!\n",
                    logIdentifier, base);
        }
        free(logIdentifier);
        base = 0;
    }
    int64_t idx =
            (strand ? base : 3 - base) * repeatSubMatrix->maximumRepeatLength * repeatSubMatrix->maximumRepeatLength +
            underlyingRepeatCount * repeatSubMatrix->maximumRepeatLength +
            observedRepeatCount;
    assert(idx < repeatSubMatrix->maxEntry);
    assert(idx >= 0);
    return &(repeatSubMatrix->logProbabilities[idx]);
}

inline double
repeatSubMatrix_getLogProb(RepeatSubMatrix *repeatSubMatrix, Symbol base, bool strand, int64_t observedRepeatCount,
                           int64_t underlyingRepeatCount) {
    double *loc = repeatSubMatrix_setLogProb(repeatSubMatrix, base, strand, observedRepeatCount, underlyingRepeatCount);
//	printf("%p\n", loc);
    return *loc;
}

void repeatSubMatrix_destruct(RepeatSubMatrix *repeatSubMatrix) {
    alphabet_destruct(repeatSubMatrix->alphabet);
    free(repeatSubMatrix->baseLogProbs_AT);
    free(repeatSubMatrix->baseLogProbs_GC);
    free(repeatSubMatrix->logProbabilities);
    free(repeatSubMatrix);
}

RepeatSubMatrix *repeatSubMatrix_constructEmpty(Alphabet *alphabet) {
    RepeatSubMatrix *repeatSubMatrix = st_calloc(1, sizeof(RepeatSubMatrix));
    repeatSubMatrix->alphabet = alphabet;
    repeatSubMatrix->maximumRepeatLength = MAXIMUM_REPEAT_LENGTH;
    repeatSubMatrix->baseLogProbs_AT = st_calloc(repeatSubMatrix->maximumRepeatLength, sizeof(double));
    repeatSubMatrix->baseLogProbs_GC = st_calloc(repeatSubMatrix->maximumRepeatLength, sizeof(double));
    repeatSubMatrix->maxEntry = repeatSubMatrix->alphabet->alphabetSize * repeatSubMatrix->maximumRepeatLength *
                                repeatSubMatrix->maximumRepeatLength;
    repeatSubMatrix->logProbabilities = st_calloc(repeatSubMatrix->maxEntry, sizeof(double));
    return repeatSubMatrix;
}

double
repeatSubMatrix_getLogProbForGivenRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
                                              stList *bamChunkReads, int64_t underlyingRepeatCount) {
    assert(underlyingRepeatCount < repeatSubMatrix->maximumRepeatLength);
    double logProb = LOG_ONE;
    for (int64_t i = 0; i < stList_length(observations); i++) {
        PoaBaseObservation *observation = stList_get(observations, i);
        BamChunkRead *read = stList_get(bamChunkReads, observation->readNo);
        int64_t observedRepeatCount = read->rleRead->repeatCounts[observation->offset];

        // Be robust to over-long repeat count observations
        observedRepeatCount = observedRepeatCount >= repeatSubMatrix->maximumRepeatLength ?
                              repeatSubMatrix->maximumRepeatLength - 1 : observedRepeatCount;

        logProb += repeatSubMatrix_getLogProb(repeatSubMatrix, base, read->forwardStrand,
                                              observedRepeatCount, underlyingRepeatCount) * observation->weight;
    }

    return logProb / PAIR_ALIGNMENT_PROB_1;
}

void repeatSubMatrix_getMinAndMaxRepeatCountObservations(RepeatSubMatrix *repeatSubMatrix, stList *observations,
                                                         stList *bamChunkReads, int64_t *minRepeatLength,
                                                         int64_t *maxRepeatLength) {
    // Get the range or repeat observations, used to avoid calculating all repeat lengths, heuristically
    *minRepeatLength = repeatSubMatrix->maximumRepeatLength;
    *maxRepeatLength = 0;
    char *maxRLReadId = NULL;
    for (int64_t i = 0; i < stList_length(observations); i++) {
        PoaBaseObservation *observation = stList_get(observations, i);
        BamChunkRead *read = stList_get(bamChunkReads, observation->readNo);
        int64_t observedRepeatCount = read->rleRead->repeatCounts[observation->offset];
        if (observedRepeatCount < *minRepeatLength) {
            *minRepeatLength = observedRepeatCount;
        }
        if (observedRepeatCount > *maxRepeatLength) {
            *maxRepeatLength = observedRepeatCount;
            maxRLReadId = read->readName;
        }
    }
    if (*maxRepeatLength >= repeatSubMatrix->maximumRepeatLength) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(
                " %s Got overlong repeat observation(s), max: %" PRIi64 " (%s), ignoring this and cutting off overlong repeat counts to max\n",
                logIdentifier, *maxRepeatLength, maxRLReadId);
        *maxRepeatLength = repeatSubMatrix->maximumRepeatLength - 1;
        free(logIdentifier);
    }
}

void repeatSubMatrix_getRepeatCountProbs(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
                                         stList *bamChunkReads, double *logProbabilities, int64_t minRepeatLength,
                                         int64_t maxRepeatLength) {
    for (int64_t i = minRepeatLength; i <= maxRepeatLength; i++) {
        logProbabilities[i - minRepeatLength] =
                repeatSubMatrix_getLogProbForGivenRepeatCount(repeatSubMatrix, base, observations, bamChunkReads, i);
    }
}

int64_t repeatSubMatrix_getMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, Symbol base, stList *observations,
                                         stList *bamChunkReads, double *logProbability) {
    int64_t minRepeatLength, maxRepeatLength;
    double logProbabilities[repeatSubMatrix->maximumRepeatLength];

    // Calculate range of repeat counts observed
    repeatSubMatrix_getMinAndMaxRepeatCountObservations(repeatSubMatrix, observations,
                                                        bamChunkReads, &minRepeatLength, &maxRepeatLength);

    if (minRepeatLength == repeatSubMatrix->maximumRepeatLength) {
        *logProbability = LOG_ZERO;
        return 0; // Case we have no valid observations, so assume repeat length of 0
    }

    // Get the prob for each repeat count
    repeatSubMatrix_getRepeatCountProbs(repeatSubMatrix, base, observations,
                                        bamChunkReads, logProbabilities, minRepeatLength, maxRepeatLength);

    return getMax(logProbabilities, maxRepeatLength - minRepeatLength + 1, logProbability) + minRepeatLength;
}

double getRepeatLengthProbForHaplotype(int64_t repeatLength, double *logProbabilitiesHap1, double *logProbabilitiesHap2,
                                       int64_t minRepeatLength, double logProbMLHap2, double logProbSubstitution) {
    double logProbHap2Same = logProbabilitiesHap2[repeatLength - minRepeatLength];
    return logProbabilitiesHap1[repeatLength - minRepeatLength] +
           ((logProbHap2Same > logProbMLHap2 + logProbSubstitution) ? logProbHap2Same : logProbMLHap2 +
                                                                                        logProbSubstitution);
}

int64_t getMax(double *values, int64_t length,
               double *maxValue) {
    // Calc the range of repeat observations
    assert(length > 0);
    double max = values[0];
    int64_t maxIndex = 0;
    for (int64_t i = 1; i < length; i++) {
        if (values[i] > max) {
            max = values[i];
            maxIndex = i;
        }
    }
    *maxValue = max;
    return maxIndex;
}

int64_t repeatSubMatrix_getPhasedMLRepeatCount(RepeatSubMatrix *repeatSubMatrix, int64_t existingRepeatCount,
        Symbol base, stList *observations, stList *bamChunkReads, double *logProbability, stSet *readsBelongingToHap1,
        stSet *readsBelongingToHap2, PolishParams *params) {
    // Calculate range of repeat counts observed
    int64_t minRepeatLength, maxRepeatLength;
    repeatSubMatrix_getMinAndMaxRepeatCountObservations(repeatSubMatrix, observations,
                                                        bamChunkReads, &minRepeatLength, &maxRepeatLength);

    if (minRepeatLength == repeatSubMatrix->maximumRepeatLength) {
        *logProbability = LOG_ZERO;
        return 0; // Case we have no valid observations, so assume repeat length of 0
    }

    // Split observations between haplotypes
    stList *observationsHap1 = stList_construct();
    stList *observationsHap2 = stList_construct();
    for (int64_t i = 0; i < stList_length(observations); i++) {
        PoaBaseObservation *observation = stList_get(observations, i);
        BamChunkRead *read = stList_get(bamChunkReads, observation->readNo);
        stList_append(stSet_search(readsBelongingToHap1, read) != NULL ? observationsHap1 : observationsHap2,
                      observation);
    }

    // Get probs for hap 1
    double logProbabilitiesHap1[repeatSubMatrix->maximumRepeatLength];
    repeatSubMatrix_getRepeatCountProbs(repeatSubMatrix, base, observationsHap1,
                                        bamChunkReads, logProbabilitiesHap1, minRepeatLength, maxRepeatLength);

    // Get probs for hap 2
    double logProbabilitiesHap2[repeatSubMatrix->maximumRepeatLength];
    repeatSubMatrix_getRepeatCountProbs(repeatSubMatrix, base, observationsHap2,
                                        bamChunkReads, logProbabilitiesHap2, minRepeatLength, maxRepeatLength);

    // Get ML prob for haplotype 2
    double logProbMLHap2;
    int64_t mLRepeatLengthHap2 =
            getMax(logProbabilitiesHap2, maxRepeatLength - minRepeatLength + 1, &logProbMLHap2) + minRepeatLength;

    // Calculate ML probability of repeat length for haplotype1
    double mlLogProb = getRepeatLengthProbForHaplotype(minRepeatLength, logProbabilitiesHap1, logProbabilitiesHap2,
                                                       minRepeatLength, logProbMLHap2,
                                                       log(params->hetRunLengthSubstitutionProbability));
    int64_t mlRepeatLength = minRepeatLength;
    for (int64_t i = minRepeatLength + 1; i <= maxRepeatLength; i++) {
        double p = getRepeatLengthProbForHaplotype(i, logProbabilitiesHap1, logProbabilitiesHap2,
                                                   minRepeatLength, logProbMLHap2,
                                                   log(params->hetRunLengthSubstitutionProbability));
        if (p >= mlLogProb) {
            mlLogProb = p;
            mlRepeatLength = i;
        }
    }
    *logProbability = mlLogProb;

    if (mlRepeatLength != existingRepeatCount) {
        st_logDebug(
                "Got %i repeat length, other hap repeat length %i, log prob:%f, log prob hap1: %f log prob hap2: %f (min rl: %i max: %i) (total obs: %i, hap1: %i, hap2 : %i)\n",
                (int) mlRepeatLength, (int) mLRepeatLengthHap2, (float) *logProbability,
                (float) logProbabilitiesHap1[mlRepeatLength - minRepeatLength],
                (float) logProbMLHap2, (int) minRepeatLength, (int) maxRepeatLength,
                (int) stList_length(observations), (int) stList_length(observationsHap1),
                (int) stList_length(observationsHap2));
    }

    // Cleanup
    stList_destruct(observationsHap1);
    stList_destruct(observationsHap2);

    return mlRepeatLength;
}