/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include <sys/stat.h>
#include <sonLibListPrivate.h>
#include <helenFeatures.h>

char *getTimeDescriptorFromSeconds(int64_t seconds) {
    int64_t minutes = (int64_t) (seconds / 60);
    int64_t hours = (int64_t) (minutes / 60);
    char *timeDescriptor;

    if (hours > 0) {
        timeDescriptor = stString_print("%"PRId64"h %"PRId64"m", hours,
                minutes - (hours * 60));
    } else if (minutes > 0) {
        timeDescriptor = stString_print("%"PRId64"m %"PRId64"s", minutes,
                seconds - (minutes * 60));
    } else {
        timeDescriptor = stString_print("%"PRId64"s", seconds);
    }
    return timeDescriptor;
}

stHash *parseReferenceSequences(char *referenceFastaFile) {
    /*
     * Get hash of reference sequence names in fasta to their sequences, doing some munging on the sequence names.
     */
    st_logCritical("> Parsing reference sequences from file: %s\n", referenceFastaFile);
    FILE *fh = fopen(referenceFastaFile, "r");
    stHash *referenceSequences = fastaReadToMap(fh);  //valgrind says blocks from this allocation are "still reachable"
    fclose(fh);
    // log names and transform (if necessary)
    stList *refSeqNames = stHash_getKeys(referenceSequences);
    int64_t origRefSeqLen = stList_length(refSeqNames);
    st_logDebug("\tReference contigs: \n");
    for (int64_t i = 0; i < origRefSeqLen; ++i) {
        char *fullRefSeqName = (char *) stList_get(refSeqNames, i);
        st_logDebug("\t\t%s\n", fullRefSeqName);
        char refSeqName[128] = "";
        if (sscanf(fullRefSeqName, "%s", refSeqName) == 1 && !stString_eq(fullRefSeqName, refSeqName)) {
            // this transformation is necessary for cases where the reference has metadata after the contig name:
            // >contig001 length=1000 date=1999-12-31
            char *newKey = stString_copy(refSeqName);
            char *refSeq = stHash_search(referenceSequences, fullRefSeqName);
            stHash_insert(referenceSequences, newKey, refSeq);
            stHash_removeAndFreeKey(referenceSequences, fullRefSeqName);
            st_logDebug("\t\t\t-> %s\n", newKey);
        }
    }
    stList_destruct(refSeqNames);

    return referenceSequences;
}

char *getFileBase(char *base, char *defawlt) {
    struct stat fileStat;
    int64_t rc = stat(base, &fileStat);
    if (S_ISDIR(fileStat.st_mode)) {
        if (optarg[strlen(base) - 1] == '/') optarg[strlen(base) - 1] = '\0';
        return stString_print("%s/%s", base, defawlt);
    } else {
        return stString_copy(base);
    }
}

RleString *bamChunk_getReferenceSubstring(BamChunk *bamChunk, stHash *referenceSequences, Params *params) {
    /*
     * Get corresponding substring of the reference for a given bamChunk.
     */
    char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
    if (fullReferenceString == NULL) {
        st_logCritical("> ERROR: Reference sequence missing from reference map: %s \n", bamChunk->refSeqName);
        return NULL;
    }
    int64_t refLen = strlen(fullReferenceString);
    char *referenceString = stString_getSubString(fullReferenceString, bamChunk->chunkBoundaryStart,
                                                  (refLen < bamChunk->chunkBoundaryEnd ? refLen : bamChunk->chunkBoundaryEnd) - bamChunk->chunkBoundaryStart);

    int64_t subRefLen = strlen(referenceString);
    // TODO: Decide where the proper place to do this is
    for (int64_t i = 0; i < subRefLen; i++) {
        referenceString[i] = (char) toupper(referenceString[i]);
    }

    RleString *rleRef = params->polishParams->useRunLengthEncoding ?
                        rleString_construct(referenceString) : rleString_construct_no_rle(referenceString);

    free(referenceString);

    return rleRef;
}


uint64_t *getPaddedHaplotypeString(uint64_t *hap, stGenomeFragment *gf, BubbleGraph *bg, Params *params) {
    /*
     * Pads a haplotype string from the genome fragment to account for any missing prefix or suffix.
     */
    uint64_t *paddedHap = bubbleGraph_getConsensusPath(bg, params->polishParams);

    for(uint64_t i=0; i<gf->length; i++) {
        paddedHap[i+gf->refStart] = hap[i];
    }

    return paddedHap;
}

stSet *bamChunkRead_to_readName(stSet *bamChunkReads) {
    stSet *readNames = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stSetIterator *itor = stSet_getIterator(bamChunkReads);
    BamChunkRead *bcr = NULL;
    while ((bcr = stSet_getNext(itor)) != NULL) {
        stSet_insert(readNames, stString_copy(bcr->readName));
    }
    stSet_destructIterator(itor);
    return readNames;
}

stList *copyListOfIntTuples(stList *toCopy) {
    stList *copy = stList_construct3(0, toCopy->destructElement);
    for (int64_t j = 0; j < stList_length(toCopy); j++) {
        stIntTuple *tupleFrom = stList_get(toCopy, j);
        int64_t n = stIntTuple_length(tupleFrom);
        stIntTuple *tupleTo = stIntTuple_constructN(n, &(tupleFrom[1]));
        stList_append(copy, tupleTo);

        //todo delete
        assert(stIntTuple_length(tupleTo) == n);
        for (int64_t k = 0; k < n; k++) {
            assert(stIntTuple_get(tupleFrom, k) == stIntTuple_get(tupleTo, k));
        }
    }
    return copy;
}

void saveAlignmentWeightsAndScores(Poa *poa, stList *reads, double *weights, double *score, int64_t length) {
    for (int64_t p = 0; p < length; p++) {
        PoaNode *node = stList_get(poa->nodes, p);

        // match/mismatch
        for (int64_t o = 0; o < stList_length(node->observations); o++) {
            PoaBaseObservation *observation = stList_get(node->observations, o);
            BamChunkRead *read = stList_get(reads, observation->readNo);

            weights[observation->readNo] += observation->weight;
            if (node->base == read->rleRead->rleString[observation->offset]) {
                /*int nrc = (int) node->repeatCount;
                int rrc = (int) read->rleRead->repeatCounts[observation->offset];
                score[observation->readNo] += (2 - abs(nrc - rrc)/((nrc+rrc)/2.0) ) * observation->weight;*/
                score[observation->readNo] += 2 * observation->weight;
            } else {
                score[observation->readNo] -= observation->weight;
            }

        }

        // insert
        for (int64_t i = 0; i < stList_length(node->inserts); i++) {
            PoaInsert *insert = stList_get(node->inserts, i);
            for (int64_t o = 0; o < stList_length(insert->observations); o++) {
                PoaBaseObservation *observation = stList_get(insert->observations, o);

                weights[observation->readNo] += observation->weight;
                score[observation->readNo] -= observation->weight;
            }
        }

        // delete
        for (int64_t d = 0; d < stList_length(node->deletes); d++) {
            PoaDelete *delete = stList_get(node->deletes, d);
            for (int64_t o = 0; o < stList_length(delete->observations); o++) {
                PoaBaseObservation *observation = stList_get(delete->observations, o);

                weights[observation->readNo] += observation->weight;
                score[observation->readNo] -= observation->weight;
            }
        }

    }

}




void rankReadPoaAlignments(stList *reads, Poa *poa_hap1, Poa *poa_hap2, stSet *hap1Reads, stSet *hap2Reads) {
    int64_t length = stList_length(reads);
    double *totalReadWeights_hap1 = st_calloc(length, sizeof(double));
    double *totalReadWeights_hap2 = st_calloc(length, sizeof(double));
    double *totalReadScore_hap1 = st_calloc(length, sizeof(double));
    double *totalReadScore_hap2 = st_calloc(length, sizeof(double));
    char *logIdentifier = getLogIdentifier();

    /*if (stList_length(poa_hap1->nodes) != stList_length(poa_hap2->nodes)) {
        st_logInfo(" %s Ranked POAs have different sizes: %"PRId64", %"PRId64"\n", logIdentifier,
                stList_length(poa_hap1->nodes), stList_length(poa_hap2->nodes));
    }
    int64_t diffNodes = 0;
    for (int64_t i = 0; i < stList_length(poa_hap1->nodes); i++) {
        PoaNode *node1 = stList_get(poa_hap1->nodes, i);
        PoaNode *node2 = stList_get(poa_hap2->nodes, i);
        char pos1[6] = "     ";
        char pos2[6] = "     ";
        pos1[0] = node1->base;
        pos2[0] = node2->base;
        pos1[1] = (char) ('0' + node1->repeatCount);
        pos2[1] = (char) ('0' + node2->repeatCount);
        pos1[3] = (char) ('0' + stList_length(node1->inserts));
        pos2[3] = (char) ('0' + stList_length(node2->inserts));
        pos1[4] = (char) ('0' + stList_length(node1->deletes));
        pos2[4] = (char) ('0' + stList_length(node2->deletes));

        if (!stString_eq(pos1, pos2)) {
            st_logInfo(" %s Ranked POA diff: %6"PRId64" %s %s\n", logIdentifier, i, pos1, pos2);
            diffNodes++;
            if (diffNodes >= 8) {
                st_logInfo(" %s Ranked POA diff: ...\n", logIdentifier);
                break;
            }
        }
    }*/

    saveAlignmentWeightsAndScores(poa_hap1, reads, totalReadWeights_hap1, totalReadScore_hap1, length);
    saveAlignmentWeightsAndScores(poa_hap2, reads, totalReadWeights_hap2, totalReadScore_hap2, length);

    int64_t scoreAndWeightAreAlignedCount = 0;
    int64_t equalWeightCount = 0;
    int64_t equalScoreCount = 0;
    int64_t unclassifiedCount = 0;
    for (int64_t i = 0; i < length; i++) {
        int hapByWeight = 0;
        int hapByScore = 0;
        if (totalReadScore_hap1[i] == totalReadScore_hap2[i]) {
            equalScoreCount++;
        } else if (totalReadScore_hap1[i] > totalReadScore_hap2[i]) {
            hapByWeight = 1;
        } else {
            hapByWeight = 2;
        }
        if (totalReadWeights_hap1[i] == totalReadWeights_hap2[i]) {
            equalWeightCount++;
        } else if (totalReadWeights_hap1[i] >= totalReadWeights_hap2[i]) {
            hapByScore = 1;
        } else {
            hapByScore = 2;
        }

        // quantify score v total weight
        if (hapByScore == hapByWeight) {
            scoreAndWeightAreAlignedCount++;
        }

        // classify by score (probably better)
        if (hapByScore == 0) {
            unclassifiedCount++;
        } else {
            stSet_insert(hapByScore == 1 ? hap1Reads : hap2Reads, ((BamChunkRead *) stList_get(reads, i))->readName);
        }
    }

    // loggit
    st_logInfo(" %s Of %"PRId64" reads, %.4f had score and weight aligned, %.4f had score equal, %.4f had equal weight, and %.4f were unclassified\n",
            logIdentifier, length, 1.0*scoreAndWeightAreAlignedCount/length, 1.0*equalScoreCount/length,
            1.0*equalWeightCount/length, 1.0*unclassifiedCount/length);

    // cleanup
    free(logIdentifier);
    free(totalReadWeights_hap1);
    free(totalReadWeights_hap2);
    free(totalReadScore_hap1);
    free(totalReadScore_hap2);
}


void rankReadPoaAlignments2(stList *reads, Poa *poa_hap1, Poa *poa_hap2, stSet *hap1Reads, stSet *hap2Reads,
                            PolishParams *polishParams) {
    int64_t length = stList_length(reads);
    double *totalReadScore_hap1 = st_calloc(length, sizeof(double));
    double *totalReadScore_hap2 = st_calloc(length, sizeof(double));
    char *logIdentifier = getLogIdentifier();
    RleString *polishedRleConsensusH1 = rleString_copy(poa_hap1->refString);
    RleString *polishedRleConsensusH2 = rleString_copy(poa_hap2->refString);
    // align POAs and rank based on misaligned nodes
    double score = 0;
    stList *alignedPairs = alignConsensusAndTruthRLEWithKmerAnchors(polishedRleConsensusH1, polishedRleConsensusH2,
                                                                    &score, polishParams);
    /*uint16_t score;
    char *polishedRawConsensusH1 = rleString_expand(polishedRleConsensusH1);
    char *polishedRawConsensusH2 = rleString_expand(polishedRleConsensusH2);
    uint64_t *polishedConsensusH1_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(polishedRleConsensusH1);
    uint64_t *polishedConsensusH2_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(polishedRleConsensusH2);
    uint16_t sswScore = 0;
    stList *alignedPairsRawSSW = alignConsensusAndTruthSSW(polishedRawConsensusH1, polishedRawConsensusH2, &sswScore);
    stList *alignedPairs = runLengthEncodeAlignment(alignedPairsRawSSW, polishedConsensusH1_nonRleToRleCoordinateMap,
            polishedConsensusH2_nonRleToRleCoordinateMap);
    free(polishedConsensusH1_nonRleToRleCoordinateMap);
    free(polishedConsensusH2_nonRleToRleCoordinateMap);
    free(polishedRawConsensusH1);
    free(polishedRawConsensusH2);*/

    printMEAAlignment2(polishedRleConsensusH1, polishedRleConsensusH2, alignedPairs);

    stListIterator *alignmentItor = stList_getIterator(alignedPairs);
    stIntTuple *currAlign = stList_getNext(alignmentItor);
    int64_t posH1 = stIntTuple_get(currAlign, 0);
    int64_t posH2 = stIntTuple_get(currAlign, 1);
    int64_t insertH1 = 0;
    int64_t insertH2 = 0;
    int64_t mismatch = 0;
    while (TRUE) {
        if (currAlign == NULL) break;
        int64_t currAlignPosHap1 = stIntTuple_get(currAlign, 0);
        int64_t currAlignPosHap2 = stIntTuple_get(currAlign, 1);
        // H2 gap / H1 insert
        if (posH1 < currAlignPosHap1) {
            insertH1++;
            PoaNode *insertNodeX = stList_get(poa_hap1->nodes, posH1);
            for (int64_t o = 0; o < stList_length(insertNodeX->observations); o++) {
                PoaBaseObservation *observation = stList_get(insertNodeX->observations, o);
                BamChunkRead *read = stList_get(reads, observation->readNo);
                if (insertNodeX->base == read->rleRead->rleString[observation->offset]) {
                    totalReadScore_hap1[observation->readNo] += observation->weight;
                } else {
                    totalReadScore_hap1[observation->readNo] -= observation->weight;
                }
            }
            posH1++;
        }
            // H1 gap / H2 insert
        else if (posH2 < currAlignPosHap2) {
            insertH2++;
            PoaNode *insertNodeY = stList_get(poa_hap1->nodes, posH2);
            for (int64_t o = 0; o < stList_length(insertNodeY->observations); o++) {
                PoaBaseObservation *observation = stList_get(insertNodeY->observations, o);
                BamChunkRead *read = stList_get(reads, observation->readNo);
                if (insertNodeY->base == read->rleRead->rleString[observation->offset]) {
                    totalReadScore_hap2[observation->readNo] += observation->weight;
                } else {
                    totalReadScore_hap2[observation->readNo] -= observation->weight;
                }
            }
            posH2++;
        }
            // match
        else if (posH1 == currAlignPosHap1 && posH2 == currAlignPosHap2) {
            PoaNode *matchNodeH1 = stList_get(poa_hap1->nodes, posH1);
            PoaNode *matchNodeH2 = stList_get(poa_hap2->nodes, posH2);
            if (matchNodeH1->base != matchNodeH2->base) {
                mismatch++;
                for (int64_t o = 0; o < stList_length(matchNodeH1->observations); o++) {
                    PoaBaseObservation *observation = stList_get(matchNodeH1->observations, o);
                    BamChunkRead *read = stList_get(reads, observation->readNo);
                    if (matchNodeH1->base == read->rleRead->rleString[observation->offset]) {
                        totalReadScore_hap1[observation->readNo] += observation->weight;
                    } else {
                        totalReadScore_hap1[observation->readNo] -= observation->weight;
                    }
                }
                for (int64_t o = 0; o < stList_length(matchNodeH2->observations); o++) {
                    PoaBaseObservation *observation = stList_get(matchNodeH2->observations, o);
                    BamChunkRead *read = stList_get(reads, observation->readNo);
                    if (matchNodeH2->base == read->rleRead->rleString[observation->offset]) {
                        totalReadScore_hap2[observation->readNo] += observation->weight;
                    } else {
                        totalReadScore_hap2[observation->readNo] -= observation->weight;
                    }
                }
            }
            posH1++;
            posH2++;
            currAlign = stList_getNext(alignmentItor);
        }
            // should never happen
        else {
            assert(FALSE);
        }
    }
    int64_t totalNoScoreLength = 0;
    int64_t noScoreCount = 0;
    int64_t unclassifiedCount = 0;
    int64_t hap1Count = 0;
    int64_t hap2Count = 0;
    for (int64_t i = 0; i < length; i++) {
        if (totalReadScore_hap1[i] < totalReadScore_hap2[i]) {
            stSet_insert(hap2Reads, ((BamChunkRead *) stList_get(reads, i))->readName);
            hap1Count++;
        } else if (totalReadScore_hap1[i] > totalReadScore_hap2[i]) {
            stSet_insert(hap1Reads, ((BamChunkRead *) stList_get(reads, i))->readName);
            hap2Count++;
        } else {
            if (totalReadScore_hap1[i] == 0) {
                totalNoScoreLength += ((BamChunkRead *) stList_get(reads, i))->rleRead->nonRleLength;
                noScoreCount++;
            }
            unclassifiedCount++;
        }
    }
    // loggit
    st_logInfo(" %s Filtered read haplotyping performed over %d aligned positions, with %"PRId64" H1 inserts, %"PRId64" H2 inserts, and %"PRId64" mismatches\n",
               logIdentifier, stList_length(alignedPairs), insertH1, insertH2, mismatch);
    st_logInfo(" %s Of %"PRId64" filtered reads: %"PRId64" (%.2f) were hap1, %"PRId64" (%.2f) were hap2, %"PRId64" (%.2f) were unclassified with %"PRId64" (%.2f) having no score (avg len %"PRId64").\n",
               logIdentifier, length, hap1Count, 1.0*hap1Count/length, hap2Count, 1.0*hap2Count/length,
               unclassifiedCount, 1.0*unclassifiedCount/length, noScoreCount,
               1.0*noScoreCount/(unclassifiedCount == 0 ? 1 : unclassifiedCount), totalNoScoreLength / noScoreCount);
    stList_destructIterator(alignmentItor);
    stList_destruct(alignedPairs);
    rleString_destruct(polishedRleConsensusH1);
    rleString_destruct(polishedRleConsensusH2);
    free(totalReadScore_hap1);
    free(totalReadScore_hap2);
}