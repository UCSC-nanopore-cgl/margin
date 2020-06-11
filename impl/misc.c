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

void assignFilteredReadsToHaplotypes2(BubbleGraph *bg, uint64_t *hap1, uint64_t *hap2, RleString *rleReference,
                                      stList *filteredReads, stList *filteredAlignments,
                                      stSet *hap1Reads, stSet *hap2Reads, Params *params) {
    // quick fail
    int64_t length = stList_length(filteredReads);
    if (length == 0) return;

    // get POA structure for filtering
    Poa *filteredPoa = poa_realign(filteredReads, filteredAlignments, rleReference, params->polishParams);
    Poa *filteredPoa_hap1 = bubbleGraph_getNewPoa(bg, hap1, filteredPoa, filteredReads, params);
    Poa *filteredPoa_hap2 = bubbleGraph_getNewPoa(bg, hap2, filteredPoa, filteredReads, params);

    // align POAs to find het sites
    RleString *polishedRleConsensusH1 = rleString_copy(filteredPoa_hap1->refString);
    RleString *polishedRleConsensusH2 = rleString_copy(filteredPoa_hap2->refString);
    double score = 0;
    stList *alignedPairs = alignConsensusAndTruthRLEWithKmerAnchors(polishedRleConsensusH1, polishedRleConsensusH2,
                                                                    &score, params->polishParams);

    // quick failure check
    char *logIdentifier = getLogIdentifier();
    if (alignedPairs == NULL || stList_length(alignedPairs) == 0) {
        st_logInfo(" %s Could not get polished consensus alignment for haplotyping filtered reads.\n", logIdentifier);

        // cleanup
        stList_destruct(alignedPairs);
        rleString_destruct(polishedRleConsensusH1);
        rleString_destruct(polishedRleConsensusH2);
        free(logIdentifier);
        return;
    }

    // score reads based on alignment at het sites
    double *totalReadScore_hap1 = st_calloc(length, sizeof(double));
    double *totalReadScore_hap2 = st_calloc(length, sizeof(double));
    stListIterator *alignmentItor = stList_getIterator(alignedPairs);
    stIntTuple *currAlign = stList_getNext(alignmentItor);
    int64_t posH1 = currAlign == NULL ? 0 : stIntTuple_get(currAlign, 0);
    int64_t posH2 = currAlign == NULL ? 0 : stIntTuple_get(currAlign, 1);
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
            PoaNode *insertNodeX = stList_get(filteredPoa_hap1->nodes, posH1);
            for (int64_t o = 0; o < stList_length(insertNodeX->observations); o++) {
                PoaBaseObservation *observation = stList_get(insertNodeX->observations, o);
                BamChunkRead *read = stList_get(filteredReads, observation->readNo);
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
            PoaNode *insertNodeY = stList_get(filteredPoa_hap2->nodes, posH2);
            for (int64_t o = 0; o < stList_length(insertNodeY->observations); o++) {
                PoaBaseObservation *observation = stList_get(insertNodeY->observations, o);
                BamChunkRead *read = stList_get(filteredReads, observation->readNo);
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
            PoaNode *matchNodeH1 = stList_get(filteredPoa_hap1->nodes, posH1);
            PoaNode *matchNodeH2 = stList_get(filteredPoa_hap2->nodes, posH2);
            if (matchNodeH1->base != matchNodeH2->base) {
                mismatch++;
                for (int64_t o = 0; o < stList_length(matchNodeH1->observations); o++) {
                    PoaBaseObservation *observation = stList_get(matchNodeH1->observations, o);
                    BamChunkRead *read = stList_get(filteredReads, observation->readNo);
                    if (matchNodeH1->base == read->rleRead->rleString[observation->offset]) {
                        totalReadScore_hap1[observation->readNo] += observation->weight;
                    } else {
                        totalReadScore_hap1[observation->readNo] -= observation->weight;
                    }
                }
                for (int64_t o = 0; o < stList_length(matchNodeH2->observations); o++) {
                    PoaBaseObservation *observation = stList_get(matchNodeH2->observations, o);
                    BamChunkRead *read = stList_get(filteredReads, observation->readNo);
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

    // get scores and save to appropriate sets
    int64_t totalNoScoreLength = 0;
    int64_t noScoreCount = 0;
    int64_t unclassifiedCount = 0;
    int64_t hap1Count = 0;
    int64_t hap2Count = 0;
    for (int64_t i = 0; i < length; i++) {
        if (totalReadScore_hap1[i] < totalReadScore_hap2[i]) {
            stSet_insert(hap2Reads, stList_get(filteredReads, i));
            hap1Count++;
        } else if (totalReadScore_hap1[i] > totalReadScore_hap2[i]) {
            stSet_insert(hap1Reads, stList_get(filteredReads, i));
            hap2Count++;
        } else {
            if (totalReadScore_hap1[i] == 0) {
                totalNoScoreLength += ((BamChunkRead *) stList_get(filteredReads, i))->rleRead->nonRleLength;
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
               1.0*noScoreCount/(unclassifiedCount == 0 ? 1 : unclassifiedCount),
               totalNoScoreLength / (noScoreCount == 0 ? 1 : noScoreCount));

    // cleanup
    stList_destructIterator(alignmentItor);
    stList_destruct(alignedPairs);
    rleString_destruct(polishedRleConsensusH1);
    rleString_destruct(polishedRleConsensusH2);
    free(totalReadScore_hap1);
    free(totalReadScore_hap2);
    poa_destruct(filteredPoa);
    poa_destruct(filteredPoa_hap1);
    poa_destruct(filteredPoa_hap2);
    free(logIdentifier);
}