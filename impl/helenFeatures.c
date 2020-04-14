//
// Created by tpesout on 3/29/19.
//

#ifdef _HDF5

#include "margin.h"
#include "htsIntegration.h"
#include "helenFeatures.h"
#include "ssw.h"
#include <hdf5.h>

PoaFeatureSimpleWeight *PoaFeature_SimpleWeight_construct(int64_t refPos, int64_t insPos) {
    PoaFeatureSimpleWeight *feature = st_calloc(1, sizeof(PoaFeatureSimpleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->label = '\0';
    feature->nextInsert = NULL;
    return feature;
}

void PoaFeature_SimpleWeight_destruct(PoaFeatureSimpleWeight *feature) {
    if (feature->nextInsert != NULL) {
        PoaFeature_SimpleWeight_destruct(feature->nextInsert);
    }
    free(feature);
}

PoaFeatureSplitRleWeight *PoaFeature_SplitRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
                                                              int64_t maxRunLength) {
    PoaFeatureSplitRleWeight *feature = st_calloc(1, sizeof(PoaFeatureSplitRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->weights = st_calloc(((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2, sizeof(double));
    return feature;
}

void PoaFeature_SplitRleWeight_destruct(PoaFeatureSplitRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_SplitRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_SplitRleWeight_destruct(feature->nextInsert);
    }
    free(feature->weights);
    free(feature);
}


PoaFeatureChannelRleWeight *PoaFeature_ChannelRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
                                                                  int64_t maxRunLength) {
    PoaFeatureChannelRleWeight *feature = st_calloc(1, sizeof(PoaFeatureChannelRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->nucleotideWeights = st_calloc((SYMBOL_NUMBER) * 2, sizeof(double));
    feature->runLengthWeights = st_calloc((SYMBOL_NUMBER - 1) * (1 + maxRunLength) * 2, sizeof(double));
    return feature;
}

void PoaFeature_ChannelRleWeight_destruct(PoaFeatureChannelRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_ChannelRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_ChannelRleWeight_destruct(feature->nextInsert);
    }
    free(feature->nucleotideWeights);
    free(feature->runLengthWeights);
    free(feature);
}

PoaFeatureDiploidRleWeight *PoaFeature_DiploidRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos,
        int64_t maxRunLength) {
    PoaFeatureDiploidRleWeight *feature = st_calloc(1, sizeof(PoaFeatureDiploidRleWeight));
    feature->refPosition = refPos;
    feature->insertPosition = insPos;
    feature->runLengthPosition = rlPos;
    feature->labelChar = '\0';
    feature->labelRunLength = 0;
    feature->nextRunLength = NULL;
    feature->nextInsert = NULL;
    feature->maxRunLength = maxRunLength;
    feature->weightsHOn = st_calloc(((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2, sizeof(double));
    feature->weightsHOff = st_calloc(((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2, sizeof(double));
    return feature;
}
void PoaFeature_DiploidRleWeight_destruct(PoaFeatureDiploidRleWeight *feature) {
    if (feature->nextRunLength != NULL) {
        PoaFeature_DiploidRleWeight_destruct(feature->nextRunLength);
    }
    if (feature->nextInsert != NULL) {
        PoaFeature_DiploidRleWeight_destruct(feature->nextInsert);
    }
    free(feature->weightsHOn);
    free(feature->weightsHOff);
    free(feature);
}


int PoaFeature_SimpleWeight_charIndex(Symbol character, bool forward) {
    int pos = character * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    return pos;
}

int PoaFeature_SimpleWeight_gapIndex(bool forward) {
    int pos = POAFEATURE_SYMBOL_GAP_POS * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    assert(pos < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    return pos;
}

int PoaFeature_SplitRleWeight_charIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    assert(runLength >= 0);
    assert(runLength <= maxRunLength);
    int pos = (character * ((int) maxRunLength + 1) + runLength) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_SplitRleWeight_gapIndex(int64_t maxRunLength, bool forward) {
    int pos = ((SYMBOL_NUMBER - 1) * ((int) maxRunLength + 1)) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_ChannelRleWeight_charNuclIndex(Symbol character, bool forward) {
    int pos = character * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_ChannelRleWeight_gapNuclIndex(bool forward) {
    int pos = (SYMBOL_NUMBER_NO_N) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_ChannelRleWeight_charRLIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    assert(runLength >= 0);
    assert(runLength <= maxRunLength);
    int pos = (character * ((int) maxRunLength + 1) + runLength) * 2 + (forward ? POS_STRAND_IDX : NEG_STRAND_IDX);
    return pos;
}

int PoaFeature_DiploidRleWeight_charIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward) {
    return PoaFeature_SplitRleWeight_charIndex(maxRunLength, character, runLength, forward);
}
int PoaFeature_DiploidRleWeight_gapIndex(int64_t maxRunLength, bool forward) {
    return PoaFeature_SplitRleWeight_gapIndex(maxRunLength, forward);
}


void handleHelenFeatures(
        // global params
        HelenFeatureType helenFeatureType, BamChunker *trueReferenceBamChunker,
        int64_t splitWeightMaxRunLength, void **helenHDF5Files, bool fullFeatureOutput,
        char *trueReferenceBam, Params *params,

        // chunk params
        char *logIdentifier, int64_t chunkIdx, BamChunk *bamChunk, Poa *poa, stList *bamChunkReads,
        char *polishedConsensusString, RleString *polishedRleConsensus) {

    st_logInfo(">%s Performing feature generation for chunk.\n", logIdentifier);

    // get filename
    char *helenFeatureOutfileBase = NULL;
    switch (helenFeatureType) {
        case HFEAT_SIMPLE_WEIGHT:
            helenFeatureOutfileBase = stString_print("simpleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        case HFEAT_SPLIT_RLE_WEIGHT:
            // name of folder, not of file
            helenFeatureOutfileBase = stString_print("splitRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        case HFEAT_CHANNEL_RLE_WEIGHT:
            // name of folder, not of file
            helenFeatureOutfileBase = stString_print("channelRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    // necessary to annotate poa with truth (if true reference BAM has been specified)
    stList *trueRefAlignment = NULL;
    RleString *trueRefRleString = NULL;
    bool validReferenceAlignment = FALSE;

    // get reference chunk
    if (trueReferenceBam != NULL) {
        // get alignment of true ref to assembly
        stList *trueRefReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *unused = stList_construct3(0, (void (*)(void *)) stList_destruct);
        // construct new chunk
        BamChunk *trueRefBamChunk = bamChunk_copyConstruct(bamChunk);
        trueRefBamChunk->parent = trueReferenceBamChunker;
        // get true ref as "read"
        uint32_t trueAlignmentCount = convertToReadsAndAlignments(trueRefBamChunk, NULL, trueRefReads, unused);

        // poor man's "do we have a unique alignment"
        if (trueAlignmentCount == 1) {
            BamChunkRead *trueRefRead = stList_get(trueRefReads, 0);
            char *trueRefExpanded = rleString_expand(trueRefRead->rleRead);

            uint16_t score;
            stList *trueRefAlignmentRawSpace = alignConsensusAndTruthSSW(polishedConsensusString, trueRefExpanded,
                                                                         &score);
            if (st_getLogLevel() == debug) {
                printMEAAlignment(polishedConsensusString, trueRefExpanded,
                                  strlen(polishedConsensusString), strlen(trueRefExpanded),
                                  trueRefAlignmentRawSpace, NULL, NULL);
            }


            // convert to rleSpace if appropriate
            if (params->polishParams->useRunLengthEncoding) {
                trueRefRleString = rleString_construct(trueRefExpanded);

                uint64_t *polishedRleConsensus_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(polishedRleConsensus);
                uint64_t *trueRefRleString_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(trueRefRleString);
                trueRefAlignment = runLengthEncodeAlignment(trueRefAlignmentRawSpace, polishedRleConsensus_nonRleToRleCoordinateMap,
                                                            trueRefRleString_nonRleToRleCoordinateMap);
                free(polishedRleConsensus_nonRleToRleCoordinateMap);
                free(trueRefRleString_nonRleToRleCoordinateMap);

                if (st_getLogLevel() == debug) {
                    printMEAAlignment(polishedRleConsensus->rleString, trueRefRleString->rleString,
                                      strlen(polishedRleConsensus->rleString),
                                      strlen(trueRefRleString->rleString),
                                      trueRefAlignment, polishedRleConsensus->repeatCounts,
                                      trueRefRleString->repeatCounts);
                }
                stList_destruct(trueRefAlignmentRawSpace);
            } else {
                trueRefRleString = rleString_construct_no_rle(trueRefExpanded);
                trueRefAlignment = trueRefAlignmentRawSpace;
            }


            // we found a single alignment of reference
            double refLengthRatio = 1.0 * strlen(trueRefExpanded) / strlen(polishedConsensusString);
            double alnLengthRatio = 1.0 * stList_length(trueRefAlignment) / polishedRleConsensus->length;
            int refLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - refLengthRatio)));
            int alnLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - alnLengthRatio)));
            if (stList_length(trueRefAlignment) > 0 && refLengthRatioHundredthsOffOne < 10 &&
                alnLengthRatioHundredthsOffOne < 10) {
                validReferenceAlignment = TRUE;
            } else {
                st_logInfo(" %s True reference alignment QC failed:  polished length %"PRId64", true ref length"
                           " ratio (true/polished) %f, aligned pairs length ratio (true/polished): %f\n",
                           logIdentifier, polishedRleConsensus->length, refLengthRatio, alnLengthRatio);
            }

            //cleanup
            free(trueRefExpanded);
        }

        stList_destruct(trueRefReads);
        stList_destruct(unused);
        bamChunk_destruct(trueRefBamChunk);
    }

    // either write it, or note that we failed to find a valid reference alignment
    if (trueReferenceBam != NULL && !validReferenceAlignment) {
        st_logInfo(" %s No valid reference alignment was found, skipping HELEN feature output.\n", logIdentifier);
    } else {
        st_logInfo(" %s Writing HELEN features with filename base: %s\n", logIdentifier, helenFeatureOutfileBase);

        // write the actual features (type dependent)
        poa_writeHelenFeatures(helenFeatureType, poa, bamChunkReads, helenFeatureOutfileBase,
                               bamChunk, trueRefAlignment, polishedRleConsensus, trueRefRleString, fullFeatureOutput,
                               splitWeightMaxRunLength, (HelenFeatureHDF5FileInfo **) helenHDF5Files);

        // write the polished chunk in fasta format
        if (fullFeatureOutput) {
            char *chunkPolishedRefFilename = stString_print("%s.fa", helenFeatureOutfileBase);
            char *chunkPolishedRefContigName = stString_print("%s\t%"PRId64"\t%"PRId64"\t%s",
                                                              bamChunk->refSeqName,
                                                              bamChunk->chunkBoundaryStart,
                                                              bamChunk->chunkBoundaryEnd,
                                                              helenFeatureOutfileBase);
            FILE *chunkPolishedRefOutFh = fopen(chunkPolishedRefFilename, "w");
            fastaWrite(polishedConsensusString, chunkPolishedRefContigName, chunkPolishedRefOutFh);
            fclose(chunkPolishedRefOutFh);
            free(chunkPolishedRefFilename);
            free(chunkPolishedRefContigName);
        }
    }

    // cleanup
    free(helenFeatureOutfileBase);
    if (trueRefAlignment != NULL) stList_destruct(trueRefAlignment);
    if (trueRefRleString != NULL) rleString_destruct(trueRefRleString);
}

void getDiploidHaplotypeAlignmentsRAW2RLE(RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
                                          RleString *trueRefRleStringA, RleString *trueRefRleStringB,
                                          RleString **trueRefRleStringToHap1, RleString **trueRefRleStringToHap2,
                                          stList **trueRefAlignmentToHap1, stList **trueRefAlignmentToHap2,
                                          Params *params, char *logIdentifier) {

    char *polishedConsensusStringH1 = rleString_expand(polishedRleConsensusH1);
    char *polishedConsensusStringH2 = rleString_expand(polishedRleConsensusH2);
    char *trueRefExpandedA = rleString_expand(trueRefRleStringA);
    char *trueRefExpandedB = rleString_expand(trueRefRleStringB);

    // align all to all
    uint16_t score_trueA_polished1, score_trueA_polished2, score_trueB_polished1, score_trueB_polished2;

    stList *trueRefAlignmentRawSpace_Polished1TrueA = alignConsensusAndTruthSSW(
            polishedConsensusStringH1, trueRefExpandedA, &score_trueA_polished1);

    stList *trueRefAlignmentRawSpace_Polished1TrueB = alignConsensusAndTruthSSW(
            polishedConsensusStringH1, trueRefExpandedB, &score_trueB_polished1);

    stList *trueRefAlignmentRawSpace_Polished2TrueA= alignConsensusAndTruthSSW(
            polishedConsensusStringH2, trueRefExpandedA, &score_trueA_polished2);

    stList *trueRefAlignmentRawSpace_Polished2TrueB = alignConsensusAndTruthSSW(
            polishedConsensusStringH2, trueRefExpandedB, &score_trueB_polished2);

    // determine best alignment
    bool use_trueA_polished1 = score_trueA_polished1 + score_trueB_polished2 > score_trueA_polished2 + score_trueB_polished1;

    // convert to rleSpace if appropriate
    if (params->polishParams->useRunLengthEncoding) {

        // hap1
        *trueRefRleStringToHap1 = rleString_copy(use_trueA_polished1 ? trueRefRleStringA : trueRefRleStringB);

        uint64_t *polishedRleConsensus1_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                polishedRleConsensusH1);
        uint64_t *trueRefRleStringToHap1_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                *trueRefRleStringToHap1);

        *trueRefAlignmentToHap1 = runLengthEncodeAlignment(
                use_trueA_polished1 ? trueRefAlignmentRawSpace_Polished1TrueA : trueRefAlignmentRawSpace_Polished1TrueB,
                polishedRleConsensus1_nonRleToRleCoordinateMap, trueRefRleStringToHap1_nonRleToRleCoordinateMap);

        free(polishedRleConsensus1_nonRleToRleCoordinateMap);
        free(trueRefRleStringToHap1_nonRleToRleCoordinateMap);

        // hap2
        *trueRefRleStringToHap2 = rleString_copy(use_trueA_polished1 ? trueRefRleStringB : trueRefRleStringA);

        uint64_t *polishedRleConsensus2_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                polishedRleConsensusH2);
        uint64_t *trueRefRleStringToHap2_nonRleToRleCoordinateMap = rleString_getNonRleToRleCoordinateMap(
                *trueRefRleStringToHap2);

        *trueRefAlignmentToHap2 = runLengthEncodeAlignment(
                use_trueA_polished1 ? trueRefAlignmentRawSpace_Polished2TrueB : trueRefAlignmentRawSpace_Polished2TrueA,
                polishedRleConsensus2_nonRleToRleCoordinateMap, trueRefRleStringToHap2_nonRleToRleCoordinateMap);

        free(polishedRleConsensus2_nonRleToRleCoordinateMap);
        free(trueRefRleStringToHap2_nonRleToRleCoordinateMap);

        //cleanup
        stList_destruct(trueRefAlignmentRawSpace_Polished1TrueA);
        stList_destruct(trueRefAlignmentRawSpace_Polished1TrueB);
        stList_destruct(trueRefAlignmentRawSpace_Polished2TrueA);
        stList_destruct(trueRefAlignmentRawSpace_Polished2TrueB);
    }

    // debugging
    if (st_getLogLevel() <= debug) {
        st_logInfo(" %s Alignment of truth to Hap1:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH1, *trueRefRleStringToHap1, *trueRefAlignmentToHap1);
        char *consensusRawH1 = rleString_expand(polishedRleConsensusH1);
        char *truthRawToH1 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H1:\n    %s\n", logIdentifier, consensusRawH1);
        st_logInfo(" %s RAW Truth Seq H1:\n    %s\n\n", logIdentifier, truthRawToH1);
        free(consensusRawH1);
        free(truthRawToH1);

        st_logInfo(" %s Alignment of truth to Hap2:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH2, *trueRefRleStringToHap2, *trueRefAlignmentToHap2);
        char *consensusRawH2 = rleString_expand(polishedRleConsensusH2);
        char *truthRawToH2 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H2:\n    %s\n", logIdentifier, consensusRawH2);
        st_logInfo(" %s RAW Truth Seq H2:\n    %s\n\n", logIdentifier, truthRawToH2);
        free(consensusRawH2);
        free(truthRawToH2);
    }
}



void getDiploidHaplotypeAlignmentsRLE(RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
                                      RleString *trueRefRleStringA, RleString *trueRefRleStringB,
                                      RleString **trueRefRleStringToHap1, RleString **trueRefRleStringToHap2,
                                      stList **trueRefAlignmentToHap1, stList **trueRefAlignmentToHap2,
                                      Params *params, char *logIdentifier) {

    // align all to all
    double score_polished1_trueA, score_polished2_trueA, score_polished1_trueB, score_polished2_trueB;

    stList *trueRefAlignmentRLESpace_Polished1TrueA = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH1, trueRefRleStringA, &score_polished1_trueA, params->polishParams);

    stList *trueRefAlignmentRLESpace_Polished1TrueB = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH1, trueRefRleStringB, &score_polished1_trueB, params->polishParams);

    stList *trueRefAlignmentRLESpace_Polished2TrueA = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH2, trueRefRleStringA, &score_polished2_trueA, params->polishParams);

    stList *trueRefAlignmentRLESpace_Polished2TrueB = alignConsensusAndTruthRLEWithKmerAnchors(
            polishedRleConsensusH2, trueRefRleStringB, &score_polished2_trueB, params->polishParams);

    // cis or trans
    bool use_polished1_trueA = score_polished1_trueA + score_polished2_trueB > score_polished2_trueA + score_polished1_trueB;

    // hap1
    *trueRefRleStringToHap1 = rleString_copy(use_polished1_trueA ? trueRefRleStringA : trueRefRleStringB);
    *trueRefAlignmentToHap1 = use_polished1_trueA ?
                              trueRefAlignmentRLESpace_Polished1TrueA : trueRefAlignmentRLESpace_Polished1TrueB;

    // hap2
    *trueRefRleStringToHap2 = rleString_copy(use_polished1_trueA ? trueRefRleStringB : trueRefRleStringA);
    *trueRefAlignmentToHap2 = use_polished1_trueA ?
                              trueRefAlignmentRLESpace_Polished2TrueB : trueRefAlignmentRLESpace_Polished2TrueA;

    // debugging
    if (st_getLogLevel() <= info) {
        st_logInfo(" %s Alignment of truth to Hap1:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH1, *trueRefRleStringToHap1, *trueRefAlignmentToHap1);
        char *consensusRawH1 = rleString_expand(polishedRleConsensusH1);
        char *truthRawToH1 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H1:\n    %s\n", logIdentifier, consensusRawH1);
        st_logInfo(" %s RAW Truth Seq H1:\n    %s\n\n", logIdentifier, truthRawToH1);
        free(consensusRawH1);
        free(truthRawToH1);

        st_logInfo(" %s Alignment of truth to Hap2:\n", logIdentifier);
        printMEAAlignment2(polishedRleConsensusH2, *trueRefRleStringToHap2, *trueRefAlignmentToHap2);
        char *consensusRawH2 = rleString_expand(polishedRleConsensusH2);
        char *truthRawToH2 = rleString_expand(*trueRefRleStringToHap1);
        st_logInfo(" %s RAW Consensus Seq H2:\n    %s\n", logIdentifier, consensusRawH2);
        st_logInfo(" %s RAW Truth Seq H2:\n    %s\n\n", logIdentifier, truthRawToH2);
        free(consensusRawH2);
        free(truthRawToH2);
    }

    //cleanup
    if (use_polished1_trueA) {
        stList_destruct(trueRefAlignmentRLESpace_Polished1TrueB);
        stList_destruct(trueRefAlignmentRLESpace_Polished2TrueA);
    } else {
        stList_destruct(trueRefAlignmentRLESpace_Polished1TrueA);
        stList_destruct(trueRefAlignmentRLESpace_Polished2TrueB);
    }
}

void handleDiploidHelenFeatures(
        // global params
        HelenFeatureType helenFeatureType, BamChunker *trueReferenceBamChunker,
        int64_t splitWeightMaxRunLength, void **helenHDF5Files, bool fullFeatureOutput,
        char *trueReferenceBamA, char *trueReferenceBamB, Params *params,

        // chunk params
        char *logIdentifier, int64_t chunkIdx, BamChunk *bamChunk, stList *bamChunkReads, Poa *poaH1, Poa *poaH2,
        stSet *readsInH1BCR, stSet *readsInH2BCR, RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
        RleString *originalReference) {

    st_logInfo(">%s Performing diploid feature generation for chunk.\n", logIdentifier);

    // get filename
    char *helenFeatureOutfileBase = NULL;
    switch (helenFeatureType) {
        case HFEAT_DIPLOID_RLE_WEIGHT:
            helenFeatureOutfileBase = stString_print("diploidRleWeight.C%05"PRId64".%s-%"PRId64"-%"PRId64,
                                                     chunkIdx, bamChunk->refSeqName,
                                                     bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
            break;
        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    // necessary to annotate poa with truth (if true reference BAM has been specified)
    stList *trueRefAlignmentToHap1 = NULL;
    stList *trueRefAlignmentToHap2 = NULL;
    RleString *trueRefRleStringToHap1 = NULL;
    RleString *trueRefRleStringToHap2 = NULL;
    bool validReferenceAlignment = FALSE;

    // get reference chunk
    if (trueReferenceBamA != NULL && trueReferenceBamB != NULL) {
        // get alignment of true ref to assembly
        stList *trueRefReadsA = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *trueRefReadsB = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *unusedA = stList_construct3(0, (void (*)(void *)) stList_destruct);
        stList *unusedB = stList_construct3(0, (void (*)(void *)) stList_destruct);
        // construct new chunk
        BamChunk *trueRefBamChunk = bamChunk_copyConstruct(bamChunk);
        trueRefBamChunk->parent = trueReferenceBamChunker;
        // get true ref as "read"
        uint32_t trueAlignmentCountA = convertToReadsAndAlignments(trueRefBamChunk, originalReference, trueRefReadsA, unusedA);
        uint32_t trueAlignmentCountB = convertToReadsAndAlignments(trueRefBamChunk, originalReference, trueRefReadsB, unusedB);

        // poor man's "do we have a unique alignment"
        if (trueAlignmentCountA != 1 || trueAlignmentCountB != 1) {
            st_logInfo(" %s Found %"PRId64", %"PRId64" true reference alignments for feature generation.\n",
                       logIdentifier, trueAlignmentCountA, trueAlignmentCountB);
            //todo maybe wrap in a if loglevel <= debug?
            stList *alignmentDesc = stList_construct3(0, free);
            for (int i = 0; i < stList_length(unusedA); i++) {
                stList *truthAlign = stList_get(unusedA, i);
                stList_append(alignmentDesc, stList_length(truthAlign) == 0 ? stString_copy("A,0-0,0") :
                                             stString_print("A,%"PRId64"-%"PRId64",%"PRId64,
                                                            stIntTuple_get(stList_get(truthAlign,0), 0),
                                                            stIntTuple_get(stList_get(truthAlign, stList_length(truthAlign) - 1), 0),
                                                            ((BamChunkRead*)stList_get(trueRefReadsA, i))->rleRead->length));
            }
            for (int i = 0; i < stList_length(unusedB); i++) {
                stList *truthAlign = stList_get(unusedB, i);
                stList_append(alignmentDesc, stList_length(truthAlign) == 0 ? stString_copy("B,0-0,0") :
                                             stString_print("B,%"PRId64"-%"PRId64",%"PRId64,
                                                            stIntTuple_get(stList_get(truthAlign,0), 0),
                                                            stIntTuple_get(stList_get(truthAlign, stList_length(truthAlign) - 1), 0),
                                                            ((BamChunkRead*)stList_get(trueRefReadsB, i))->rleRead->length));
            }
            char *alnSum = stString_join2("; ", alignmentDesc);
            st_logInfo(" %s   Alignment summary: %s\n", logIdentifier, alnSum);
            free(alnSum);
            stList_destruct(alignmentDesc);

        } else {
            // get reads
            BamChunkRead *trueRefReadA = stList_get(trueRefReadsA, 0);
            BamChunkRead *trueRefReadB = stList_get(trueRefReadsB, 0);

            RleString *trueRefRleStringA = trueRefReadA->rleRead;
            RleString *trueRefRleStringB = trueRefReadB->rleRead;

            // assign correct haplotypes to the trueRefXXXToHapX
            getDiploidHaplotypeAlignmentsRLE(polishedRleConsensusH1, polishedRleConsensusH2, trueRefRleStringA, trueRefRleStringB,
                    &trueRefRleStringToHap1, &trueRefRleStringToHap2, &trueRefAlignmentToHap1, &trueRefAlignmentToHap2,
                    params, logIdentifier);

            // we found a single alignment of reference
            double refLengthRatio = 1.0 * (trueRefRleStringA->length + trueRefRleStringB->length) /
                    (polishedRleConsensusH1->length + polishedRleConsensusH2->length);
            double alnLengthRatio = 1.0 * (stList_length(trueRefAlignmentToHap1) + stList_length(trueRefAlignmentToHap2)) /
                    (polishedRleConsensusH1->length + polishedRleConsensusH2->length);
            int refLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - refLengthRatio)));
            int alnLengthRatioHundredthsOffOne = abs((int) (100 * (1.0 - alnLengthRatio)));
            if (stList_length(trueRefAlignmentToHap1) > 0 && stList_length(trueRefAlignmentToHap2) > 0 &&
                    refLengthRatioHundredthsOffOne < 10 && alnLengthRatioHundredthsOffOne < 10) {
                validReferenceAlignment = TRUE;
            } else {
                st_logInfo(" %s True reference alignment QC failed:  avg polished length %"PRId64", true ref length"
                           " ratio (true/polished) %f, aligned pairs length ratio (true/polished): %f\n",
                           logIdentifier, (polishedRleConsensusH1->length + polishedRleConsensusH2->length) / 2,
                           refLengthRatio, alnLengthRatio);
            }
        }

        stList_destruct(trueRefReadsA);
        stList_destruct(unusedA);
        stList_destruct(trueRefReadsB);
        stList_destruct(unusedB);
        bamChunk_destruct(trueRefBamChunk);
    }

    // either write it, or note that we failed to find a valid reference alignment
    if (trueReferenceBamA != NULL && trueReferenceBamB != NULL && !validReferenceAlignment) {
        st_logInfo(" %s No valid reference alignment was found, skipping HELEN feature output.\n", logIdentifier);
    } else {
        st_logInfo(" %s Writing HELEN features with filename base: %s\n", logIdentifier, helenFeatureOutfileBase);

        // write the actual features (type dependent)
        stSet *readIdsInH1 = bamChunkRead_to_readName(readsInH1BCR);
        stSet *readIdsInH2 = bamChunkRead_to_readName(readsInH2BCR);
        poa_writeDiploidHelenFeatures(helenFeatureType, bamChunkReads, helenFeatureOutfileBase,
                bamChunk, poaH1, poaH2, readIdsInH1, readIdsInH2, trueRefAlignmentToHap1, trueRefAlignmentToHap2,
                trueRefRleStringToHap1,  trueRefRleStringToHap2, splitWeightMaxRunLength,
                (HelenFeatureHDF5FileInfo**) helenHDF5Files);
        stSet_destruct(readIdsInH1);
        stSet_destruct(readIdsInH2);
    }

    // cleanup
    free(helenFeatureOutfileBase);
    if (trueRefAlignmentToHap1 != NULL) stList_destruct(trueRefAlignmentToHap1);
    if (trueRefAlignmentToHap1 != NULL) stList_destruct(trueRefAlignmentToHap2);
    if (trueRefRleStringToHap1 != NULL) rleString_destruct(trueRefRleStringToHap1);
    if (trueRefRleStringToHap2 != NULL) rleString_destruct(trueRefRleStringToHap2);
}

stList *poa_getSimpleWeightFeatures(Poa *poa, stList *bamChunkReads) {

    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_SimpleWeight_destruct);
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_SimpleWeight_construct(i - 1, 0));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for (int64_t i = 0; i < stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureSimpleWeight *feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i +
                                               1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // examine each observation
        stList *observations = node->observations;
        for (int64_t o = 0; o < stList_length(observations); o++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, o);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            Symbol symbol = poa->alphabet->convertCharToSymbol(bamChunkRead->rleRead->rleString[observation->offset]);
            bool forward = bamChunkRead->forwardStrand;

            // save weight
            feature->weights[PoaFeature_SimpleWeight_charIndex(symbol, forward)] += observation->weight;
        }

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t j = 0; j < stList_length(node->deletes); j++) {
                PoaDelete *delete = stList_get(node->deletes, j);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureSimpleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->weights[PoaFeature_SimpleWeight_gapIndex(TRUE)] += delete->weightForwardStrand;
                    delFeature->weights[PoaFeature_SimpleWeight_gapIndex(FALSE)] += delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t j = 0; j < stList_length(node->inserts); j++) {
                PoaInsert *insert = stList_get(node->inserts, j);

                // get feature iterator
                PoaFeatureSimpleWeight *prevFeature = feature;
                for (int64_t k = 0; k < strlen(insert->insert->rleString); k++) {
                    // get current feature (or create if necessary)
                    PoaFeatureSimpleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_SimpleWeight_construct(i, k + 1);
                        prevFeature->nextInsert = currFeature;
                    }

                    Symbol c = poa->alphabet->convertCharToSymbol(insert->insert->rleString[k]);

                    // add weights
                    currFeature->weights[PoaFeature_SimpleWeight_charIndex(c, TRUE)] += insert->weightForwardStrand;
                    currFeature->weights[PoaFeature_SimpleWeight_charIndex(c, FALSE)] += insert->weightReverseStrand;

                    // iterate
                    prevFeature = currFeature;
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}


void poa_addSplitRunLengthFeaturesForObservations(Poa *poa, PoaFeatureSplitRleWeight *baseFeature, stList *observations,
                                                  stList *bamChunkReads, const int64_t maxRunLength,
                                                  int64_t observationOffset) {


    PoaFeatureSplitRleWeight *currFeature = baseFeature;
    int64_t currentRunLengthIndex = 0;
    bool beforeMaxObservedRunLength = TRUE;

    while (beforeMaxObservedRunLength) {
        beforeMaxObservedRunLength = FALSE;

        // examine each observation
        for (int64_t i = 0; i < stList_length(observations); i++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            RleString *rleString = bamChunkRead->rleRead;
            Symbol symbol = poa->alphabet->convertCharToSymbol(
                    rleString->rleString[observation->offset + observationOffset]);
            int64_t runLength = rleString->repeatCounts[observation->offset + observationOffset];
            bool forward = bamChunkRead->forwardStrand;

            // get correct run length
            runLength -= currentRunLengthIndex * maxRunLength;
            if (runLength < 0) {
                runLength = 0;
            } else if (runLength > maxRunLength) {
                runLength = maxRunLength;
                beforeMaxObservedRunLength = TRUE;
            }

            int64_t pos = PoaFeature_SplitRleWeight_charIndex(maxRunLength, symbol, runLength, forward);
            currFeature->weights[pos] += observation->weight;

        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
                PoaFeatureSplitRleWeight *prevFeature = currFeature;
                currFeature = PoaFeature_SplitRleWeight_construct(baseFeature->refPosition, baseFeature->insertPosition,
                                                                  currentRunLengthIndex, maxRunLength);
                prevFeature->nextRunLength = currFeature;
                currFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE)] =
                        baseFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE)];
                currFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE)] =
                        baseFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE)];
            }
        }
    }
}


stList *poa_getSplitRleWeightFeatures(Poa *poa, stList *bamChunkReads, const int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_SplitRleWeight_destruct);
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_SplitRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for (int64_t i = 0; i < stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureSplitRleWeight *feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i +
                                               1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addSplitRunLengthFeaturesForObservations(poa, feature, node->observations, bamChunkReads,
                                                     maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureSplitRleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength,
                                                                           TRUE)] += delete->weightForwardStrand;
                    delFeature->weights[PoaFeature_SplitRleWeight_gapIndex(maxRunLength,
                                                                           FALSE)] += delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureSplitRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert->rleString); o++) {

                    // get feature iterator
                    PoaFeatureSplitRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_SplitRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addSplitRunLengthFeaturesForObservations(poa, currFeature, insert->observations, bamChunkReads,
                                                                 maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}

void poa_addChannelRunLengthFeaturesForObservations(Poa *poa, PoaFeatureChannelRleWeight *baseFeature, stList *observations,
                                                    stList *bamChunkReads, const int64_t maxRunLength,
                                                    int64_t observationOffset) {

    PoaFeatureChannelRleWeight *currFeature = baseFeature;
    int64_t currentRunLengthIndex = 0;
    bool beforeMaxObservedRunLength = TRUE;

    while (beforeMaxObservedRunLength) {
        beforeMaxObservedRunLength = FALSE;

        // examine each observation
        for (int64_t i = 0; i < stList_length(observations); i++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            RleString *rleString = bamChunkRead->rleRead;
            Symbol symbol = poa->alphabet->convertCharToSymbol(
                    rleString->rleString[observation->offset + observationOffset]);
            int64_t runLength = rleString->repeatCounts[observation->offset + observationOffset];
            bool forward = bamChunkRead->forwardStrand;

            // get correct run length
            runLength -= currentRunLengthIndex * maxRunLength;
            if (runLength < 0) {
                runLength = 0;
            } else if (runLength > maxRunLength) {
                runLength = maxRunLength;
                beforeMaxObservedRunLength = TRUE;
            }

            // save weight for both nucleotide totals and runLenght totals
            int64_t nuclPos = PoaFeature_ChannelRleWeight_charNuclIndex(symbol, forward);
            int64_t runLengthPos = PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, symbol, runLength, forward);
            currFeature->nucleotideWeights[nuclPos] += observation->weight;
            currFeature->runLengthWeights[runLengthPos] += observation->weight;

        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
                PoaFeatureChannelRleWeight *prevFeature = currFeature;
                currFeature = PoaFeature_ChannelRleWeight_construct(baseFeature->refPosition,
                                                                    baseFeature->insertPosition,
                                                                    currentRunLengthIndex, maxRunLength);
                prevFeature->nextRunLength = currFeature;
                currFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)] =
                        baseFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)];
                currFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)] =
                        baseFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)];
            }
        }
    }
}

stList *poa_getChannelRleWeightFeatures(Poa *poa, stList *bamChunkReads, int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_ChannelRleWeight_destruct);
    for (int64_t i = 1; i < stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_ChannelRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for (int64_t i = 0; i < stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureChannelRleWeight *feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i +
                                               1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addChannelRunLengthFeaturesForObservations(poa, feature, node->observations, bamChunkReads,
                                                       maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logCritical(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureChannelRleWeight *delFeature = stList_get(featureList, i + k);

                    delFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)] +=
                            delete->weightForwardStrand;
                    delFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)] +=
                            delete->weightReverseStrand;
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureChannelRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert->rleString); o++) {

                    // get feature iterator
                    PoaFeatureChannelRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_ChannelRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addChannelRunLengthFeaturesForObservations(poa, currFeature, insert->observations,
                                                                   bamChunkReads,
                                                                   maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}


void poa_addDiploidRunLengthFeaturesForObservations(Poa *poa, PoaFeatureDiploidRleWeight *baseFeature, stList *observations,
                                                  stList *bamChunkReads, stSet *onHapReads, const int64_t maxRunLength,
                                                  int64_t observationOffset) {


    PoaFeatureDiploidRleWeight *currFeature = baseFeature;
    int64_t currentRunLengthIndex = 0;
    bool beforeMaxObservedRunLength = TRUE;

    while (beforeMaxObservedRunLength) {
        beforeMaxObservedRunLength = FALSE;

        // examine each observation
        for (int64_t i = 0; i < stList_length(observations); i++) {

            // get data
            PoaBaseObservation *observation = stList_get(observations, i);
            BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
            RleString *rleString = bamChunkRead->rleRead;
            Symbol symbol = poa->alphabet->convertCharToSymbol(rleString->rleString[observation->offset + observationOffset]);
            int64_t runLength = rleString->repeatCounts[observation->offset + observationOffset];
            bool forward = bamChunkRead->forwardStrand;

            // get correct run length
            runLength -= currentRunLengthIndex * maxRunLength;
            if (runLength < 0) {
                runLength = 0;
            } else if (runLength > maxRunLength) {
                runLength = maxRunLength;
                beforeMaxObservedRunLength = TRUE;
            }

            int64_t pos = PoaFeature_DiploidRleWeight_charIndex(maxRunLength, symbol, runLength, forward);
            if (stSet_search(onHapReads, bamChunkRead->readName) != NULL) {
                currFeature->weightsHOn[pos] += observation->weight;
            } else {
                currFeature->weightsHOff[pos] += observation->weight;
            }
        }

        // update currFeature if we're going ot run again
        if (beforeMaxObservedRunLength) {
            currentRunLengthIndex++;
            if (currFeature->nextRunLength != NULL) {
                currFeature = currFeature->nextRunLength;
            } else {
                PoaFeatureDiploidRleWeight *prevFeature = currFeature;
                currFeature = PoaFeature_DiploidRleWeight_construct(baseFeature->refPosition, baseFeature->insertPosition,
                                                                  currentRunLengthIndex, maxRunLength);
                prevFeature->nextRunLength = currFeature;

                currFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)] =
                        baseFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)];
                currFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)] =
                        baseFeature->weightsHOn[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)];

                currFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)] =
                        baseFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE)];
                currFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)] =
                        baseFeature->weightsHOff[PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE)];
            }
        }
    }
}



stList *poa_getDiploidRleWeightFeatures(Poa *poa, stList *bamChunkReads, stSet *onHapReads, const int64_t maxRunLength) {
    // initialize feature list
    stList *featureList = stList_construct3(0, (void (*)(void *)) PoaFeature_DiploidRleWeight_destruct);
    for(int64_t i=1; i<stList_length(poa->nodes); i++) {
        stList_append(featureList, PoaFeature_DiploidRleWeight_construct(i - 1, 0, 0, maxRunLength));
    }

    // for logging (of errors)
    char *logIdentifier = getLogIdentifier();

    // iterate over all positions
    for(int64_t i=0; i<stList_length(featureList); i++) {

        // get feature and node
        PoaFeatureDiploidRleWeight* feature = stList_get(featureList, i);
        PoaNode *node = stList_get(poa->nodes, i + 1); //skip the first poa node, as it's always an 'N', so featureIdx and poaIdx are off by one

        // save run length nodes
        poa_addDiploidRunLengthFeaturesForObservations(poa, feature, node->observations, bamChunkReads, onHapReads,
                                                       maxRunLength, 0);

        // Deletes
        if (stList_length(node->deletes) > 0) {

            // iterate over all deletes
            for (int64_t d = 0; d < stList_length(node->deletes); d++) {
                PoaDelete *delete = stList_get(node->deletes, d);

                // Deletes start AFTER the current position, need to add counts/weights to nodes after the current node
                for (int64_t k = 1; k < delete->length; k++) {
                    if (i + k >= stList_length(poa->nodes)) {
                        st_logInfo(" %s Encountered DELETE that occurs after end of POA!\n", logIdentifier);
                        break;
                    }
                    PoaFeatureDiploidRleWeight *delFeature = stList_get(featureList, i + k);
                    stList *observations = delete->observations;

                    // examine each observation for hap and direction
                    for (int64_t o = 0; o < stList_length(observations); o++) {

                        // get observation data
                        PoaBaseObservation *observation = stList_get(observations, o);
                        BamChunkRead *bamChunkRead = stList_get(bamChunkReads, observation->readNo);
                        bool forward = bamChunkRead->forwardStrand;

                        int64_t pos = PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, forward);
                        if (stSet_search(onHapReads, bamChunkRead->readName) != NULL) {
                            delFeature->weightsHOn[pos] += observation->weight;
                        } else {
                            delFeature->weightsHOff[pos] += observation->weight;
                        }
                    }
                }
            }
        }

        // Inserts
        if (stList_length(node->inserts) > 0) {

            // iterate over all inserts
            for (int64_t n = 0; n < stList_length(node->inserts); n++) {
                PoaInsert *insert = stList_get(node->inserts, n);

                // handle each insert base
                PoaFeatureDiploidRleWeight *prevFeature = feature;
                for (int64_t o = 0; o < strlen(insert->insert->rleString); o++) {

                    // get feature iterator
                    PoaFeatureDiploidRleWeight *currFeature = prevFeature->nextInsert;
                    if (currFeature == NULL) {
                        currFeature = PoaFeature_DiploidRleWeight_construct(i, o + 1, 0, maxRunLength);
                        prevFeature->nextInsert = currFeature;
                    }

                    // save insert run lengths
                    poa_addDiploidRunLengthFeaturesForObservations(poa, currFeature, insert->observations,
                            bamChunkReads, onHapReads, maxRunLength, o);
                }
            }
        }
    }

    free(logIdentifier);
    return featureList;
}





void printMEAAlignment2(RleString *X, RleString *Y, stList *alignedPairs) {
    printMEAAlignment(X->rleString, Y->rleString, X->length, Y->length, alignedPairs, X->repeatCounts, Y->repeatCounts);
}

void printMEAAlignment(char *X, char *Y, int64_t lX, int64_t lY, stList *alignedPairs, uint64_t *Xrl, uint64_t *Yrl) {
    if (stList_length(alignedPairs) == 0) {
        return;
    }

    // should we do run lengths
    bool handleRunLength = Xrl != NULL && Yrl != NULL;

    // build strings to print
    int64_t bufferLen = lX + lY;
    char *alnXStr = st_calloc(bufferLen, sizeof(char));
    char *alnYStr = st_calloc(bufferLen, sizeof(char));
    char *alnDesc = st_calloc(bufferLen, sizeof(char));
    char *rlXStr = NULL;
    char *rlYStr = NULL;
    if (handleRunLength) {
        rlXStr = st_calloc(bufferLen, sizeof(char));
        rlYStr = st_calloc(bufferLen, sizeof(char));
    }

    // stats to track
    int64_t nuclMatches = 0;
    int64_t nuclMismatches = 0;
    int64_t nuclXInserts = 0;
    int64_t nuclYInserts = 0;
    int64_t rlMatches = 0;
    int64_t rlMismatches = 0;


    // iterate over alignment
    stListIterator *alignmentItor = stList_getIterator(alignedPairs);
    stIntTuple *prevAlign = NULL;
    stIntTuple *currAlign = stList_getNext(alignmentItor);
    int64_t posX = stIntTuple_get(currAlign, 0);
    int64_t posY = stIntTuple_get(currAlign, 1);
    int64_t outStrPos = 0;

    while (TRUE) {
        if (currAlign == NULL) break;
        int64_t currAlignPosX = stIntTuple_get(currAlign, 0);
        int64_t currAlignPosY = stIntTuple_get(currAlign, 1);
        // Y gap / X insert
        if (posX < currAlignPosX) {
            alnXStr[outStrPos] = X[posX];
            alnYStr[outStrPos] = '_';
            alnDesc[outStrPos] = ' ';
            if (handleRunLength) {
                rlXStr[outStrPos] = (char) ('0' + Xrl[posX]);
                rlYStr[outStrPos] = ' ';
            }
            posX++;
            nuclXInserts++;
        }

            // X gap / Y insert
        else if (posY < currAlignPosY) {
            alnXStr[outStrPos] = '_';
            alnYStr[outStrPos] = Y[posY];
            alnDesc[outStrPos] = ' ';
            if (handleRunLength) {
                rlXStr[outStrPos] = ' ';
                rlYStr[outStrPos] = (char) ('0' + Yrl[posY]);
            }
            posY++;
            nuclYInserts++;
        }

            // match
        else if (posX == currAlignPosX && posY == currAlignPosY) {
            alnXStr[outStrPos] = X[posX];
            alnYStr[outStrPos] = Y[posY];
            if (handleRunLength) {
                rlXStr[outStrPos] = (char) ('0' + Xrl[posX]);
                rlYStr[outStrPos] = (char) ('0' + Yrl[posY]);
            }
            if (X[posX] == Y[posY]) {
                nuclMatches++;
                alnDesc[outStrPos] = '|';
                if (handleRunLength) {
                    if (Xrl[posX] == Yrl[posY]) {
                        rlMatches++;
                    } else {
                        rlMismatches++;
                        alnDesc[outStrPos] = ':';
                    }
                }
            } else {
                alnDesc[outStrPos] = ' ';
                nuclMismatches++;
            }
            posX++;
            posY++;
            prevAlign = currAlign;
            currAlign = stList_getNext(alignmentItor);
        }

            // should never happen
        else {
            assert(FALSE);
        }

        outStrPos++;
    }

    // print
    fprintf(stderr, "\n");
    if (handleRunLength) fprintf(stderr, "  %s\n", rlXStr);
    fprintf(stderr, "  %s\n", alnXStr);
    fprintf(stderr, "  %s\n", alnDesc);
    fprintf(stderr, "  %s\n", alnYStr);
    if (handleRunLength) fprintf(stderr, "  %s\n", rlYStr);
    fprintf(stderr, "  Matches:    %"PRId64"\n", nuclMatches);
    if (handleRunLength) {
        fprintf(stderr, "    RL Match: %"PRId64"\n", rlMatches);
        fprintf(stderr, "    RL Miss:  %"PRId64"\n", rlMismatches);
    }
    fprintf(stderr, "  Mismatches: %"PRId64"\n", nuclMismatches);
    fprintf(stderr, "  X Inserts:  %"PRId64"\n", nuclXInserts);
    fprintf(stderr, "  Y Inserts:  %"PRId64"\n", nuclYInserts);
    fprintf(stderr, "\n");

    // cleanup
    free(alnXStr);
    free(alnYStr);
    free(alnDesc);
    if (handleRunLength) {
        free(rlXStr);
        free(rlYStr);
    }
    stList_destructIterator(alignmentItor);
}


void poa_annotateHelenFeaturesWithTruth(stList *features, HelenFeatureType featureType, stList *trueRefAlignment,
                                        RleString *trueRefRleString, int64_t *firstMatchedFeaure,
                                        int64_t *lastMatchedFeature) {
    /*
     Each index in features represents a character in the final consensus string (which in turn was a single node in the
     final POA which was used to generate the consensus).  This consensus was aligned against the true reference (which
     is reflected in the trueRefRleString).  Items in the true refAlignment are stIntTuples with index 0 being the
     alignment weight (discarded), index 1 being the position in the consensus string (and features), and index 2 being
     the position in the trueRefRleString.  So we can iterate over features and the true alignment to assign truth
     labels to each feature.
     */
    static int FEATURE_POS = 0;
    static int REFERENCE_POS = 1;
    *firstMatchedFeaure = -1;
    *lastMatchedFeature = -1;
    char *logIdentifier = getLogIdentifier();

    // iterate over true ref alignment
    stListIterator *trueRefAlignItor = stList_getIterator(trueRefAlignment);
    stIntTuple *currTrueRefAlign = stList_getNext(trueRefAlignItor);

    // iterate over features
    int64_t trueRefPos = stIntTuple_get(currTrueRefAlign, REFERENCE_POS);
    for (int64_t featureRefPos = 0; featureRefPos < stList_length(features); featureRefPos++) {
        void *feature = stList_get(features, featureRefPos);
        int64_t trueRunLength = -1;
        PoaFeatureSplitRleWeight* srlFeature = NULL;
        PoaFeatureChannelRleWeight* crlFeature = NULL;
        PoaFeatureDiploidRleWeight* drlFeature = NULL;

        int64_t featureInsPos = 0;
        while (feature != NULL) {

            // no more ref bases, everything is gaps
            if (currTrueRefAlign == NULL) {
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = '_';
                        feature = ((PoaFeatureSimpleWeight *) feature)->nextInsert;
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = '_';
                            srlFeature->labelRunLength = 0;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureSplitRleWeight *) feature)->nextInsert;
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = '_';
                            crlFeature->labelRunLength = 0;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureChannelRleWeight *) feature)->nextInsert;
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = '_';
                            drlFeature->labelRunLength = 0;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        feature = ((PoaFeatureDiploidRleWeight*)feature)->nextInsert;
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                continue;
            }

            // sanity checks
            assert(stIntTuple_get(currTrueRefAlign, FEATURE_POS) >= featureRefPos &&
                   stIntTuple_get(currTrueRefAlign, REFERENCE_POS) >= trueRefPos);

            // match
            if (stIntTuple_get(currTrueRefAlign, FEATURE_POS) == featureRefPos &&
                stIntTuple_get(currTrueRefAlign, REFERENCE_POS) == trueRefPos) {
                st_logDebug(
                        " %s LABEL MATCH  %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                        logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // save label (based on feature type)
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                srlFeature->labelRunLength = 0;
                            } else if (trueRunLength > srlFeature->maxRunLength) {
                                srlFeature->labelRunLength = srlFeature->maxRunLength;
                            } else {
                                srlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= srlFeature->maxRunLength;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                crlFeature->labelRunLength = 0;
                            } else if (trueRunLength > crlFeature->maxRunLength) {
                                crlFeature->labelRunLength = crlFeature->maxRunLength;
                            } else {
                                crlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= crlFeature->maxRunLength;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                drlFeature->labelRunLength = 0;
                            } else if (trueRunLength > drlFeature->maxRunLength) {
                                drlFeature->labelRunLength = drlFeature->maxRunLength;
                            } else {
                                drlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= drlFeature->maxRunLength;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }

                // iterate
                trueRefPos++;
                currTrueRefAlign = stList_getNext(trueRefAlignItor);
                // handle first and last match
                if (featureInsPos == 0) {
                    if (*firstMatchedFeaure == -1) {
                        *firstMatchedFeaure = featureRefPos;
                    }
                    *lastMatchedFeature = featureRefPos;
                }
            }

                // insert
            else if (trueRefPos < stIntTuple_get(currTrueRefAlign, REFERENCE_POS)) {
                st_logDebug(
                        " %s LABEL INSERT %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                        logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = trueRefRleString->rleString[trueRefPos];
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                srlFeature->labelRunLength = 0;
                            } else if (trueRunLength > srlFeature->maxRunLength) {
                                srlFeature->labelRunLength = srlFeature->maxRunLength;
                            } else {
                                srlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= srlFeature->maxRunLength;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                crlFeature->labelRunLength = 0;
                            } else if (trueRunLength > crlFeature->maxRunLength) {
                                crlFeature->labelRunLength = crlFeature->maxRunLength;
                            } else {
                                crlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= crlFeature->maxRunLength;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        trueRunLength = trueRefRleString->repeatCounts[trueRefPos];
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = trueRefRleString->rleString[trueRefPos];
                            if (trueRunLength <= 0) {
                                drlFeature->labelRunLength = 0;
                            } else if (trueRunLength > drlFeature->maxRunLength) {
                                drlFeature->labelRunLength = drlFeature->maxRunLength;
                            } else {
                                drlFeature->labelRunLength = trueRunLength;
                            }
                            trueRunLength -= drlFeature->maxRunLength;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
                trueRefPos++;
            }

                // delete
            else if (featureRefPos < stIntTuple_get(currTrueRefAlign, FEATURE_POS)) {
                st_logDebug(
                        " %s LABEL DELETE %c trueRefPos:%"PRId64" featureRefPos:%"PRId64" featureInsPos:%"PRId64"\n",
                        logIdentifier, featureInsPos == 0 ? ' ' : 'I', trueRefPos, featureRefPos, featureInsPos);
                // apply label
                switch (featureType) {
                    case HFEAT_SIMPLE_WEIGHT:
                        ((PoaFeatureSimpleWeight *) feature)->label = '_';
                        break;
                    case HFEAT_SPLIT_RLE_WEIGHT:
                        srlFeature = ((PoaFeatureSplitRleWeight *) feature);
                        while (srlFeature != NULL) {
                            srlFeature->labelChar = '_';
                            srlFeature->labelRunLength = 0;
                            srlFeature = srlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_CHANNEL_RLE_WEIGHT:
                        crlFeature = ((PoaFeatureChannelRleWeight *) feature);
                        while (crlFeature != NULL) {
                            crlFeature->labelChar = '_';
                            crlFeature->labelRunLength = 0;
                            crlFeature = crlFeature->nextRunLength;
                        }
                        break;
                    case HFEAT_DIPLOID_RLE_WEIGHT:
                        drlFeature = ((PoaFeatureDiploidRleWeight*)feature);
                        while (drlFeature != NULL) {
                            drlFeature->labelChar = '_';
                            drlFeature->labelRunLength = 0;
                            drlFeature = drlFeature->nextRunLength;
                        }
                        break;
                    default:
                        st_errAbort("Unhandled FeatureType!\n");
                }
            }

                // programmer error
            else {
                st_errAbort("Unhandled case annotating features with true reference characters!\n");
            }

            // always iterate over insert features
            switch (featureType) {
                case HFEAT_SIMPLE_WEIGHT:
                    feature = ((PoaFeatureSimpleWeight *) feature)->nextInsert;
                    break;
                case HFEAT_SPLIT_RLE_WEIGHT:
                    feature = ((PoaFeatureSplitRleWeight *) feature)->nextInsert;
                    break;
                case HFEAT_CHANNEL_RLE_WEIGHT:
                    feature = ((PoaFeatureChannelRleWeight *) feature)->nextInsert;
                    break;
                case HFEAT_DIPLOID_RLE_WEIGHT:
                    feature = ((PoaFeatureDiploidRleWeight*)feature)->nextInsert;
                    break;
                default:
                    st_errAbort("Unhandled FeatureType!\n");
            }
            featureInsPos++;
        }

        // this catches any true inserts which are not present in the poa / feature list
        while (currTrueRefAlign != NULL &&
               featureRefPos < stIntTuple_get(currTrueRefAlign, FEATURE_POS) &&
               trueRefPos < stIntTuple_get(currTrueRefAlign, REFERENCE_POS)) {
            trueRefPos++;
        }
    }

    stList_destructIterator(trueRefAlignItor);
    free(logIdentifier);
}

void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads,
                            char *outputFileBase, BamChunk *bamChunk, stList *trueRefAlignment, RleString *consensusRleString,
                            RleString *trueRefRleString, bool fullFeatureOutput, int64_t maxRunLength,
                            HelenFeatureHDF5FileInfo** helenHDF5Files) {
    // prep
    int64_t firstMatchedFeature = -1;
    int64_t lastMatchedFeature = -1;
    stList *features = NULL;
    bool outputLabels = trueRefAlignment != NULL && trueRefRleString != NULL;

# ifdef _OPENMP
    int64_t threadIdx = omp_get_thread_num();
# else
    int64_t threadIdx = 0;
# endif

    // handle differently based on type
    switch (type) {
        case HFEAT_SIMPLE_WEIGHT :
            // get features
            features = poa_getSimpleWeightFeatures(poa, bamChunkReads);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                                   &firstMatchedFeature, &lastMatchedFeature);
            }

            writeSimpleWeightHelenFeaturesHDF5(poa->alphabet, helenHDF5Files[threadIdx], outputFileBase, bamChunk,
                                               outputLabels, features, firstMatchedFeature, lastMatchedFeature);

            break;

        case HFEAT_SPLIT_RLE_WEIGHT:
            // get features
            features = poa_getSplitRleWeightFeatures(poa, bamChunkReads, maxRunLength);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                                   &firstMatchedFeature, &lastMatchedFeature);
            }

            writeSplitRleWeightHelenFeaturesHDF5(poa->alphabet, helenHDF5Files[threadIdx],
                                                 outputFileBase, bamChunk, outputLabels, features, firstMatchedFeature, lastMatchedFeature,
                                                 maxRunLength);
            break;

        case HFEAT_CHANNEL_RLE_WEIGHT:
            // get features
            features = poa_getChannelRleWeightFeatures(poa, bamChunkReads, maxRunLength);
            firstMatchedFeature = 0;
            lastMatchedFeature = stList_length(features) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(features, type, trueRefAlignment, trueRefRleString,
                                                   &firstMatchedFeature, &lastMatchedFeature);
            }

            writeChannelRleWeightHelenFeaturesHDF5(poa->alphabet, helenHDF5Files[threadIdx],
                                                   outputFileBase, bamChunk, outputLabels, features, firstMatchedFeature,
                                                   lastMatchedFeature, maxRunLength);
            break;

        default:
            st_errAbort("Unhandled HELEN feature type!\n");
    }

    //cleanup
    stList_destruct(features);
}




void poa_writeDiploidHelenFeatures(HelenFeatureType type, stList *bamChunkReads, char *outputFileBase,
        BamChunk *bamChunk, Poa *poaH1, Poa *poaH2,  stSet *readsInH1, stSet *readsInH2,
        stList *trueRefAlignmentToH1, stList *trueRefAlignmentToH2,
        RleString *trueRefRleStringToH1, RleString *trueRefRleStringToH2,
        int64_t maxRunLength, HelenFeatureHDF5FileInfo** helenHDF5Files) {

    // prep
    int64_t firstMatchedFeatureH1 = -1;
    int64_t firstMatchedFeatureH2 = -1;
    int64_t lastMatchedFeatureH1 = -1;
    int64_t lastMatchedFeatureH2 = -1;
    stList *featuresH1 = NULL;
    stList *featuresH2 = NULL;
    bool outputLabels = trueRefAlignmentToH1 != NULL && trueRefAlignmentToH2 != NULL &&
            trueRefRleStringToH1 != NULL && trueRefRleStringToH1 != NULL;

# ifdef _OPENMP
    int64_t threadIdx = omp_get_thread_num();
# else
    int64_t threadIdx = 0;
# endif

    // handle differently based on type
    switch (type) {
        case HFEAT_DIPLOID_RLE_WEIGHT :
            // get features
            featuresH1 = poa_getDiploidRleWeightFeatures(poaH1, bamChunkReads, readsInH1, maxRunLength);
            featuresH2 = poa_getDiploidRleWeightFeatures(poaH2, bamChunkReads, readsInH2, maxRunLength);

            firstMatchedFeatureH1 = 0;
            firstMatchedFeatureH2 = 0;
            lastMatchedFeatureH1 = stList_length(featuresH1) - 1;
            lastMatchedFeatureH2 = stList_length(featuresH2) - 1;

            // get truth (if we have it)
            if (outputLabels) {
                poa_annotateHelenFeaturesWithTruth(featuresH1, type, trueRefAlignmentToH1, trueRefRleStringToH1,
                                                   &firstMatchedFeatureH1, &lastMatchedFeatureH1);
                poa_annotateHelenFeaturesWithTruth(featuresH2, type, trueRefAlignmentToH2, trueRefRleStringToH2,
                                                   &firstMatchedFeatureH2, &lastMatchedFeatureH2);
            }

            writeDiploidRleWeightHelenFeaturesHDF5(poaH1->alphabet, helenHDF5Files[threadIdx], outputFileBase, bamChunk,
                                                   outputLabels, featuresH1, firstMatchedFeatureH1, lastMatchedFeatureH1,
                                                   featuresH2, firstMatchedFeatureH2, lastMatchedFeatureH2, maxRunLength);

            break;

        default:
            st_errAbort("Unhandled HELEN feature type for diploid!\n");
    }

    //cleanup
    stList_destruct(featuresH1);
    stList_destruct(featuresH2);
}



stList *alignConsensusAndTruthRLE(RleString *consensusStr, RleString *truthStr, double *score, PolishParams *polishParams) {


    // Symbol strings
    SymbolString sX = rleString_constructSymbolString(consensusStr, 0, consensusStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    SymbolString sY = rleString_constructSymbolString(truthStr, 0, truthStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    uint16_t apScore = 0;

    // Run the alignment
    stList *alignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapXPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapYPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);

    getAlignedPairsWithIndels(polishParams->stateMachineForForwardStrandRead, sX, sY,
            polishParams->p, &alignedPairs, &gapXPairs, &gapYPairs, TRUE, TRUE);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                          sX.length, sY.length, score, polishParams->p);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                                                               stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    symbolString_destruct(sX);
    symbolString_destruct(sY);

    return finalAlignedPairs;
}


stList *alignConsensusAndTruthRLEWithKmerAnchors(RleString *consensusStr, RleString *truthStr, double *score,
                                                 PolishParams *polishParams) {


    // Symbol strings
    SymbolString sX = rleString_constructSymbolString(consensusStr, 0, consensusStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    SymbolString sY = rleString_constructSymbolString(truthStr, 0, truthStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    uint16_t apScore = 0;

    // Run the alignment
    stList *alignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapXPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapYPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *anchorPairs = getKmerAlignmentAnchors(sX, sY, (uint64_t) polishParams->p->diagonalExpansion);
    for (int64_t i = 0; i<stList_length(anchorPairs); i++) {
        int64_t *ap = stList_get(anchorPairs, i);
        ap[3] = PAIR_ALIGNMENT_PROB_1;
    }

    getAlignedPairsWithIndelsUsingAnchors(polishParams->stateMachineForForwardStrandRead, sX, sY, anchorPairs,
                                          polishParams->p, &alignedPairs, &gapXPairs, &gapYPairs, TRUE, TRUE);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                          sX.length, sY.length, score, polishParams->p);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                                                               stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    stList_destruct(anchorPairs);
    symbolString_destruct(sX);
    symbolString_destruct(sY);

    return finalAlignedPairs;
}


stList *alignConsensusAndTruthRLEWithSSWAnchors(RleString *consensusStr, RleString *truthStr, double *score,
                                                 PolishParams *polishParams) {


    // Symbol strings
    SymbolString sX = rleString_constructSymbolString(consensusStr, 0, consensusStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    SymbolString sY = rleString_constructSymbolString(truthStr, 0, truthStr->length, polishParams->alphabet,
                                                      TRUE, (uint64_t) polishParams->repeatSubMatrix->maximumRepeatLength);
    uint16_t apScore = 0;

    // Run the alignment
    stList *alignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapXPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapYPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *anchorPairs = alignConsensusAndTruthSSW(consensusStr->rleString, truthStr->rleString, &apScore);
    for (int64_t i = 0; i<stList_length(anchorPairs); i++) {
        int64_t *ap = stList_get(anchorPairs, i);
        ap[3] = PAIR_ALIGNMENT_PROB_1;
    }

    getAlignedPairsWithIndelsUsingAnchors(polishParams->stateMachineForForwardStrandRead, sX, sY, anchorPairs,
                                          polishParams->p, &alignedPairs, &gapXPairs, &gapYPairs, TRUE, TRUE);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
                                                                          sX.length, sY.length, score, polishParams->p);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                                                               stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    stList_destruct(anchorPairs);
    symbolString_destruct(sX);
    symbolString_destruct(sY);

    return finalAlignedPairs;
}

stList *alignConsensusAndTruthCPECAN(char *consensusStr, char *truthStr, double *score, PolishParams *polishParams) {


    // Symbol strings
    SymbolString sX = symbolString_construct(consensusStr, 0, strlen(consensusStr), polishParams->alphabet);
    SymbolString sY = symbolString_construct(truthStr, 0, strlen(truthStr), polishParams->alphabet);
    uint16_t apScore = 0;

    // Use default state machine for alignment
    StateMachine *sM = stateMachine3_constructNucleotide(threeState);

    // Run the alignment
    stList *alignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapXPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    stList *gapYPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);

    getAlignedPairsWithIndels(sM, sX, sY, polishParams->p, &alignedPairs, &gapXPairs,
            &gapYPairs, TRUE, TRUE);
    stList *meaAlignedPairs = getMaximalExpectedAccuracyPairwiseAlignment(alignedPairs, gapXPairs, gapYPairs,
            sX.length, sY.length, score, polishParams->p);

    // refactor
    stList *finalAlignedPairs = stList_construct3(0, (void(*)(void*))stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(meaAlignedPairs); i++) {
        stIntTuple *ap = stList_get(meaAlignedPairs, i);
        stList_append(finalAlignedPairs, stIntTuple_construct3(stIntTuple_get(ap, 1),
                stIntTuple_get(ap, 2), stIntTuple_get(ap, 0)));
    }

    // Cleanup
    stList_destruct(meaAlignedPairs);
    stList_destruct(alignedPairs);
    stList_destruct(gapXPairs);
    stList_destruct(gapYPairs);
    symbolString_destruct(sX);
    symbolString_destruct(sY);
    stateMachine_destruct(sM);

    return finalAlignedPairs;
}




// this function taken from https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/master/src/example.c
stList *alignConsensusAndTruthSSW(char *consensusStr, char *truthStr, uint16_t *score) {

    int64_t l, m, k;
    uint8_t match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;    // default parameters for genome sequence alignment

    int64_t consensusLen = strlen(consensusStr);
    int64_t truthLen = strlen(truthStr);

    s_profile *profile;

    int8_t *num = st_calloc(consensusLen, sizeof(int8_t));    // the read sequence represented in numbers
    int8_t *ref_num = st_calloc(truthLen, sizeof(int8_t));  // the read sequence represented in numbers
    s_align *result;

    /* This table is used to transform nucleotide letters into numbers. */
    static const int8_t nt_table[128] = {
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };

    // initialize scoring matrix for genome sequences
    //  A  C  G  T	N (or other ambiguous code)
    //  2 -2 -2 -2 	0	A
    // -2  2 -2 -2 	0	C
    // -2 -2  2 -2 	0	G
    // -2 -2 -2  2 	0	T
    //	0  0  0  0  0	N (or other ambiguous code)
    int8_t *mat = (int8_t *) calloc(25, sizeof(int8_t));
    for (l = k = 0; l < 4; ++l) {
        for (m = 0; m < 4; ++m)
            mat[k++] = (int8_t) (l == m ? match : -mismatch);    /* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base: no penalty
    }
    for (m = 0; m < 5; ++m) mat[k++] = 0;

    for (m = 0; m < consensusLen; ++m) num[m] = nt_table[(int) consensusStr[m]];
    profile = ssw_init(num, (int32_t) consensusLen, mat, 5, 2);
    for (m = 0; m < truthLen; ++m) ref_num[m] = nt_table[(int) truthStr[m]];

    // Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
    result = ssw_align (profile, ref_num, (int32_t) truthLen, gap_open, gap_extension, 1, 0, 0, 15);
    *score = result->score1;

    // Convert from cigar to aligned pairs
    int32_t consensusPos = result->read_begin1;
    int32_t truthPos = result->ref_begin1;
    stList *alignedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    if (result->cigar) {
        for (int32_t cigIdx = 0; cigIdx < result->cigarLen; cigIdx++) {
            char letter = cigar_int_to_op(result->cigar[cigIdx]);
            uint32_t length = cigar_int_to_len(result->cigar[cigIdx]);

            if (letter == 'M') {
                for (int32_t matchIdx = 0; matchIdx < length; matchIdx++) {
                    stList_append(alignedPairs, stIntTuple_construct3(consensusPos, truthPos, 0));
                    consensusPos++;
                    truthPos++;
                }
            } else if (letter == 'I') {
                consensusPos += length;
            } else if (letter == 'D') {
                truthPos += length;
            } else {
                st_errAbort("Got unexpected alignment character '%c' aligning truth and consensus", letter);
            }
        }
    }

    // cleanup
    align_destroy(result);
    init_destroy(profile);
    free(mat);
    free(ref_num);
    free(num);

    return alignedPairs;
}

#define HDF5_FEATURE_SIZE 1000

double **getTwoDArrayDouble(int64_t rowCount, int64_t columnCount, bool zeroValues) {
    double **array = st_calloc(rowCount, sizeof(double *));
    array[0] = (double *) st_calloc(columnCount * rowCount, sizeof(double));
    if (zeroValues) {
        for (int64_t i = 0; i < columnCount * rowCount; i++) {
            (*array)[i] = 0.0;
        }
    }
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

float **getTwoDArrayFloat(int64_t rowCount, int64_t columnCount) {
    float **array = st_calloc(rowCount, sizeof(float *));
    array[0] = (float *) st_calloc(columnCount * rowCount, sizeof(float));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

uint32_t **getTwoDArrayUInt32(int64_t rowCount, int64_t columnCount) {
    uint32_t **array = st_calloc(rowCount, sizeof(uint32_t *));
    array[0] = (uint32_t *) st_calloc(columnCount * rowCount, sizeof(uint32_t));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

uint8_t **getTwoDArrayUInt8(int64_t rowCount, int64_t columnCount) {
    uint8_t **array = st_calloc(rowCount, sizeof(uint8_t *));
    array[0] = (uint8_t *) st_calloc(columnCount * rowCount, sizeof(uint8_t));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}

uint8_t ***getThreeDArrayUInt8(int64_t depthCount, int64_t rowCount, int64_t columnCount) {
    uint8_t ***array = st_calloc(depthCount, sizeof(uint8_t **));
    array[0] = (uint8_t **) st_calloc(depthCount * rowCount, sizeof(uint8_t *));
    array[0][0] = (uint8_t *) st_calloc(depthCount * rowCount * columnCount, sizeof(uint8_t));
    for (int64_t i = 0; i < depthCount; i++) {
        array[i] = array[0] + i * rowCount;
        for (int64_t j = 0; j < rowCount; j++) {
            array[i][j] = (array[0][0] + i * rowCount * columnCount + j * columnCount);
        }
    }
    return array;
}

char **getTwoDArrayChar(int64_t rowCount, int64_t columnCount) {
    char **array = st_calloc(rowCount, sizeof(char *));
    array[0] = (char *) st_calloc(columnCount * rowCount, sizeof(char));
    for (int64_t i = 1; i < rowCount; i++) {
        array[i] = array[0] + i * columnCount;
    }
    return array;
}


// todo rethink this
#define MAX_TOTAL_WEIGHT 256.0

uint8_t convertTotalWeightToUInt8(double totalWeight) {
    // convert to "depth space"
    totalWeight /= PAIR_ALIGNMENT_PROB_1;
    // cap at 64
    if (totalWeight > MAX_TOTAL_WEIGHT) totalWeight = MAX_TOTAL_WEIGHT;
    // convert to uint8
    return (uint8_t) (totalWeight / MAX_TOTAL_WEIGHT * (UINT8_MAX - 1));
}

uint8_t normalizeWeightToUInt8(double totalWeight, double weight) {
    return (uint8_t) (weight / totalWeight * (UINT8_MAX - 1));
}


void
writeSimpleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                   BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx,
                                   int64_t featureEndIdxInclusive) {

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            featureCount++;
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount,
                   HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 2);
    int64_t columnCount = POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE;
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **imageData = getTwoDArrayUInt8(featureCount, columnCount);
    char **labelCharacterData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayChar(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSimpleWeight *feature = stList_get(features, i);

        // get total weight
        double totalWeight = 0;
        for (int64_t j = 0; j < columnCount; j++) {
            totalWeight += feature->weights[j];
        }

        while (feature != NULL) {
            //positions
            positionData[featureCount][0] = (uint32_t) feature->refPosition;
            positionData[featureCount][1] = (uint32_t) feature->insertPosition;

            // normalization and image data
            normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);
            for (int64_t j = 0; j < columnCount; j++) {
                imageData[featureCount][j] = normalizeWeightToUInt8(totalWeight, feature->weights[j]);
            }

            // nucleotide weights
            int64_t pos;
            for (int64_t symbol = 0; symbol < SYMBOL_NUMBER_NO_N; symbol++) {
                pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, TRUE);
                imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
                pos = PoaFeature_SimpleWeight_charIndex((Symbol) symbol, FALSE);
                imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
            }
            // gap weights
            pos = PoaFeature_SimpleWeight_gapIndex(TRUE);
            imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);
            pos = PoaFeature_SimpleWeight_gapIndex(FALSE);
            imageData[featureCount][pos] = normalizeWeightToUInt8(totalWeight, feature->weights[pos]);

            // potentially labels
            if (outputLabels) {
                Symbol label = alphabet->convertCharToSymbol(feature->label);
                labelCharacterData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ? 0 :
                                                                 label + 1);
            }

            // increment
            featureCount++;
            feature = feature->nextInsert;

        }
    }

    /*
     * Get hdf5 data set up
     */

    hid_t status;
    hid_t stringType = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);

    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 2};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t imageDimension[2] = {featureSize, (hsize_t) columnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t imageSpace = H5Screate_simple(2, imageDimension, NULL);

    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles =
            (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) /
                                   (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }

    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace,
                                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type,
                           H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace,
                                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type,
                           H5S_ALL, H5S_ALL, H5P_DEFAULT, &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type,
                           H5S_ALL, H5S_ALL, H5P_DEFAULT, positionData[chunkFeatureStartIdx]);

        // write normalization data
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationData[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate(group, "image", hdf5FileInfo->uint8Type, imageSpace,
                                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           imageData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterData[chunkFeatureStartIdx]);
            status |= H5Dclose(labelCharacterDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Dclose(imageDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }


    // cleanup
    free(normalizationData[0]);
    free(normalizationData);
    free(imageData[0]);
    free(imageData);
    free(positionData[0]);
    free(positionData);
    status |= H5Tclose(stringType);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(imageSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(labelCharacterSpace);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}

void
writeSplitRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                     BamChunk *bamChunk,
                                     bool outputLabels, stList *features, int64_t featureStartIdx,
                                     int64_t featureEndIdxInclusive, const int64_t maxRunLength) {

    herr_t status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSplitRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            PoaFeatureSplitRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount,
                   HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 3);
    int64_t rleNucleotideColumnCount = ((SYMBOL_NUMBER - 1) * (maxRunLength + 1) + 1) * 2;
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **imageData = getTwoDArrayUInt8(featureCount, rleNucleotideColumnCount);
    uint8_t **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayUInt8(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureSplitRleWeight *refFeature = stList_get(features, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < rleNucleotideColumnCount; j++) {
            totalWeight += refFeature->weights[j];
        }

        // iterate over all insert features
        PoaFeatureSplitRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureSplitRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionData[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionData[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionData[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t j = 0; j < rleNucleotideColumnCount; j++) {
                    imageData[featureCount][j] =
                            normalizeWeightToUInt8(totalWeight, rlFeature->weights[j]);
                }

                // labels
                if (outputLabels) {
                    Symbol label = alphabet->convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : label + 1);
                    labelRunLengthData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : rlFeature->labelRunLength);
                    if (labelRunLengthData[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthData[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t imageDimension[2] = {featureSize, (hsize_t) rleNucleotideColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t imageSpace = H5Screate_simple(2, imageDimension, NULL);

    hid_t stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles =
            (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) /
                                   (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           positionData[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate(group, "image", hdf5FileInfo->uint8Type, imageSpace, H5P_DEFAULT, H5P_DEFAULT,
                                       H5P_DEFAULT);
        status |= H5Dwrite(imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           imageData[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT,
                                                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate(group, "label_run_length", hdf5FileInfo->uint8Type,
                                                    labelRunLengthSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose(labelCharacterDataset);
            status |= H5Dclose(labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(imageDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }

    // cleanup
    free(imageData[0]);
    free(imageData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(imageSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(labelRunLengthSpace);
    status |= H5Sclose(labelCharacterSpace);
    status |= H5Tclose(stringType);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}


void
writeChannelRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                       BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx,
                                       int64_t featureEndIdxInclusive, const int64_t maxRunLength) {

    herr_t status = 0;

    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureChannelRleWeight *feature = stList_get(features, i);
        while (feature != NULL) {
            PoaFeatureChannelRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }
    if (featureCount < HDF5_FEATURE_SIZE && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" less than minimum of %d\n", logIdentifier, featureCount,
                   HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // sizes
    int64_t nucleotideColumnCount = SYMBOL_NUMBER * 2; //ACGTGap x {fwd,bwd}
    int64_t runLengthColumnCount = (maxRunLength + 1) * 2; //(runLenght + 0) x {fwd,bwd}

    // get all feature data into an array
    uint32_t **positionData = getTwoDArrayUInt32(featureCount, 3);
    uint8_t **normalizationData = getTwoDArrayUInt8(featureCount, 1);
    uint8_t **nucleotideData = getTwoDArrayUInt8(featureCount, nucleotideColumnCount);
    uint8_t ***runLengthData = getThreeDArrayUInt8(featureCount, runLengthColumnCount, (SYMBOL_NUMBER - 1));
    //uint8_t ***runLengthData = getThreeDArrayUInt8(featureCount, (SYMBOL_NUMBER - 1), runLengthColumnCount);
    uint8_t **labelCharacterData = NULL;
    uint8_t **labelRunLengthData = NULL;
    if (outputLabels) {
        labelCharacterData = getTwoDArrayUInt8(featureCount, 1);
        labelRunLengthData = getTwoDArrayUInt8(featureCount, 1);
    }

    // add all data to features
    featureCount = 0;
    for (int64_t i = featureStartIdx; i <= featureEndIdxInclusive; i++) {
        PoaFeatureChannelRleWeight *refFeature = stList_get(features, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < nucleotideColumnCount; j++) {
            totalWeight += refFeature->nucleotideWeights[j];
        }

        // iterate over all insert features
        PoaFeatureChannelRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureChannelRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionData[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionData[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionData[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationData[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t c = 0; c < SYMBOL_NUMBER - 1; c++) {
                    // overall nucl count
                    nucleotideData[featureCount][c * 2 + POS_STRAND_IDX] = normalizeWeightToUInt8(totalWeight,
                                                                                                  rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_charNuclIndex(
                                                                                                          c, TRUE)]);
                    nucleotideData[featureCount][c * 2 + NEG_STRAND_IDX] = normalizeWeightToUInt8(totalWeight,
                                                                                                  rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_charNuclIndex(
                                                                                                          c, FALSE)]);

                    // run length counts
                    for (int64_t r = 0; r <= maxRunLength; r++) {
                        runLengthData[featureCount][r * 2 + POS_STRAND_IDX][c] =
                                normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                                        PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, TRUE)]);
                        runLengthData[featureCount][r * 2 + NEG_STRAND_IDX][c] =
                                normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                                        PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, FALSE)]);

                        //runLengthData[featureCount][c][r * 2 + POS_STRAND_IDX] =
                        //        normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                        //                PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, TRUE)]);
                        //runLengthData[featureCount][c][r * 2 + NEG_STRAND_IDX] =
                        //        normalizeWeightToUInt8(totalWeight, rlFeature->runLengthWeights[
                        //                PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, c, r, FALSE)]);
                    }
                }
                // gap counts
                nucleotideData[featureCount][SYMBOL_NUMBER_NO_N * 2 + 0 + POS_STRAND_IDX] = normalizeWeightToUInt8(
                        totalWeight, rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE)]);
                nucleotideData[featureCount][SYMBOL_NUMBER_NO_N * 2 + NEG_STRAND_IDX] = normalizeWeightToUInt8(
                        totalWeight, rlFeature->nucleotideWeights[PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE)]);

                // labels
                if (outputLabels) {
                    Symbol label = alphabet->convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : label + 1);
                    labelRunLengthData[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : rlFeature->labelRunLength);
                    if (labelRunLengthData[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthData[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */


    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t nucleotideDimension[2] = {featureSize, (hsize_t) nucleotideColumnCount};
    hsize_t runLengthDimension[3] = {featureSize, (hsize_t) runLengthColumnCount, (hsize_t) SYMBOL_NUMBER - 1};
//    hsize_t runLengthDimension[3] = {featureSize, (hsize_t) SYMBOL_NUMBER - 1, (hsize_t) runLengthColumnCount};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t nucleotideSpace = H5Screate_simple(2, nucleotideDimension, NULL);
    hid_t runLengthSpace = H5Screate_simple(3, runLengthDimension, NULL);

    hid_t stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);


    /*
     * Write features to files
     */

    // each file must have exactly 1000 features
    int64_t totalFeatureFiles =
            (int64_t) (featureCount / HDF5_FEATURE_SIZE) + (featureCount % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffset = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffset = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFiles - featureCount) /
                                   (int64_t) (featureCount / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFiles; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffset * featureIndex);
        if (featureIndex + 1 == totalFeatureFiles && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           positionData[chunkFeatureStartIdx]);

        // write nucl, rl, and norm data
        hid_t nucleotideDataset = H5Dcreate(group, "nucleotide", hdf5FileInfo->uint8Type, nucleotideSpace, H5P_DEFAULT,
                                            H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(nucleotideDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           nucleotideData[chunkFeatureStartIdx]);
        hid_t runLengthDataset = H5Dcreate(group, "runLengths", hdf5FileInfo->uint8Type, runLengthSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(runLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           runLengthData[chunkFeatureStartIdx][0]);
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationData[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT,
                                                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterData[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate(group, "label_run_length", hdf5FileInfo->uint8Type,
                                                    labelRunLengthSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelRunLengthData[chunkFeatureStartIdx]);

            status |= H5Dclose(labelCharacterDataset);
            status |= H5Dclose(labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(nucleotideDataset);
        status |= H5Dclose(runLengthDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }

    // cleanup
    free(nucleotideData[0]);
    free(nucleotideData);
    free(runLengthData[0][0]);
    free(runLengthData[0]);
    free(runLengthData);
    free(normalizationData[0]);
    free(normalizationData);
    free(positionData[0]);
    free(positionData);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(nucleotideSpace);
    status |= H5Sclose(runLengthSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(labelRunLengthSpace);
    status |= H5Sclose(labelCharacterSpace);
    status |= H5Tclose(stringType);
    if (outputLabels) {
        free(labelCharacterData[0]);
        free(labelCharacterData);
        free(labelRunLengthData[0]);
        free(labelRunLengthData);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}

void writeDiploidRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo* hdf5FileInfo,
        char *outputFileBase, BamChunk *bamChunk, bool outputLabels,
        stList *featuresH1, int64_t featureStartIdxH1, int64_t featureEndIdxInclusiveH1,
        stList *featuresH2, int64_t featureStartIdxH2, int64_t featureEndIdxInclusiveH2,
        const int64_t maxRunLength) {

    herr_t status = 0;
    /*
     * Get feature data set up
     */

    // count features, create feature array
    uint64_t featureCountH1 = 0;
    uint64_t featureCountH2 = 0;
    for (int64_t i = featureStartIdxH1; i <= featureEndIdxInclusiveH1; i++) {
        PoaFeatureDiploidRleWeight *feature = stList_get(featuresH1, i);
        while (feature != NULL) {
            PoaFeatureDiploidRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCountH1++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }
    for (int64_t i = featureStartIdxH2; i <= featureEndIdxInclusiveH2; i++) {
        PoaFeatureDiploidRleWeight *feature = stList_get(featuresH2, i);
        while (feature != NULL) {
            PoaFeatureDiploidRleWeight *rlFeature = feature;
            while (rlFeature != NULL) {
                featureCountH2++;
                rlFeature = rlFeature->nextRunLength;
            }
            feature = feature->nextInsert;
        }
    }

    // don't write small feature sets when training
    if ((featureCountH1 < HDF5_FEATURE_SIZE || featureCountH2 < HDF5_FEATURE_SIZE) && outputLabels) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Feature count %"PRId64" / %"PRId64" less than minimum of %d\n", logIdentifier,
                featureCountH1, featureCountH2, HDF5_FEATURE_SIZE);
        free(logIdentifier);
        return;
    }

    // get all feature data into an array
    uint32_t **positionDataH1 = getTwoDArrayUInt32(featureCountH1, 3);
    uint32_t **positionDataH2 = getTwoDArrayUInt32(featureCountH2, 3);
    int64_t rleNucleotideStrandColumnCount_singleHap = ((SYMBOL_NUMBER - 1) * (maxRunLength + 1) + 1) * 2; //(nucl*(rl+0)+gap)*strand
    uint8_t **normalizationDataH1 = getTwoDArrayUInt8(featureCountH1, 1);
    uint8_t **normalizationDataH2 = getTwoDArrayUInt8(featureCountH2, 1);
    uint8_t **imageDataH1 = getTwoDArrayUInt8(featureCountH1, rleNucleotideStrandColumnCount_singleHap * 2);
    uint8_t **imageDataH2 = getTwoDArrayUInt8(featureCountH2, rleNucleotideStrandColumnCount_singleHap * 2);
    uint8_t **labelCharacterDataH1 = NULL;
    uint8_t **labelCharacterDataH2 = NULL;
    uint8_t **labelRunLengthDataH1 = NULL;
    uint8_t **labelRunLengthDataH2 = NULL;
    if (outputLabels) {
        labelCharacterDataH1 = getTwoDArrayUInt8(featureCountH1, 1);
        labelCharacterDataH2 = getTwoDArrayUInt8(featureCountH2, 1);
        labelRunLengthDataH1 = getTwoDArrayUInt8(featureCountH1, 1);
        labelRunLengthDataH2 = getTwoDArrayUInt8(featureCountH2, 1);
    }


    /*
     *   Set up feature data in HDF5 format
     */

    // add feature data H1
    int64_t featureCount = 0;
    for (int64_t i = featureStartIdxH1; i <= featureEndIdxInclusiveH1; i++) {
        PoaFeatureDiploidRleWeight *refFeature = stList_get(featuresH1, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < rleNucleotideStrandColumnCount_singleHap; j++) {
            totalWeight += refFeature->weightsHOn[j] + refFeature->weightsHOff[j];
        }

        // iterate over all insert features
        PoaFeatureDiploidRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureDiploidRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionDataH1[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionDataH1[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionDataH1[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationDataH1[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t j = 0; j < rleNucleotideStrandColumnCount_singleHap; j++) {
                    imageDataH1[featureCount][j] = normalizeWeightToUInt8(totalWeight, rlFeature->weightsHOn[j]);
                    imageDataH1[featureCount][j+rleNucleotideStrandColumnCount_singleHap] =
                            normalizeWeightToUInt8(totalWeight, rlFeature->weightsHOff[j]);
                }

                // labels
                if (outputLabels) {
                    Symbol label = alphabet->convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterDataH1[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : label + 1);
                    labelRunLengthDataH1[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : rlFeature->labelRunLength);
                    if (labelRunLengthDataH1[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in H1 chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthDataH1[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    // add feature data H2
    featureCount = 0;
    for (int64_t i = featureStartIdxH2; i <= featureEndIdxInclusiveH2; i++) {
        PoaFeatureDiploidRleWeight *refFeature = stList_get(featuresH2, i);

        // total weight is calculated for the very first refPos feature, used for all inserts and run lengths
        double totalWeight = 0;
        for (int64_t j = 0; j < rleNucleotideStrandColumnCount_singleHap; j++) {
            totalWeight += refFeature->weightsHOn[j] + refFeature->weightsHOff[j];
        }

        // iterate over all insert features
        PoaFeatureDiploidRleWeight *insFeature = refFeature;
        while (insFeature != NULL) {

            // iterate over all run length features
            PoaFeatureDiploidRleWeight *rlFeature = insFeature;
            while (rlFeature != NULL) {
                // position
                positionDataH2[featureCount][0] = (uint32_t) rlFeature->refPosition;
                positionDataH2[featureCount][1] = (uint32_t) rlFeature->insertPosition;
                positionDataH2[featureCount][2] = (uint32_t) rlFeature->runLengthPosition;

                // normalization
                normalizationDataH2[featureCount][0] = convertTotalWeightToUInt8(totalWeight);

                // copy weights over (into normalized uint8 space)
                for (int64_t j = 0; j < rleNucleotideStrandColumnCount_singleHap; j++) {
                    imageDataH2[featureCount][j] = normalizeWeightToUInt8(totalWeight, rlFeature->weightsHOn[j]);
                    imageDataH2[featureCount][j+rleNucleotideStrandColumnCount_singleHap] =
                            normalizeWeightToUInt8(totalWeight, rlFeature->weightsHOff[j]);
                }

                // labels
                if (outputLabels) {
                    Symbol label = alphabet->convertCharToSymbol(rlFeature->labelChar);
                    labelCharacterDataH2[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : label + 1);
                    labelRunLengthDataH2[featureCount][0] = (uint8_t) (alphabet->convertSymbolToChar(label) == 'N' ?
                                                                     0 : rlFeature->labelRunLength);
                    if (labelRunLengthDataH2[featureCount][0] > maxRunLength) {
                        st_errAbort("Encountered run length of %d (max %"PRId64") in H2 chunk %s:%"PRId64"-%"PRId64,
                                    labelRunLengthDataH2[featureCount][0], maxRunLength, bamChunk->refSeqName,
                                    bamChunk->chunkBoundaryStart, bamChunk->chunkBoundaryEnd);
                    }
                }

                // iterate
                featureCount++;
                rlFeature = rlFeature->nextRunLength;
            }
            insFeature = insFeature->nextInsert;
        }
    }

    /*
     * Get hdf5 data set up
     */
    // so that we can produce chunks smaller than HDF5_FEATURE_SIZE (not used during training)
    hsize_t featureSize = (hsize_t) (featureCount < HDF5_FEATURE_SIZE ? featureCount : HDF5_FEATURE_SIZE);

    hsize_t metadataDimension[1] = {1};
    hsize_t postionDimension[2] = {featureSize, 3};
    hsize_t labelCharacterDimension[2] = {featureSize, 1};
    hsize_t labelRunLengthDimension[2] = {featureSize, 1};
    hsize_t normalizationDimension[2] = {featureSize, 1};
    hsize_t imageDimension[2] = {featureSize, (hsize_t) rleNucleotideStrandColumnCount_singleHap * 2};

    hid_t metadataSpace = H5Screate_simple(1, metadataDimension, NULL);
    hid_t positionSpace = H5Screate_simple(2, postionDimension, NULL);
    hid_t labelCharacterSpace = H5Screate_simple(2, labelCharacterDimension, NULL);
    hid_t labelRunLengthSpace = H5Screate_simple(2, labelRunLengthDimension, NULL);
    hid_t normalizationSpace = H5Screate_simple(2, normalizationDimension, NULL);
    hid_t imageSpace = H5Screate_simple(2, imageDimension, NULL);

    hid_t stringType = H5Tcopy(H5T_C_S1);
    H5Tset_size(stringType, strlen(bamChunk->refSeqName) + 1);


    /*
     *   Write features to files
     */

    // Write for hap1
    int64_t totalFeatureFilesH1 =
            (int64_t) (featureCountH1 / HDF5_FEATURE_SIZE) + (featureCountH1 % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffsetH1 = 0;
    if (featureCount >= HDF5_FEATURE_SIZE) {
        featureOffsetH1 = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFilesH1 - featureCountH1) /
                                   (int64_t) (featureCountH1 / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFilesH1; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffsetH1 * featureIndex);
        if (featureIndex + 1 == totalFeatureFilesH1 && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.H1.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);
        hid_t haplotypeIndexDataset = H5Dcreate(group, "haplotype", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        int64_t haplotype = 1;
        status |= H5Dwrite(haplotypeIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &haplotype);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           positionDataH1[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate(group, "image", hdf5FileInfo->uint8Type, imageSpace, H5P_DEFAULT, H5P_DEFAULT,
                                       H5P_DEFAULT);
        status |= H5Dwrite(imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           imageDataH1[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationDataH1[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT,
                                                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterDataH1[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate(group, "label_run_length", hdf5FileInfo->uint8Type,
                                                    labelRunLengthSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelRunLengthDataH1[chunkFeatureStartIdx]);

            status |= H5Dclose(labelCharacterDataset);
            status |= H5Dclose(labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(haplotypeIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(imageDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }

    // Write for hap2
    int64_t totalFeatureFilesH2 =
            (int64_t) (featureCountH2 / HDF5_FEATURE_SIZE) + (featureCountH2 % HDF5_FEATURE_SIZE == 0 ? 0 : 1);
    int64_t featureOffsetH2 = 0;
    if (featureCountH2 >= HDF5_FEATURE_SIZE) {
        featureOffsetH2 = (int64_t) ((HDF5_FEATURE_SIZE * totalFeatureFilesH2 - featureCountH2) /
                                   (int64_t) (featureCountH2 / HDF5_FEATURE_SIZE));
    }
    for (int64_t featureIndex = 0; featureIndex < totalFeatureFilesH2; featureIndex++) {
        // get start pos
        int64_t chunkFeatureStartIdx = (HDF5_FEATURE_SIZE * featureIndex) - (featureOffsetH2 * featureIndex);
        if (featureIndex + 1 == totalFeatureFilesH2 && featureCount >= HDF5_FEATURE_SIZE) {
            chunkFeatureStartIdx = featureCount - HDF5_FEATURE_SIZE;
        }

        // create group
        char *outputGroup = stString_print("images/%s.H2.%"PRId64, outputFileBase, featureIndex);
        hid_t group = H5Gcreate(hdf5FileInfo->file, outputGroup, hdf5FileInfo->groupPropertyList, H5P_DEFAULT,
                                H5P_DEFAULT);

        // write metadata
        hid_t contigDataset = H5Dcreate(group, "contig", stringType, metadataSpace, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT);
        status |= H5Dwrite(contigDataset, stringType, H5S_ALL, H5S_ALL, H5P_DEFAULT, bamChunk->refSeqName);
        hid_t contigStartDataset = H5Dcreate(group, "contig_start", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                             H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigStartDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryStart);
        hid_t contigEndDataset = H5Dcreate(group, "contig_end", hdf5FileInfo->int64Type, metadataSpace, H5P_DEFAULT,
                                           H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(contigEndDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &bamChunk->chunkBoundaryEnd);
        hid_t chunkIndexDataset = H5Dcreate(group, "feature_chunk_idx", hdf5FileInfo->int64Type, metadataSpace,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(chunkIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &featureIndex);
        hid_t haplotypeIndexDataset = H5Dcreate(group, "haplotype", hdf5FileInfo->int64Type, metadataSpace,
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        int64_t haplotype = 2;
        status |= H5Dwrite(haplotypeIndexDataset, hdf5FileInfo->int64Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &haplotype);

        // write position info
        hid_t positionDataset = H5Dcreate(group, "position", hdf5FileInfo->uint32Type, positionSpace, H5P_DEFAULT,
                                          H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(positionDataset, hdf5FileInfo->uint32Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           positionDataH2[chunkFeatureStartIdx]);

        // write rle data
        hid_t imageDataset = H5Dcreate(group, "image", hdf5FileInfo->uint8Type, imageSpace, H5P_DEFAULT, H5P_DEFAULT,
                                       H5P_DEFAULT);
        status |= H5Dwrite(imageDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           imageDataH2[chunkFeatureStartIdx]);
        hid_t normalizationDataset = H5Dcreate(group, "normalization", hdf5FileInfo->uint8Type, normalizationSpace,
                                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status |= H5Dwrite(normalizationDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           normalizationDataH2[chunkFeatureStartIdx]);

        // if labels, add all these too
        if (outputLabels) {
            hid_t labelCharacterDataset = H5Dcreate(group, "label_base", hdf5FileInfo->uint8Type, labelCharacterSpace,
                                                    H5P_DEFAULT,
                                                    H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelCharacterDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelCharacterDataH2[chunkFeatureStartIdx]);
            hid_t labelRunLengthDataset = H5Dcreate(group, "label_run_length", hdf5FileInfo->uint8Type,
                                                    labelRunLengthSpace,
                                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status |= H5Dwrite(labelRunLengthDataset, hdf5FileInfo->uint8Type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                               labelRunLengthDataH2[chunkFeatureStartIdx]);

            status |= H5Dclose(labelCharacterDataset);
            status |= H5Dclose(labelRunLengthDataset);
        }

        // cleanup
        status |= H5Dclose(contigDataset);
        status |= H5Dclose(contigStartDataset);
        status |= H5Dclose(contigEndDataset);
        status |= H5Dclose(chunkIndexDataset);
        status |= H5Dclose(haplotypeIndexDataset);
        status |= H5Dclose(positionDataset);
        status |= H5Dclose(imageDataset);
        status |= H5Dclose(normalizationDataset);
        status |= H5Gclose(group);
        free(outputGroup);
    }

    // cleanup
    free(imageDataH1[0]);
    free(imageDataH2[0]);
    free(imageDataH1);
    free(imageDataH2);
    free(normalizationDataH1[0]);
    free(normalizationDataH2[0]);
    free(normalizationDataH1);
    free(normalizationDataH2);
    free(positionDataH1[0]);
    free(positionDataH2[0]);
    free(positionDataH1);
    free(positionDataH2);
    status |= H5Sclose(metadataSpace);
    status |= H5Sclose(positionSpace);
    status |= H5Sclose(imageSpace);
    status |= H5Sclose(normalizationSpace);
    status |= H5Sclose(labelRunLengthSpace);
    status |= H5Sclose(labelCharacterSpace);
    status |= H5Tclose(stringType);
    if (outputLabels) {
        free(labelCharacterDataH1[0]);
        free(labelCharacterDataH2[0]);
        free(labelCharacterDataH1);
        free(labelCharacterDataH2);
        free(labelRunLengthDataH1[0]);
        free(labelRunLengthDataH2[0]);
        free(labelRunLengthDataH1);
        free(labelRunLengthDataH2);
    }

    if (status) {
        char *logIdentifier = getLogIdentifier();
        st_logInfo(" %s Error writing HELEN features to HDF5 files: %s\n", logIdentifier, outputFileBase);
        free(logIdentifier);
    }
}


HelenFeatureHDF5FileInfo *HelenFeatureHDF5FileInfo_construct(char *filename) {
    HelenFeatureHDF5FileInfo *fileInfo = st_calloc(1, sizeof(HelenFeatureHDF5FileInfo));
    fileInfo->filename = stString_copy(filename);
    fileInfo->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    fileInfo->int64Type = H5Tcopy(H5T_NATIVE_UINT32);
    H5Tset_order(fileInfo->int64Type, H5T_ORDER_LE);
    fileInfo->uint32Type = H5Tcopy(H5T_NATIVE_UINT32);
    H5Tset_order(fileInfo->uint32Type, H5T_ORDER_LE);
    fileInfo->uint8Type = H5Tcopy(H5T_NATIVE_UINT8);
    H5Tset_order(fileInfo->uint8Type, H5T_ORDER_LE);
    fileInfo->floatType = H5Tcopy(H5T_NATIVE_FLOAT);
    H5Tset_order(fileInfo->floatType, H5T_ORDER_LE);
    fileInfo->groupPropertyList = H5Pcreate(H5P_LINK_CREATE);
    H5Pset_create_intermediate_group(fileInfo->groupPropertyList, 1);
    return fileInfo;
}

void HelenFeatureHDF5FileInfo_destruct(HelenFeatureHDF5FileInfo *fileInfo) {
    free(fileInfo->filename);
    H5Tclose(fileInfo->int64Type);
    H5Tclose(fileInfo->uint32Type);
    H5Tclose(fileInfo->uint8Type);
    H5Tclose(fileInfo->floatType);
    H5Pclose(fileInfo->groupPropertyList);
    H5Fclose(fileInfo->file);
    free(fileInfo);
}

HelenFeatureHDF5FileInfo **openHelenFeatureHDF5FilesByThreadCount(char *filenameBase, int64_t threadCount) {
    HelenFeatureHDF5FileInfo **infoArray = st_calloc(threadCount, sizeof(HelenFeatureHDF5FileInfo *));
    for (int64_t i = 0; i < threadCount; i++) {
        char *filename = stString_print("%s.T%02"PRId64".h5", filenameBase, i);
        infoArray[i] = HelenFeatureHDF5FileInfo_construct(filename);
        free(filename);
    }
    return infoArray;
}

#endif