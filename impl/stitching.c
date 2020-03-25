//
// Created by Benedict Paten on 3/14/20.
//
// Code for stitching together "chunks" of inferred sequence
//

#include "margin.h"

/*
 * ChunkToStitch
 */

typedef struct _chunkToStitch {
    /*
     * Object for managing the output of a polished sequence.
     */
    bool startOfSequence; // Indicates if it is the first chunk in a sequence
    int64_t chunkOrdinal; // The index of the chunk in the sequence of chunks
    char *seqName1; // The name of the primary sequence being polished
    char *seqHap1; // The primary sequence

    char *seqName2; // The name of secondary haplotype sequence being polished, may be null.
    char *seqHap2; // The secondary haplotype, may be null.

    // Following from the CSV files, each line corresponds to a position in the seqHap1

    // Both these will only be present if phasing
    stList *readsHap1Lines; // Reads from primary sequence
    stList *readsHap2Lines; // Reads from the secondary sequence

    // Lines from the POA
    stList *poaHap1StringsLines;
    stList *poaHap2StringsLines;

    // Lines from the repeat count file
    stList *repeatCountLines;
} ChunkToStitch;

int chunkToStitch_cmp(ChunkToStitch *chunk1, ChunkToStitch *chunk2) {
    /*
     * Compares to chunks by their ordinal in the output order.
     */
    return chunk1->chunkOrdinal < chunk2->chunkOrdinal ? -1 : chunk1->chunkOrdinal > chunk2->chunkOrdinal ? 1 : 0;
}

void chunkToStitch_destruct(ChunkToStitch *chunkToStitch) {
    free(chunkToStitch->seqName1);
    free(chunkToStitch->seqHap1);

    // Second sequence and remaining fields are optional
    if(chunkToStitch->seqName2 != NULL) {
        free(chunkToStitch->seqName2);
        free(chunkToStitch->seqHap2);
    }

    if(chunkToStitch->readsHap1Lines != NULL) {
        stList_destruct(chunkToStitch->readsHap1Lines);
    }
    if(chunkToStitch->readsHap2Lines != NULL) {
        stList_destruct(chunkToStitch->readsHap2Lines);
    }

    if(chunkToStitch->poaHap1StringsLines != NULL) {
        stList_destruct(chunkToStitch->poaHap1StringsLines);
    }
    if(chunkToStitch->poaHap2StringsLines != NULL) {
        stList_destruct(chunkToStitch->poaHap2StringsLines);
    }

    if(chunkToStitch->repeatCountLines != NULL) {
        stList_destruct(chunkToStitch->repeatCountLines);
    }

    free(chunkToStitch);
}

stList *readChunk(FILE *fh, char **seqName, int64_t *chunkOrdinal) {
    /*
     * Reads "chunks" from a file. Each chunk starts with a line formatted as:
     * SEQ_NAME\tCHUNK_ORDINAL\tLINE_NUMBER\n
     * where SEQ_NAME is the name of the string, containing any characters other than white space,
     * CHUNK_ORDINAL is the order of the chunk in the output,
     * and LINE_NUMBER is an integer >= 0 that gives the remaining number of lines in the chunk to read.
     * Returns the chunk as a list of lines, in order, and initializes the seqName argument to be a string representing
     * seqName.
     *
     * If the EOF is reached will return NULL, set seqName to NULL and chunkOrdinal to -1.
     *
     * Newlines characters are omitted from the ends of each line string.
     */
    char *headerLine = stFile_getLineFromFile(fh);
    if(headerLine == NULL) {
        *chunkOrdinal = -1;
        *seqName = NULL; // Seq seqName to NULL;
        return NULL;
    }

    stList *tokens = stString_splitByString(headerLine, "\t");

    if(stList_length(tokens) != 3) {
        st_errAbort("Expected three tokens in header line, got %" PRIi64 "", stList_length(tokens));
    }

    *seqName = stList_removeFirst(tokens); // Set seqName

    *chunkOrdinal = strtol(stList_get(tokens, 0), NULL, 10); // Get chunk ordinal

    int64_t lineNo = strtol(stList_peek(tokens), NULL, 10); // Get line number

    stList *lines = stList_construct3(0, free);
    for(int64_t i=0; i<lineNo; i++) {
        char *bodyLine = stFile_getLineFromFile(fh);
        if(bodyLine == NULL) {
            st_errAbort("Failed to read body line from chunk, line %" PRIi64 " of %" PRIi64 " lines\n", i, lineNo);
        }
        stList_append(lines, bodyLine);
    }

    // Cleanup
    stList_destruct(tokens);
    free(headerLine);

    return lines;
}

stList *readChunk2(FILE *fh, char *expectedSequenceName, int64_t expectedChunkOrdinal) {
    /*
     * As readChunk, but creates an error if no chunk is found or the read seqName is not equal to expectedSequenceName.
     */
    char *seqName;
    int64_t chunkOrdinal;
    stList *lines = readChunk(fh, &seqName, &chunkOrdinal);
    if(lines == NULL) {
        st_errAbort("Got no chunk when one was expected\n");
    }
    if(!stString_eq(seqName, expectedSequenceName)) {
        st_errAbort("Got an unexpected sequence name: %s in reading chunk (expected: %s)\n", seqName, expectedSequenceName);
    }
    if(expectedChunkOrdinal != chunkOrdinal) {
        st_errAbort("Got an unexpected chunk ordinal (%" PRIi64 ") in reading chunk (expected: %" PRIi64 ")\n", chunkOrdinal, expectedChunkOrdinal);
    }
    free(seqName);
    return lines;
}

bool chunkToStitch_readSequenceChunk(FILE *fh, ChunkToStitch *chunk, bool phased) {
    /*
     * Reads a "chunk" from a sequence containing file, adding it to chunk. Returns non-zero if no more chunks remain.
     */

    // Read the next set of lines from the file
    stList *lines = readChunk(fh, &chunk->seqName1, &chunk->chunkOrdinal);

    // If we get nothing we have exhausted the file
    if(lines == NULL) {
        return 0;
    }

    chunk->seqHap1 = stString_join2("", lines);   // Concatenate the lines to make the sequence
    stList_destruct(lines); // Cleanup

    if(phased) {
        int64_t i;
        lines = readChunk(fh, &chunk->seqName2, &i);
        if(lines == NULL) {
            st_errAbort("Error trying get alt sequence from chunk");
        }
        if(i != chunk->chunkOrdinal) {
            st_errAbort("Got an unexpected chunk ordinal (%" PRIi64 ") in reading second haplotype chunk (expected: %" PRIi64 ")\n", i, chunk->chunkOrdinal);
        }
        chunk->seqHap2 = stString_join2("", lines);
        stList_destruct(lines); // Cleanup
    }

    return 0;
}

void chunkToStitch_readPoaChunk(FILE *fh, ChunkToStitch *chunk, bool phased) {
    /*
     * Reads a chunk from a POA containing file.
     */
    chunk->poaHap1StringsLines = readChunk2(fh, chunk->seqName1, chunk->chunkOrdinal);
    if(phased) {
        chunk->poaHap2StringsLines = readChunk2(fh, chunk->seqName2, chunk->chunkOrdinal);
    }
}

void chunkToStitch_readReadPhasingChunk(FILE *fh, ChunkToStitch *chunk) {
    /*
     * Reads a read phasing chunk.
     */
    chunk->readsHap1Lines = readChunk2(fh, chunk->seqName1, chunk->chunkOrdinal);
    chunk->readsHap2Lines = readChunk2(fh, chunk->seqName2, chunk->chunkOrdinal);
}

void chunkToStitch_readRepeatCountChunk(FILE *fh, ChunkToStitch *chunk) {
    /*
     * Reads repeat counts chunk
     */
    chunk->repeatCountLines = readChunk2(fh, chunk->seqName1, chunk->chunkOrdinal);
}

stSet *getReadNames(stList *readPartitionLines) {
    /*
     * Parse the names of the reads from the lines of output representing the relative read phasing and return as a set
     * of strings.
     */
    stSet *readNames = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    for(int64_t i=1; i<stList_length(readPartitionLines); i++) {
        char *line = stList_get(readPartitionLines, i);
        stList *tokens = stString_splitByString(line, ",");
        stSet_insert(readNames, stList_removeFirst(tokens));
        stList_destruct(tokens);
    }
    return readNames;
}

static void swap(void **a, void **b) {
    void *c = *a;
    *a = *b;
    *b = c;
}

void chunkToStitch_phaseAdjacentChunks(ChunkToStitch *pChunk, ChunkToStitch *chunk) {
    /*
     * Phases chunk so that hap1 in chunk corresponds to hap1 in pChunk (as best as we can tell).
     */

    // Get the names of the reads in the different read sets
    stSet *pChunkHap1Reads = getReadNames(pChunk->readsHap1Lines);
    stSet *pChunkHap2Reads = getReadNames(pChunk->readsHap2Lines);
    stSet *chunkHap1Reads = getReadNames(chunk->readsHap1Lines);
    stSet *chunkHap2Reads = getReadNames(chunk->readsHap2Lines);

    // Calculate the intersection between reads shared between the chunks
    int64_t i = stSet_sizeOfIntersection(pChunkHap1Reads, chunkHap1Reads);
    int64_t j = stSet_sizeOfIntersection(pChunkHap1Reads, chunkHap2Reads);
    int64_t k = stSet_sizeOfIntersection(pChunkHap2Reads, chunkHap1Reads);
    int64_t l = stSet_sizeOfIntersection(pChunkHap2Reads, chunkHap2Reads);

    // Calculate support for the cis (keeping the current relative phasing) and the trans (switching the phasing) configurations
    int64_t cisPhase = i + l;
    int64_t transPhase = j + k;

    // Log the support for the phasing
    st_logDebug("In phasing between two chunks got %" PRIi64 " hap1 and %" PRIi64 " hap2 reads in pChunk, and got %"
    PRIi64 " hap1 and %" PRIi64 " hap2 reads in chunk and %" PRIi64 " reads in intersections\n",
                stSet_size(pChunkHap1Reads), stSet_size(pChunkHap2Reads), stSet_size(chunkHap1Reads), stSet_size(chunkHap2Reads), i+j+k+l);
    st_logDebug(" Got %" PRIi64 " (%f%) in h11-h21 intersection, %" PRIi64 " (%f%) in h11-h22 intersection\n", i, j, (float)i/(i+j), (float)j/(i+j));
    st_logDebug(" Got %" PRIi64 " (%f%) in h12-h21 intersection, %" PRIi64 " (%f%) in h12-h22 intersection\n", k, l, (float)k/(k+l), (float)l/(k+l));

    // Switch the relative phasing is the trans phase configuration has more support
    if(cisPhase < transPhase) {
        st_logDebug("Flipping phase of chunk\n");
        swap((void *)&chunk->seqName1, (void *)&chunk->seqName2);
        swap((void *)&chunk->seqHap1, (void *)&chunk->seqHap2);
        swap((void *)&chunk->poaHap1StringsLines, (void *)&chunk->poaHap2StringsLines);
        swap((void *)&chunk->readsHap1Lines, (void *)&chunk->readsHap2Lines);
    }

    //TODO: Remove reads from conflicting partition output
}

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

    if (stList_length(alignedPairs) == 0 && st_getLogLevel() >= info) {
        st_logInfo("    Failed to find good overlap. Suffix-string: %s\n", &(prefixString[i]));
        st_logInfo("    Failed to find good overlap. Prefix-string: %s\n", suffixString);
    }

    // Remove the suffix crop
    suffixString[j] = c;

    // Pick the median point
    stIntTuple *maxPair = NULL;
    for (int64_t k = 0; k < stList_length(alignedPairs); k++) {
        stIntTuple *aPair = stList_get(alignedPairs, k);
        if (maxPair == NULL || stIntTuple_get(aPair, 0) > stIntTuple_get(maxPair, 0)) {
            maxPair = aPair;
        }
    }
    if (maxPair == NULL) {
        st_logCritical("    Failed to find any aligned pairs between overlapping strings, not "
                       "doing any trimming (approx overlap: %i, len x: %i, len y: %i)\n", approxOverlap,
                       prefixStringLength, suffixStringLength);
        *prefixStringCropEnd = prefixStringLength;
        *suffixStringCropStart = 0;
    } else {
        *prefixStringCropEnd = stIntTuple_get(maxPair, 1) + i; // Exclusive
        *suffixStringCropStart = stIntTuple_get(maxPair, 2);  // Inclusive
    }

    int64_t overlapWeight = maxPair == NULL ? 0 : stIntTuple_get(maxPair, 0);

    stList_destruct(alignedPairs);

    return overlapWeight;
}

void stList_removeInterval(stList *l, int64_t start, int64_t length) {

}

void chunkToStitch_trimAdjacentChunks2(char *pSeq, char **seq,
        stList *pPoa, stList *poa, stList *pRepeatCounts, stList *repeatCounts,
        Params *params) {
    // Get the trim factor
    int64_t prefixStringCropEnd, suffixStringCropStart;
    int64_t overlapMatchWeight = removeOverlap(pSeq, *seq, params->polishParams->chunkBoundary * 2,
            params->polishParams, &prefixStringCropEnd, &suffixStringCropStart);

    // Log
    st_logInfo("Removed overlap between neighbouring chunks. Approx overlap size: %i, overlap-match weight: %f, "
               "left-trim: %i, right-trim: %i:\n", (int) params->polishParams->chunkBoundary * 2,
               (float) overlapMatchWeight / PAIR_ALIGNMENT_PROB_1,
               strlen(pSeq) - prefixStringCropEnd, suffixStringCropStart);

    // Trim the sequences

    // Crop the suffix of the previous sequence
    pSeq[prefixStringCropEnd] = '\0';

    // Crop the the prefix of the current chunk's sequence
    char *c = *seq;
    *seq = stString_copy(&(*seq[suffixStringCropStart]));
    free(c);

    // Trim the remaining stuff

    // Poa
    stList_removeInterval(pPoa, prefixStringCropEnd, stList_length(pPoa)-prefixStringCropEnd);
    stList_removeInterval(poa, 1, suffixStringCropStart-1);

    // Repeat counts
    stList_removeInterval(pRepeatCounts, prefixStringCropEnd, stList_length(pPoa));
    stList_removeInterval(repeatCounts, 1, suffixStringCropStart);
}

void chunkToStitch_trimAdjacentChunks(ChunkToStitch *pChunk, ChunkToStitch *chunk, Params *params) {
    /*
     * Trims the right end of pChunk and the left end of chunk so that they do not overlap, but are directly contiguous.
     */
    // Checks that they are part of the same sequence
    assert(stString_eq(pChunk->seqName1, chunk->seqName1));
    assert(stString_eq(pChunk->seqName2, chunk->seqName2));

    // Trim haplotype 1 sequences
    chunkToStitch_trimAdjacentChunks2(pChunk->seqHap1, &chunk->seqHap1,
                                      pChunk->poaHap1StringsLines, chunk->poaHap1StringsLines,
                                      pChunk->repeatCountLines, chunk->repeatCountLines, params);

    // Trim haplotype 2 sequences
    chunkToStitch_trimAdjacentChunks2(pChunk->seqHap2, &chunk->seqHap2,
                                      pChunk->poaHap2StringsLines, chunk->poaHap2StringsLines, NULL, NULL, params);
}

/*
 * OutputChunker
 */

typedef struct _outputChunker {
    /*
     * Object for managing the output of a polished sequence.
     */
    // Sequence file
    char *outputSequenceFile;
    FILE *outputSequenceFileHandle;
    // Poa file
    char *outputPoaFile;
    FILE *outputPoaFileHandle;
    // Repeat count file
    char *outputRepeatCountFile;
    FILE *outputRepeatCountFileHandle;
    // Read partition file
    char *outputReadPartitionFile;
    FILE *outputReadPartitionFileHandle;

    Params *params;
} OutputChunker;

void outputChunker_open(OutputChunker *outputChunker, char *openStr) {
    /*
     * Open the files.
     */
    outputChunker->outputSequenceFileHandle = fopen(outputChunker->outputSequenceFile, openStr);
    outputChunker->outputPoaFileHandle = outputChunker->outputPoaFile != NULL ?
                                         fopen(outputChunker->outputPoaFile, openStr) : NULL;
    outputChunker->outputRepeatCountFileHandle = outputChunker->outputRepeatCountFile != NULL ?
                                                 fopen(outputChunker->outputRepeatCountFile, openStr) : NULL;
    outputChunker->outputReadPartitionFileHandle = outputChunker->outputReadPartitionFile != NULL ?
                                                   fopen(outputChunker->outputReadPartitionFile, openStr) : NULL;
}

OutputChunker *outputChunker_construct(Params *params, char *outputSequenceFile, char *outputPoaFile,
                                       char *outputReadPartitionFile, char *outputRepeatCountFile) {
    /*
     * Create an OutputChunker object, ready to write chunks of output to the given output files.
     */
    OutputChunker *outputChunker = st_calloc(1, sizeof(OutputChunker));

    // Initialize variables
    outputChunker->outputSequenceFile = outputSequenceFile;
    outputChunker->outputPoaFile = outputPoaFile;
    outputChunker->outputRepeatCountFile = outputRepeatCountFile;
    outputChunker->outputReadPartitionFile = outputReadPartitionFile;
    outputChunker->params = params;

    // Open files for writing
    outputChunker_open(outputChunker, "w");

    return outputChunker;
}

void
outputChunker_processChunkSequence(OutputChunker *outputChunker, int64_t chunkOrdinal, Poa *poa, BamChunk *bamChunk,
                                   stList *reads) {
    // Create chunk name
    char *headerLinePrefix = stString_print("%s\t%" PRIi64 "\t", bamChunk->refSeqName, chunkOrdinal);

    // Do run-length decoding
    char *outputSequence = rleString_expand(poa->refString);

    // Output the sequence
    fprintf(outputChunker->outputSequenceFileHandle, "%s1\n%s\n", headerLinePrefix, outputSequence);

    // Write any optional outputs about repeat count and POA, etc.

    // Poa
    if (outputChunker->outputPoaFileHandle != NULL) {
        fprintf(outputChunker->outputPoaFileHandle, "%s%" PRIi64 "\n", headerLinePrefix, stList_length(poa->nodes)+1);
        poa_printCSV(poa, outputChunker->outputPoaFileHandle, reads, outputChunker->params->polishParams->repeatSubMatrix, 5);
    }

    // Now repeat counts
    if (outputChunker->outputRepeatCountFileHandle != NULL) {
        fprintf(outputChunker->outputRepeatCountFileHandle, "%s%" PRIi64 "\n", headerLinePrefix, stList_length(poa->nodes)+1);
        poa_printRepeatCounts(poa, outputChunker->outputRepeatCountFileHandle, reads);
    }

    // Cleanup
    free(headerLinePrefix);
    free(outputSequence);
}

void
outputChunker_processChunkSequencePhased2(OutputChunker *outputChunker, char *headerLinePrefix,
                                          Poa *poa, stList *reads, stSet *readsBelongingToHap1, stSet *readsBelongingToHap2) {
    // Output the sequence
    char *outputSequence = rleString_expand(poa->refString);  // Do run-length decoding
    fprintf(outputChunker->outputSequenceFileHandle, "%s1\n%s\n", headerLinePrefix, outputSequence);
    free(outputSequence); // Cleanup

    // Write any optional outputs about repeat count and POA, etc.

    // Poa
    if (outputChunker->outputPoaFileHandle != NULL) {
        fprintf(outputChunker->outputPoaFileHandle, "%s%" PRIi64 "\n", headerLinePrefix, stList_length(poa->nodes)+1);
        poa_printPhasedCSV(poa, outputChunker->outputPoaFileHandle, reads, readsBelongingToHap1, readsBelongingToHap2,
                            outputChunker->params->polishParams->repeatSubMatrix, 5);
    }
}

void
outputChunker_processChunkSequencePhased(OutputChunker *outputChunker, int64_t chunkOrdinal, Poa *poaHap1, Poa *poaHap2,
                                         BamChunk *bamChunk, stList *reads, stSet *readsBelongingToHap1,
                                         stSet *readsBelongingToHap2, stGenomeFragment *gF) {
    // Create chunk name
    char *headerLinePrefix = stString_print("%s\t%" PRIi64 "\t", bamChunk->refSeqName, chunkOrdinal);

    outputChunker_processChunkSequencePhased2(outputChunker, headerLinePrefix,
                                            poaHap1, reads, readsBelongingToHap1, readsBelongingToHap2);

    outputChunker_processChunkSequencePhased2(outputChunker, headerLinePrefix,
                                              poaHap2, reads, readsBelongingToHap2, readsBelongingToHap1);

    // Output the read partition
    if (outputChunker->outputReadPartitionFileHandle != NULL) {
        fprintf(outputChunker->outputReadPartitionFileHandle, "%s%" PRIi64 "\n", headerLinePrefix,
                stSet_size(gF->reads1)+stSet_size(gF->reads2)+1);
        stGenomeFragment_printPartitionAsCSV(gF, outputChunker->outputReadPartitionFileHandle);
    }

    // Output the repeat counts
    if (outputChunker->outputRepeatCountFileHandle != NULL) {
        fprintf(outputChunker->outputRepeatCountFileHandle, "%s%" PRIi64 "\n", headerLinePrefix, stList_length(poaHap1->nodes)+1);
        poa_printRepeatCounts(poaHap1, outputChunker->outputRepeatCountFileHandle, reads);
    }

    // Cleanup
    free(headerLinePrefix);
}

ChunkToStitch *outputChunker_readChunk(OutputChunker *outputChunker, bool phased) {
    /*
     * Read a chunk of output from the outputChunker.
     */
    ChunkToStitch *chunk = st_calloc(1, sizeof(ChunkToStitch));

    if (!chunkToStitch_readSequenceChunk(outputChunker->outputSequenceFileHandle, chunk, phased)) {
        free(chunk);
        return NULL;
    }
    if (outputChunker->outputPoaFile != NULL) {
        chunkToStitch_readPoaChunk(outputChunker->outputPoaFileHandle, chunk, phased);
    }
    if (phased && outputChunker->outputReadPartitionFile != NULL) {
        chunkToStitch_readReadPhasingChunk(outputChunker->outputReadPartitionFileHandle, chunk);
    }
    if (outputChunker->outputRepeatCountFile != NULL) {
        chunkToStitch_readRepeatCountChunk(outputChunker->outputRepeatCountFileHandle, chunk);
    }

    return chunk;
}

void outputChunker_close(OutputChunker *rSeq) {
    // Cleanup the sequence output file
    if (rSeq->outputSequenceFileHandle != NULL) {
        fclose(rSeq->outputSequenceFileHandle);
        rSeq->outputSequenceFileHandle = NULL;
    }

    // Cleanup repeat count file
    if (rSeq->outputRepeatCountFileHandle != NULL) {
        fclose(rSeq->outputRepeatCountFileHandle);
        rSeq->outputRepeatCountFileHandle = NULL;
    }

    // Cleanup poa file
    if (rSeq->outputPoaFileHandle != NULL) {
        fclose(rSeq->outputPoaFileHandle);
        rSeq->outputPoaFileHandle = NULL;
    }

    // Cleanup read partition file
    if (rSeq->outputReadPartitionFileHandle != NULL) {
        fclose(rSeq->outputReadPartitionFileHandle);
        rSeq->outputReadPartitionFileHandle = NULL;
    }
}

void outputChunker_destruct(OutputChunker *rSeq) {
    // Close any open file descriptors
    outputChunker_close(rSeq);

    // Cleanup the sequence output file
    free(rSeq->outputSequenceFile);

    // Cleanup repeat count file
    if (rSeq->outputRepeatCountFile != NULL) {
        free(rSeq->outputRepeatCountFile);
    }

    // Cleanup poa file
    if (rSeq->outputPoaFile != NULL) {
        free(rSeq->outputPoaFile);
    }

    // Cleanup read partition file
    if (rSeq->outputReadPartitionFile != NULL) {
        free(rSeq->outputReadPartitionFile);
    }

    // Cleanup residual
    free(rSeq);
}

/*
 * OutputChunkers
 */

struct _outputChunkers {
    stList *tempFileChunkers;
    int64_t tempFileChunkerCounter;
    int64_t chunkOrderNo;
    OutputChunker *outputChunkerHap1;
    OutputChunker *outputChunkerHap2;
    Params *params;
};

static char *printFileName(char *fileName, int64_t index, char *suffix) {
    if(fileName == NULL) {
        return NULL;
    }
    return stString_print("%s%" PRIi64 "%s", fileName, index, suffix);
}

static char *printFileName2(char *fileName, char *suffix, char *fileSuffix) {
    if(fileName == NULL) {
        return NULL;
    }
    return stString_print("%s%s.%s", fileName, suffix, fileSuffix);
}

OutputChunkers *outputChunkers_construct(int64_t noOfOutputChunkers, Params *params,
                                         char *outputSequenceFile, char *outputPoaFile,
                                         char *outputReadPartitionFile, char *outputRepeatCountFile,
                                         char *hap1Suffix, char *hap2Suffix) {
    OutputChunkers *outputChunkers = st_calloc(1, sizeof(OutputChunkers));
    outputChunkers->params = params;
    // Make the temporary, parallel chunkers
    outputChunkers->tempFileChunkers = stList_construct3(0, (void (*)(void *)) outputChunker_destruct);
    for (int64_t i = 0; i < noOfOutputChunkers; i++) {
        stList_append(outputChunkers->tempFileChunkers, outputChunker_construct(params,
                                                                                printFileName(outputSequenceFile, i, ".temp.fa"),
                                                                                printFileName(outputPoaFile, i, ".temp.csv"),
                                                                                printFileName(outputReadPartitionFile, i, ".temp.csv"),
                                                                                printFileName(outputRepeatCountFile, i, ".temp.csv")));
    }
    // Make the final output chunkers
    outputChunkers->outputChunkerHap1 = outputChunker_construct(params,
                                                                printFileName2(outputSequenceFile, hap1Suffix, "fa"),
                                                                printFileName2(outputPoaFile, hap1Suffix, "csv"),
                                                                printFileName2(outputReadPartitionFile, hap1Suffix, "csv"),
                                                                printFileName2(outputRepeatCountFile, hap1Suffix, "csv"));

    if(hap2Suffix != NULL) {
        outputChunkers->outputChunkerHap2 = outputChunker_construct(params,
                                                                    printFileName2(outputSequenceFile, hap2Suffix,
                                                                                   "fa"),
                                                                    printFileName2(outputPoaFile, hap2Suffix, "csv"),
                                                                    printFileName2(outputReadPartitionFile, hap2Suffix,
                                                                                   "csv"),
                                                                    printFileName2(outputRepeatCountFile, hap2Suffix,
                                                                                   "csv"));
    }

    return outputChunkers;
}

void
outputChunkers_processChunkSequence(OutputChunkers *outputChunkers, int64_t chunker, int64_t chunkOrdinal, Poa *poa,
                                    stList *reads, BamChunk *bamChunk) {
    outputChunker_processChunkSequence(stList_get(outputChunkers->tempFileChunkers, chunker), chunkOrdinal, poa, bamChunk,
                                       reads);
}

void outputChunkers_processChunkSequencePhased(OutputChunkers *outputChunkers, int64_t chunker, int64_t chunkOrdinal,
                                               Poa *poaHap1, Poa *poaHap2, BamChunk *bamChunk, stList *reads,
                                               stSet *readsBelongingToHap1, stSet *readsBelongingToHap2,
                                               stGenomeFragment *gF) {
    outputChunker_processChunkSequencePhased(stList_get(outputChunkers->tempFileChunkers, chunker), chunkOrdinal,
                                             poaHap1, poaHap2, bamChunk, reads, readsBelongingToHap1,
                                             readsBelongingToHap2, gF);
}


void outputChunkers_close(OutputChunkers *outputChunkers) {
    /*
     * Closes the file handles used.
     */
    for (int64_t i = 0; i < stList_length(outputChunkers->tempFileChunkers); i++) {
        outputChunker_close(stList_get(outputChunkers->tempFileChunkers, i));
    }
    outputChunker_close(outputChunkers->outputChunkerHap1);
    if(outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_close(outputChunkers->outputChunkerHap2);
    }
}

void outputChunkers_openForStitching(OutputChunkers *outputChunkers) {
    /*
     * Open the output chunkers for stitching
     */
    outputChunkers_close(outputChunkers);
    for (int64_t i = 0; i < stList_length(outputChunkers->tempFileChunkers); i++) {
        outputChunker_open(stList_get(outputChunkers->tempFileChunkers, i), "r");
    }
    outputChunker_open(outputChunkers->outputChunkerHap1, "w");
    if(outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_open(outputChunkers->outputChunkerHap2, "w");
    }
}

static ChunkToStitch *outputChunkers_readChunk(OutputChunkers *outputChunkers, bool phased) {
    for (int64_t i = 0; i < stList_length(outputChunkers->tempFileChunkers); i++) {
        ChunkToStitch *chunk = outputChunker_readChunk(stList_get(outputChunkers->tempFileChunkers,
                                                                   outputChunkers->tempFileChunkerCounter++ %
                                                                   stList_length(outputChunkers->tempFileChunkers)), phased);
        if (chunk != NULL) {
            return chunk;
        }
    }
    return NULL;
}

static ChunkToStitch *outputChunkers_getNextChunkInSequence(OutputChunkers *outputChunkers,
                                                     stSortedSet *orderedChunks, bool phased) {
    while (1) {
        ChunkToStitch *chunk = outputChunkers_readChunk(outputChunkers, phased);
        if (chunk == NULL) {
            if (stSortedSet_size(orderedChunks) != 0) {
                st_errAbort("Ran out of chunks while still having some out");
            }
            return NULL;
        }
        stSortedSet_insert(orderedChunks, chunk);
        chunk = stSortedSet_getFirst(orderedChunks);
        if (chunk->chunkOrdinal == outputChunkers->chunkOrderNo) {
            stSortedSet_remove(orderedChunks, chunk);
            outputChunkers->chunkOrderNo++;
            return chunk;
        }
    }
}

void writeLines(FILE *fh, stList *lines) {
    for(int64_t i=0; i<stList_length(lines); i++) {
        fprintf(fh, "%s\n", stList_get(lines, i));
    }
}

void outputChunker_writeChunkToFinalOutput(OutputChunker *outputChunker,
        char *seqName, char *seq, stList *poaLines, stList *repeatCountLines,
        stList *readPartitionLines, bool startOfSequence) {
    /*
     * Writes the chunk to the final output files
     */

    // Write the sequence
    if(startOfSequence) {
        fastaWrite(seq, seqName, outputChunker->outputSequenceFileHandle);
    }
    else {
        fprintf(outputChunker->outputSequenceFileHandle, "%s\n", seq);
    }

    // Write the POA
    if(outputChunker->outputPoaFile != NULL) {
        writeLines(outputChunker->outputPoaFileHandle, poaLines);
    }

    // Write the repeat counts
    if(outputChunker->outputRepeatCountFile != NULL && repeatCountLines != NULL) {
        writeLines(outputChunker->outputRepeatCountFileHandle, repeatCountLines);
    }

    // Write the read
    if(outputChunker->outputReadPartitionFile != NULL) {
        writeLines(outputChunker->outputReadPartitionFileHandle, readPartitionLines);
    }
}

void outputChunkers_writeChunk(OutputChunkers *outputChunkers, ChunkToStitch *chunk) {
    /*
     * Writes the chunk to the final output files
     */
    outputChunker_writeChunkToFinalOutput(outputChunkers->outputChunkerHap1,
                                        chunk->seqName1, chunk->seqHap1, chunk->poaHap1StringsLines,
                                        chunk->repeatCountLines, chunk->readsHap1Lines, chunk->startOfSequence);
    if(outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_writeChunkToFinalOutput(outputChunkers->outputChunkerHap2,
                                              chunk->seqName2, chunk->seqHap2, chunk->poaHap2StringsLines,
                                              NULL, chunk->readsHap2Lines, chunk->startOfSequence);
    }
}

void outputChunkers_stitch(OutputChunkers *outputChunkers, bool phased) {
    /*
     * Stitch together the outputs
     */

    // Setup to write out the chunks
    outputChunkers_openForStitching(outputChunkers);

    // Create a cache to hold the chunks, ordered by their ordinal
    stSortedSet *orderedChunks = stSortedSet_construct3((int (*)(const void *, const void *)) chunkToStitch_cmp, NULL);

    // Get the first chunk
    ChunkToStitch *pChunk = outputChunkers_getNextChunkInSequence(outputChunkers, orderedChunks, phased), *chunk;

    // Set the start of sequence flag
    if(pChunk != NULL) {
        pChunk->startOfSequence = 1;
    }

    // Get each successive chunk and stitch and phase progressively
    while ((chunk = outputChunkers_getNextChunkInSequence(outputChunkers, orderedChunks, phased)) != NULL) {
        assert(pChunk != NULL);

        // If phased, ensure the chunks phasing is consistent
        if (phased) {
            chunkToStitch_phaseAdjacentChunks(pChunk, chunk);
        }

        // Set the flag determining if this is the start of a new sequence
        chunk->startOfSequence = stString_eq(pChunk->seqName1, chunk->seqName1);

        // Trim the boundaries of the two chunks so that they don't overlap
        if(!chunk->startOfSequence) {
            chunkToStitch_trimAdjacentChunks(pChunk, chunk, outputChunkers->params);
        }

        // Write out the pChunk
        outputChunkers_writeChunk(outputChunkers, pChunk);

        // Cleanup the pChunk
        chunkToStitch_destruct(pChunk);

        // Set the new pChunk
        pChunk = chunk;
    }

    if (pChunk != NULL) {
        // Write out the pChunk
        outputChunkers_writeChunk(outputChunkers, pChunk);

        // Cleanup the pChunk
        chunkToStitch_destruct(pChunk);
    }

    // Cleanup
    if(stSortedSet_size(orderedChunks) != 0) {
        st_errAbort("Got chunks left over after writing out chunks");
    }
    stSortedSet_destruct(orderedChunks);
}

void outputChunkers_destruct(OutputChunkers *outputChunkers) {
    // Close file streams
    outputChunkers_close(outputChunkers);


    stList_destruct(outputChunkers->tempFileChunkers);
    outputChunker_destruct(outputChunkers->outputChunkerHap1);
    if(outputChunkers->outputChunkerHap2 != NULL) {
        outputChunker_destruct(outputChunkers->outputChunkerHap2);
    }
    free(outputChunkers);
}

