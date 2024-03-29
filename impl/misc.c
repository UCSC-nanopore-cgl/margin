/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "margin.h"
#include <sys/stat.h>
#include <sonLibListPrivate.h>
#include <helenFeatures.h>
#include <htsIntegration.h>

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
    FILE *fh = safe_fopen(referenceFastaFile, "r");
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
    if (access( base, F_OK ) != -1  && S_ISDIR(fileStat.st_mode)) {
        if (optarg[strlen(base) - 1] == '/') optarg[strlen(base) - 1] = '\0';
        return stString_print("%s/%s", base, defawlt);
    } else {
        return stString_copy(base);
    }
}

FILE *safe_fopen(char *file, char *openStr) {
    FILE *fp = fopen(file, openStr);
    if (fp == NULL) {
        st_logInfo("> File %s failed to open! Attempting again. \n", file);
        usleep(100000); // 100ms
        fp = fopen(file, openStr);
        if (fp == NULL) {
            st_errAbort("Could not open file %s for '%s'\n", file, openStr);
        }
    }
    return fp;
}

RleString *bamChunk_getReferenceSubstring(BamChunk *bamChunk, char *referenceFile, Params *params) {
    /*
     * Get corresponding substring of the reference for a given bamChunk.
     */
    char *referenceString = getSequenceFromReference(referenceFile, bamChunk->refSeqName, bamChunk->chunkOverlapStart,
                                                     bamChunk->chunkOverlapEnd);
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
    }
    return copy;
}

double fromLog(double theLog) {
    return pow(10.0, theLog);
}

double toPhred(double prob) {
    return -10 * log10(prob <= 0.1 ? 0.000001 : prob >= 0.999999 ? 0.999999 : prob);
}

double fromPhred(double phred) {
    return pow(10, phred / -10.0);
}


void removeReadsStartingAfterChunkEnd(BamChunk *bamChunk, stList *reads, stList *alignments, char *logIdentifier) {
    int64_t chunkEnd = bamChunk->chunkEnd - bamChunk->chunkOverlapStart;
    stList *indicesToRemove = stList_construct();
    for (int64_t i = stList_length(reads) - 1; i >= 0; i--) {
        stList *alignment = stList_get(alignments, i);
        if (stList_length(alignment) == 0) continue;
        int64_t firstAlignPos = stIntTuple_get(stList_get(alignment, 0), 0);
        if (firstAlignPos >= chunkEnd) {
            stList_append(indicesToRemove, (void*) i);
        }
    }
    st_logInfo(" %s Removing %"PRIu64" of %"PRIu64" reads for starting after chunkEnd\n", logIdentifier,
               stList_length(indicesToRemove), stList_length(reads));
    for (int64_t i = 0; i < stList_length(indicesToRemove); i++) {
        int64_t indexToRemove = (int64_t) stList_get(indicesToRemove, i);
        bamChunkRead_destruct(stList_get(reads, indexToRemove));
        stList_destruct(stList_get(alignments, indexToRemove));
        stList_remove(reads, indexToRemove);
        stList_remove(alignments, indexToRemove);
    }
    stList_destruct(indicesToRemove);
}

void removeReadsOnlyInChunkBoundary(BamChunk *bamChunk, stList *reads, stList *alignments, char *logIdentifier) {
    int64_t chunkEnd = bamChunk->chunkEnd - bamChunk->chunkOverlapStart;
    int64_t chunkStart = bamChunk->chunkStart - bamChunk->chunkOverlapStart;
    stList *indicesToRemove = stList_construct();
    for (int64_t i = stList_length(reads) - 1; i >= 0; i--) {
        stList *alignment = stList_get(alignments, i);
        if (stList_length(alignment) == 0) continue;
        int64_t firstAlignPos = stIntTuple_get(stList_get(alignment, 0), 0);
        int64_t lastAlignPos = stIntTuple_get(stList_get(alignment, stList_length(alignment) - 1), 0);
        if (lastAlignPos < chunkStart || firstAlignPos >= chunkEnd) {
            stList_append(indicesToRemove, (void*) i);
        }
    }
    st_logInfo(" %s Removing %"PRIu64" of %"PRIu64" reads for falling outside chunk boundaries\n", logIdentifier,
               stList_length(indicesToRemove), stList_length(reads));
    for (int64_t i = 0; i < stList_length(indicesToRemove); i++) {
        int64_t indexToRemove = (int64_t) stList_get(indicesToRemove, i);
        bamChunkRead_destruct(stList_get(reads, indexToRemove));
        stList_destruct(stList_get(alignments, indexToRemove));
        stList_remove(reads, indexToRemove);
        stList_remove(alignments, indexToRemove);
    }
    stList_destruct(indicesToRemove);
}

void writePhasedReadInfoJSON(BamChunk *bamChunk, stList *primaryReads, stList *primaryAlignments, stList *filteredReads,
        stList *filteredAlignments, stSet *readsInHap1, stSet *readsInHap2, uint64_t *reference_rleToNonRleCoordMap,
        FILE *out) {


    fprintf(out, ",\n \"reads\": [");

    // get scores and save to appropriate sets
    for (int64_t i = 0; i < stList_length(primaryReads); i++) {
        BamChunkRead *read = stList_get(primaryReads, i);
        stList *alignment = stList_get(primaryAlignments, i);
        stIntTuple *firstAlign = stList_get(alignment, 0);
        stIntTuple *lastAlign = stList_get(alignment, stList_length(alignment) - 1);
        int hap = 0;
        if (stSet_search(readsInHap1, read)) {
            hap = 1;
        } else if (stSet_search(readsInHap2, read)) {
            hap = 2;
        }

        //write
        if (i != 0) {
            fprintf(out, ",");
        }
        fprintf(out, "\n  {\n");
        fprintf(out, "     \"name\": \"%s\",\n", read->readName);
        fprintf(out, "     \"strand\": \"%s\",\n", read->forwardStrand ? "+" : "-");
        fprintf(out, "     \"startPos\": %"PRId64",\n", bamChunk->chunkOverlapStart + reference_rleToNonRleCoordMap[stIntTuple_get(firstAlign, 0)]);
        fprintf(out, "     \"endPos\": %"PRId64",\n", bamChunk->chunkOverlapStart + reference_rleToNonRleCoordMap[stIntTuple_get(lastAlign, 0)]);
        fprintf(out, "     \"hap\": %d\n", hap);
        fprintf(out, "  }");
    }

    // filtered reads
    for (int64_t i = 0; i < stList_length(filteredReads); i++) {
        BamChunkRead *read = stList_get(filteredReads, i);
        stList *alignment = stList_get(filteredAlignments, i);
        stIntTuple *firstAlign = stList_get(alignment, 0);
        stIntTuple *lastAlign = stList_get(alignment, stList_length(alignment) - 1);
        int hap = 0;
        if (stSet_search(readsInHap1, read)) {
            hap = 1;
        } else if (stSet_search(readsInHap2, read)) {
            hap = 2;
        }

        //write
        fprintf(out, ",\n  {\n");
        fprintf(out, "     \"name\": \"%s\",\n", read->readName);
        fprintf(out, "     \"strand\": \"%s\",\n", read->forwardStrand ? "+" : "-");
        fprintf(out, "     \"startPos\": %"PRId64",\n", bamChunk->chunkOverlapStart + reference_rleToNonRleCoordMap[stIntTuple_get(firstAlign, 0)]);
        fprintf(out, "     \"endPos\": %"PRId64",\n", bamChunk->chunkOverlapStart + reference_rleToNonRleCoordMap[stIntTuple_get(lastAlign, 0)]);
        fprintf(out, "     \"hap\": %d\n", hap);
        fprintf(out, "  }");
    }

    // write
    if (out != NULL) {
        fprintf(out, "\n ]");
    }
}


stList *produceVcfEntriesFromBubbleGraph(BamChunk *bamChunk, BubbleGraph *bg, stHash *readsToPSeqs,
                                         stGenomeFragment *gF, double strandSkewThreshold,
                                         double readSkewThreshold) {

    stList *vcfEntries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);

    char *logIdentifier = getLogIdentifier();
    int64_t passes = 0;
    int64_t total = 0;
    int64_t failsStrandSkew = 0;
    int64_t failsReadSkew = 0;
    for (uint64_t i = 0; i < gF->length; i++) {

        // bubble and hap info
        Bubble *b = &bg->bubbles[gF->refStart + i];
        RleString *hap1 = b->alleles[gF->haplotypeString1[i]];
        RleString *hap2 = b->alleles[gF->haplotypeString2[i]];

        // we only care about hets
        if (hap1 == hap2) continue;

        // read info
        int64_t hap1AlleleNo = gF->haplotypeString1[i];
        int64_t hap2AlleleNo = gF->haplotypeString2[i];
        int64_t totalReads = 0;
        int64_t hap1Reads = 0;
        for (uint64_t j = 0; j < b->readNo; j++) {

            double readHap1Support = b->alleleReadSupports[hap1AlleleNo * b->readNo + j];
            double readHap2Support = b->alleleReadSupports[hap2AlleleNo * b->readNo + j];

            if (readHap1Support != readHap2Support) {
                totalReads++;
                if (readHap1Support > readHap2Support) {
                    hap1Reads++;
                }
            }
        }

        // bubble skew info
        double strandSkew = bubble_phasedStrandSkew(b, readsToPSeqs, gF);
        double readSupportSkew = binomialPValue(totalReads, hap1Reads); // 1-bpv = prob of having a less extreme

        bool pass = TRUE;
        if (strandSkew < strandSkewThreshold) {
            failsStrandSkew++;
            pass = FALSE;
        }
        if (readSupportSkew < readSkewThreshold) {
            failsReadSkew++;
            pass = FALSE;
        }

        // maybe save to vcf entry
        uint64_t bubblePos = (uint64_t ) (b->refStart + b->bubbleLength / 2);
        if (pass) {
            for (int64_t cvp = 0; cvp < stList_length(b->variantPositionOffsets); cvp++) {
                stList *alleles = stList_construct();
                stList_append(alleles, rleString_copy(b->refAllele));
                int hap1Allele, hap2Allele;
                if (b->refAllele == hap1) {
                    hap1Allele = 0;
                } else {
                    hap1Allele = 1;
                    stList_append(alleles, rleString_copy(hap1));
                }
                if (b->refAllele == hap2) {
                    hap2Allele = 0;
                } else {
                    hap2Allele = hap1Allele + 1;
                    stList_append(alleles, rleString_copy(hap2));
                }
                stList_append(vcfEntries, vcfEntry_construct(bamChunk->refSeqName,
                        b->refStart + (int64_t) stList_get(b->variantPositionOffsets, cvp), -1, -1,
                        hap1->nonRleLength != hap2->nonRleLength, FALSE, alleles, hap1Allele, hap2Allele));
            }
            passes++;
        }
        total++;

        //loggit
        st_logDebug(" %s Bubble at %"PRIu64" (%"PRIu64"+%"PRIu64") has strand skew %.5f%s and read support skew %.5f%s (%2"PRId64":%-2"PRId64" ): %s\n",
                logIdentifier, bubblePos, b->refStart, b->bubbleLength, strandSkew, strandSkew < strandSkewThreshold ? "*" : " ",
                readSupportSkew, readSupportSkew < readSkewThreshold ? "*" : " ", hap1Reads, totalReads-hap1Reads,
                pass ? "PASS" : "FAIL");
    }

    //loggit
    st_logInfo(" %s Kept %"PRId64" of %"PRId64" bubbles after quality filtering. Failures: %"PRId64" strand, %"PRId64" reads\n",
            logIdentifier, passes, total, failsStrandSkew, failsReadSkew);
    free(logIdentifier);
    return vcfEntries;
}

ChunkTruthHaplotypes **chunkTruthHaplotypes_construct(int64_t length) {
    ChunkTruthHaplotypes **cths = st_calloc(length, sizeof(ChunkTruthHaplotypes *));
    for (int i = 0; i < length; i++) {
        ChunkTruthHaplotypes *cth = st_calloc(1, sizeof(ChunkTruthHaplotypes));
        cth->chunkIdx = i;
        cth->hap1Reads = stList_construct3(0, free);
        cth->hap2Reads = stList_construct3(0, free);
        cths[i] = cth;
    }
    return cths;
}
void chunkTruthHaplotypes_destruct(ChunkTruthHaplotypes **chunkTruthHaplotypes, int64_t length) {
    for (int i = 0; i < length; i++) {
        ChunkTruthHaplotypes *cth = chunkTruthHaplotypes[i];
        stList_destruct(cth->hap1Reads);
        stList_destruct(cth->hap2Reads);
        free(cth);
    }
    free(chunkTruthHaplotypes);
}
bool chunkTruthHaplotypes_isChunkTruthRead(char *readId) {
    if (strlen(readId) <= CHUNK_TRUTH_READ_ID_LEN) return FALSE;
    for (int64_t i = 0; i < CHUNK_TRUTH_READ_ID_LEN; i++) {
        if (CHUNK_TRUTH_READ_ID[i] != readId[i]) return FALSE;
    }
    if (readId[CHUNK_TRUTH_READ_ID_LEN] != CHUNK_TRUTH_READ_ID_SEP[0]) return FALSE;
    return TRUE;
}
void *chunkTruthHaplotypes_print(stList *readsInHap1, stList *readsInHap2, stList *bamChunks, int64_t length, char *filename) {
    ChunkTruthHaplotypes **cths = chunkTruthHaplotypes_construct(length);

    // get reads
    for (int hap = 1; hap <= 2; hap++) {
        stList *readNames = hap == 1 ? readsInHap1 : readsInHap2;
        stListIterator *itor = stList_getIterator(readNames);
        char *cthReadName = NULL;
        while ((cthReadName = stList_getNext(itor)) != NULL) {
            if (!chunkTruthHaplotypes_isChunkTruthRead(cthReadName)) continue;

            // get parts
            stList *parts = stString_splitByString(cthReadName, CHUNK_TRUTH_READ_ID_SEP);
            assert(stList_length(parts) >= 3);

            // get chunk idx
            char *ptr;
            int64_t chunkIdx = strtol(stList_get(parts, 1), &ptr, 10);

            // get read name
            stList_removeInterval(parts, 0, 2);
            char *readName = stString_join2(CHUNK_TRUTH_READ_ID_SEP, parts);

            // save it
            ChunkTruthHaplotypes *cth = cths[chunkIdx];
            stList_append((hap == 1 ? cth->hap1Reads : cth->hap2Reads), readName);

            // cleanup
            stList_destruct(parts);
        }
        stList_destructIterator(itor);
    }

    // write
    FILE *out = safe_fopen(filename, "w");
    fprintf(out, "#contig\tstartPos\tendPos\toverlapStart\toverlapEnd\thap\tsequenceName\n");
    for (int64_t chunkIdx = 0; chunkIdx < length; chunkIdx++) {
        ChunkTruthHaplotypes *cth = cths[chunkIdx];
        BamChunk *bc = stList_get(bamChunks, chunkIdx);
        for (int64_t readIdx = 0; readIdx < stList_length(cth->hap1Reads); readIdx++) {
            fprintf(out, "%s\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t1\t%s\n", bc->refSeqName, bc->chunkStart,
                    bc->chunkEnd, bc->chunkOverlapStart, bc->chunkOverlapEnd,
                    (char*) stList_get(cth->hap1Reads, readIdx));
        }
        for (int64_t readIdx = 0; readIdx < stList_length(cth->hap2Reads); readIdx++) {
            fprintf(out, "%s\t%"PRId64"\t%"PRId64"\t%"PRId64"\t%"PRId64"\t2\t%s\n", bc->refSeqName, bc->chunkStart,
                    bc->chunkEnd, bc->chunkOverlapStart, bc->chunkOverlapEnd,
                    (char*) stList_get(cth->hap2Reads, readIdx));
        }
    }

    // cleanup
    fclose(out);
    chunkTruthHaplotypes_destruct(cths, length);
}


void chunkTruthHaplotypes_addTruthReadsToFilteredReadSet(BamChunk *bamChunk, BamChunker *bamChunker,
        stList *readsToAdd, stList *alignmentsToAdd, RleString *rleReference, Params *params, char *logIdentifier) {

    stList *reads = stList_construct();
    stList *alignments = stList_construct();
    BamChunker *tmp = bamChunk->parent;
    bamChunk->parent = bamChunker;

    // we want supplementary truth aligments
    PolishParams paramsCopy;
    memcpy(&paramsCopy, params->polishParams, sizeof(PolishParams));
    paramsCopy.includeSupplementaryAlignments = TRUE;

    // get truth "reads"
    convertToReadsAndAlignments(bamChunk, rleReference, reads, alignments, &paramsCopy);
    bamChunk->parent = tmp;
    st_logInfo(" %s Saving %"PRId64" truth reads from file: %s\n", logIdentifier, stList_length(reads),
               bamChunker->bamFile);

    // rename reads for later retrieval
    for (int64_t i = 0; i < stList_length(reads); i++) {
        BamChunkRead *bcr = stList_get(reads, i);

        // just giving them a new name should make the rest of the infrastructure work
        char *newName = stString_print("%s%s%"PRId64"%s%s", CHUNK_TRUTH_READ_ID, CHUNK_TRUTH_READ_ID_SEP,
                                       bamChunk->chunkIdx, CHUNK_TRUTH_READ_ID_SEP, bcr->readName);
        free(bcr->readName);
        bcr->readName = newName;
        stList_append(readsToAdd, bcr);
        stList_append(alignmentsToAdd, stList_get(alignments, i));
    }

    // cleanup
    stList_destruct(reads);
    stList_destruct(alignments);
}
