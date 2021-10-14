/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <getopt.h>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>
#include <hashTableC.h>
#include <unistd.h>
#include <time.h>
#include "marginVersion.h"

#include "margin.h"
#include "htsIntegration.h"
#include "helenFeatures.h"


/*
 * Main functions
 */

void usage() {
    fprintf(stderr, "usage: tagFromPhasedVcf <ALIGN_BAM> <REFERENCE_FASTA> <VARIANT_VCF> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Tags reads in ALIGN_BAM using already-phased variants in VARIANT_VCF.\n");
    fprintf(stderr, "    This tool does not attempt to phase variants!\n");
    fprintf(stderr, "    This tool does not have phaseset awareness, so the VCF should avoid overlapping phasesets\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    ALIGN_BAM is the alignment of reads to the reference.\n");
    fprintf(stderr, "    REFERENCE_FASTA is the reference sequence BAM file in fasta format.\n");
    fprintf(stderr, "    VARIANT_VCF is the set of variants to use for phasing.\n");
    fprintf(stderr, "    PARAMS is the file with margin parameters.\n");

    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help                : Print this help screen\n");
    fprintf(stderr, "    -a --logLevel            : Set the log level [default = info]\n");
# ifdef _OPENMP
    fprintf(stderr, "    -t --threads             : Set number of concurrent threads [default = 1]\n");
#endif
    fprintf(stderr, "    -o --outputBase          : Name to use for output files [default = 'output']\n");
    fprintf(stderr, "    -r --region              : If set, will only compute for given chromosomal region\n");
    fprintf(stderr, "                                 Format: chr:start_pos-end_pos (chr3:2000-3000)\n");
    fprintf(stderr, "    -k --tempFilesToDisk     : Write temporary files to disk (for --diploid or supplementary output)\n");

    fprintf(stderr, "\n");
}


int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *logLevelString = stString_copy("critical");
    char *bamInFile = NULL;
    char *paramsFile = NULL;
    char *referenceFastaFile = NULL;
    char *outputBase = stString_copy("output");
    char *regionStr = NULL;
    char *vcfFile = NULL;
    int numThreads = 1;
    bool inMemory = TRUE;

    if (argc < 4) {
        free(outputBase);
        free(logLevelString);
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    referenceFastaFile = stString_copy(argv[2]);
    vcfFile = stString_copy(argv[3]);
    paramsFile = stString_copy(argv[4]);

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "help", no_argument, 0, 'h' },
                { "logLevel", required_argument, 0, 'a' },
# ifdef _OPENMP
                { "threads", required_argument, 0, 't'},
#endif
                { "outputBase", required_argument, 0, 'o'},
                { "region", required_argument, 0, 'r'},
                { "tempFilesToDisk", no_argument, 0, 'k'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "ha:t:o:r:k", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'a':
            free(logLevelString);
            logLevelString = stString_copy(optarg);
            break;
        case 'h':
            usage();
            return 0;
        case 'o':
            free(outputBase);
            outputBase = getFileBase(optarg, "output");
            break;
        case 'r':
            regionStr = stString_copy(optarg);
            break;
        case 't':
            numThreads = atoi(optarg);
            if (numThreads <= 0) {
                st_errAbort("Invalid thread count: %d", numThreads);
            }
            break;
        case 'k':
            inMemory = FALSE;
            break;
        default:
            usage();
            free(outputBase);
            free(logLevelString);
            free(bamInFile);
            free(referenceFastaFile);
            free(paramsFile);
            return 0;
        }
    }

    // sanity check (verify files exist)
    if (access(bamInFile, R_OK) != 0) {
        st_errAbort("Could not read from input bam file: %s\n", bamInFile);
        char *idx = stString_print("%s.bai", bamInFile);
        if (access(idx, R_OK) != 0) {
            st_errAbort("BAM does not appear to be indexed: %s\n", bamInFile);
        }
        free(idx);
    }
    if (access(referenceFastaFile, R_OK) != 0) {
        st_errAbort("Could not read from reference fastafile: %s\n", referenceFastaFile);
    }
    if (access(vcfFile, R_OK) != 0) {
        st_errAbort("Could not read from vcf file: %s\n", vcfFile);
    }
    if (access(paramsFile, R_OK) != 0) {
        st_errAbort("Could not read from params file: %s\n", paramsFile);
    }

    // Initialization from arguments
    time_t startTime = time(NULL);
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);
    if (st_getLogLevel() >= info) {
        st_setCallocDebug(true);
    }
# ifdef _OPENMP
    if (numThreads <= 0) {
        numThreads = 1;
    }
    omp_set_num_threads(numThreads);
    st_logCritical("Running OpenMP with %d threads.\n", omp_get_max_threads());
    # endif

    // Parse parameters
    st_logCritical("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // Print a report of the parsed parameters
    if (st_getLogLevel() == debug) {
        params_printParameters(params, stderr);
    }

    // get vcf entries (if set)
    stHash *vcfEntries = NULL;
    if (vcfFile != NULL) {
        vcfEntries = parseVcf2(vcfFile, regionStr, params);
    }

    // get valid contigs (to help bam chunker construction)
    stList *vcfContigsTmp = stHash_getKeys(vcfEntries);
    stSet *vcfContigs = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    for (int64_t i = 0; i < stList_length(vcfContigsTmp); i++) {
        stSet_insert(vcfContigs, stList_get(vcfContigsTmp, i));
    }

    // get chunker for bam.  if regionStr is NULL, it will be ignored
    time_t chunkingStart = time(NULL);
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, vcfContigs, params->polishParams, TRUE);
    char *regionStrInformative = regionStr != NULL ? stString_copy(regionStr) : stString_join2(",", vcfContigsTmp);
    st_logCritical(
            "> Set up bam chunker in %"PRId64"s with chunk size %i and overlap %i (for region=%s), resulting in %i total chunks\n",
            time(NULL) - chunkingStart, (int) bamChunker->chunkSize, (int) bamChunker->chunkBoundary,
            regionStrInformative, bamChunker->chunkCount);
    if (bamChunker->chunkCount == 0) {
        st_errAbort("> Found no valid reads!\n");
    }
    free(regionStrInformative);
    stList_destruct(vcfContigsTmp);
    stSet_destruct(vcfContigs);

    // output chunker tracks intermediate output files
    OutputChunkers *outputChunkers = outputChunkers_construct(numThreads, params, NULL, NULL, NULL, NULL,
            ".hap1", ".hap2", inMemory);

    // (may) need to shuffle chunks
    stList *chunkOrder = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        stList_append(chunkOrder, stIntTuple_construct1(i));
    }
    if (params->polishParams->shuffleChunks) {
        switch (params->polishParams->shuffleChunksMethod) {
            case SCM_SIZE_DESC:
                st_logCritical("> Ordering chunks by estimated depth\n");
                stList_sort2(chunkOrder, compareBamChunkDepthByIndexInList, bamChunker->chunks);
                stList_reverse(chunkOrder);
                break;
            case SCM_RANDOM:
                st_logCritical("> Randomly shuffling chunks\n");
                stList_shuffle(chunkOrder);
                break;
        }
    }

    // multiproccess the chunks, save to results
    st_logCritical("> Setup complete, beginning run\n");
    int64_t lastReportedPercentage = 0;
    time_t polishStartTime = time(NULL);

    # ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    # endif
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        int64_t chunkIdx = stIntTuple_get(stList_get(chunkOrder, i), 0);
        // Time all chunks
        time_t chunkStartTime = time(NULL);

        // Get chunk
        BamChunk *bamChunk = bamChunker_getChunk(bamChunker, chunkIdx);

        // logging
        char *logIdentifier;
        bool logProgress = FALSE;
        int64_t currentPercentage = (int64_t) (100 * i / bamChunker->chunkCount);
        # ifdef _OPENMP
        int64_t threadIdx = omp_get_thread_num();
        logIdentifier = stString_print(" T%02d_C%05"PRId64, threadIdx, chunkIdx);
        if (threadIdx == 0) {
            if (currentPercentage != lastReportedPercentage) {
                logProgress = TRUE;
                lastReportedPercentage = currentPercentage;
            }
        }
        # else
        int64_t threadIdx = 0;
        logIdentifier = stString_copy("");
        if (currentPercentage != lastReportedPercentage) {
            logProgress = TRUE;
            lastReportedPercentage = currentPercentage;
        }
        # endif

        // prints percentage complete and estimated time remaining
        if (logProgress) {
            // log progress
            int64_t timeTaken = (int64_t) (time(NULL) - polishStartTime);
            int64_t secondsRemaining = (int64_t) floor(1.0 * timeTaken / currentPercentage * (100 - currentPercentage));
            char *timeDescriptor = (secondsRemaining == 0 && currentPercentage <= 50 ?
                                    stString_print("unknown") : getTimeDescriptorFromSeconds(secondsRemaining));
            st_logCritical("> Polishing %2"PRId64"%% complete (%"PRId64"/%"PRId64").  Estimated time remaining: %s\n",
                           currentPercentage, i, bamChunker->chunkCount, timeDescriptor);
            free(timeDescriptor);
        }

        // Get reference string for chunk of alignment
        char *chunkReference = getSequenceFromReference(referenceFastaFile, bamChunk->refSeqName,
                bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);
        st_logInfo(">%s Going to process a chunk for reference sequence: %s, starting at: %i and ending at: %i\n",
                   logIdentifier, bamChunk->refSeqName, (int) bamChunk->chunkOverlapStart, bamChunk->chunkOverlapEnd);

        // get VCF string
        stList *chunkVcfEntries = getVcfEntriesForRegion(vcfEntries, NULL, bamChunk->refSeqName,
                                                 bamChunk->chunkOverlapStart,  bamChunk->chunkOverlapEnd, params);
        updateVcfEntriesWithSubstringsAndPositions(chunkVcfEntries, chunkReference, strlen(chunkReference),
                FALSE, params);

        // Convert bam lines into corresponding reads and alignments
        st_logInfo(" %s Parsing input reads from file: %s\n", logIdentifier, bamInFile);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        extractReadSubstringsAtVariantPositions(bamChunk, chunkVcfEntries, reads, NULL, params->polishParams);

        time_t primaryPhasingStart = time(NULL);

        // iteratively find bubbles
        int64_t bubbleFindingIteration = 0;
        BubbleGraph *bg = NULL;
        stSet *readsBelongingToHap1 = stSet_construct(), *readsBelongingToHap2 = stSet_construct();
        stList *vcfEntriesToBubbles = NULL;

        // Get the bubble graph representation
        bg =  bubbleGraph_constructFromVCFAndBamChunkReadVcfEntrySubstrings(reads, chunkVcfEntries, params,
                &vcfEntriesToBubbles);

        bubbleGraph_partitionFilteredReadsFromPhasedVcfEntries(reads, bg, vcfEntriesToBubbles, readsBelongingToHap1,
                readsBelongingToHap2, params, logIdentifier);

        // Output
        stGenomeFragment *gF = stGenomeFragment_constructEmpty(NULL, 0, 1,
                stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free),
                stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free));
        outputChunkers_processChunkSequencePhased(outputChunkers, threadIdx, chunkIdx, bamChunk->refSeqName,
                                                  NULL, NULL, reads, readsBelongingToHap1, readsBelongingToHap2, gF,
                                                  params);

        // Cleanup
        if (chunkVcfEntries != NULL) stList_destruct(chunkVcfEntries);
        stSet_destruct(readsBelongingToHap1);
        stSet_destruct(readsBelongingToHap2);
        bubbleGraph_destruct(bg);
        stGenomeFragment_destruct(gF);
        stList_destruct(vcfEntriesToBubbles);
        free(chunkReference);

        // report timing
        if (st_getLogLevel() >= info) {
            st_logInfo(">%s Chunk with ~%"PRId64" reads processed in %d sec\n",
                       logIdentifier, stList_length(reads), (int) (time(NULL) - chunkStartTime));
        }

        // final post-completion logging cleanup
        stList_destruct(reads);
        free(logIdentifier);
    }

    // for writing haplotyped chunks
    stList *allReadIdsHap1 = stList_construct3(0, free);
    stList *allReadIdsHap2 = stList_construct3(0, free);

    // merge chunks
    time_t mergeStartTime = time(NULL);
    st_logCritical("> Starting merge\n");
    outputChunkers_stitchAndTrackExtraData(outputChunkers, TRUE, bamChunker->chunkCount, allReadIdsHap1, allReadIdsHap2,
            NULL);
    time_t mergeEndTime = time(NULL);
    char *tds = getTimeDescriptorFromSeconds((int) mergeEndTime - mergeStartTime);
    st_logCritical("> Merging took %s\n", tds);
    outputChunkers_destruct(outputChunkers);
    free(tds);
    tds = getTimeDescriptorFromSeconds((int) time(NULL) - mergeEndTime);
    st_logCritical("> Merge cleanup took %s\n", tds);
    free(tds);

    // maybe write final haplotyped bams
    // logging
    time_t hapBamStart = time(NULL);
    st_logInfo("> Writing final haplotyped BAMs\n");

    // get all reads
    stSet *allReadIdsForHaplotypingHap1 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    stSet *allReadIdsForHaplotypingHap2 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    for (int64_t i = 0; i < stList_length(allReadIdsHap1); i++) {
        stSet_insert(allReadIdsForHaplotypingHap1, stList_get(allReadIdsHap1, i));
    }
    for (int64_t i = 0; i < stList_length(allReadIdsHap2); i++) {
        stSet_insert(allReadIdsForHaplotypingHap2, stList_get(allReadIdsHap2, i));
    }

    // write it
    writeHaplotaggedBam(bamChunker->bamFile, outputBase, regionStr,
                        allReadIdsForHaplotypingHap1, allReadIdsForHaplotypingHap2, NULL, params, "");

    // loggit
    char *hapBamTDS = getTimeDescriptorFromSeconds(time(NULL) - hapBamStart);
    st_logCritical("> Wrote haplotyped bams in %s\n", hapBamTDS);

    // cleanup
    free(hapBamTDS);
    stSet_destruct(allReadIdsForHaplotypingHap1);
    stSet_destruct(allReadIdsForHaplotypingHap2);

    // cleanup
    bamChunker_destruct(bamChunker);
    params_destruct(params);
    if (regionStr != NULL) free(regionStr);
    stList_destruct(chunkOrder);
    free(vcfFile);
    stHash_destruct(vcfEntries);
    if (allReadIdsHap1 != NULL) stList_destruct(allReadIdsHap1);
    if (allReadIdsHap2 != NULL) stList_destruct(allReadIdsHap2);
    free(outputBase);
    free(bamInFile);
    free(referenceFastaFile);
    free(paramsFile);

    // log completion
    char *timeDescriptor = getTimeDescriptorFromSeconds(time(NULL) - startTime);
    st_logCritical("> Finished phasing in %s.\n", timeDescriptor);
    free(timeDescriptor);

//    while(1); // Use this for testing for memory leaks

    return 0;
}

