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
    fprintf(stderr, "usage: marginPolish <BAM_FILE> <ASSEMBLY_FASTA> <PARAMS> [options]\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Polishes the ASSEMBLY_FASTA using alignments in BAM_FILE.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    BAM_FILE is the alignment of reads to the assembly (or reference).\n");
    fprintf(stderr, "    ASSEMBLY_FASTA is the reference sequence BAM file in fasta format.\n");
    fprintf(stderr, "    PARAMS is the file with marginPolish parameters.\n");

    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help                : Print this help screen\n");
    fprintf(stderr, "    -a --logLevel            : Set the log level [default = info]\n");
# ifdef _OPENMP
    fprintf(stderr, "    -t --threads             : Set number of concurrent threads [default = 1]\n");
#endif
    fprintf(stderr, "    -o --outputBase          : Name to use for output files [default = 'output']\n");
    fprintf(stderr, "    -r --region              : If set, will only compute for given chromosomal region.\n");
    fprintf(stderr, "                                 Format: chr:start_pos-end_pos (chr3:2000-3000).\n");
    fprintf(stderr, "    -p --depth               : Will override the downsampling depth set in PARAMS.\n");
    fprintf(stderr, "    -2 --diploid             : Will perform diploid phasing.\n");
    fprintf(stderr, "    -P --minPartitionPhred   : Min Phred-scale liklihood for partition inclusion (diploid)\n");
    fprintf(stderr, "    -k --tempFilesToDisk     : Write temporary files to disk (for --diploid or supplementary output).\n");

# ifdef _HDF5
    fprintf(stderr, "\nHELEN feature generation options:\n");
    fprintf(stderr, "    -f --produceFeatures     : output splitRleWeight or diploidRleWeight (based on -2 flag) features for HELEN.\n");
    fprintf(stderr, "    -F --featureType         : output specific feature type for HELEN (overwrites -f).  Valid types:\n");
    fprintf(stderr, "                                 splitRleWeight:   [default] run lengths split into chunks\n");
    fprintf(stderr, "                                 channelRleWeight: run lengths split into per-nucleotide channels\n");
    fprintf(stderr, "                                 simpleWeight:     weighted likelihood from POA nodes (non-RLE)\n");
    fprintf(stderr, "                                 diploidRleWeight: [default] produces diploid features \n");
    fprintf(stderr, "    -L --splitRleWeightMaxRL : max run length (for RLE feature types) \n");
    fprintf(stderr, "                                 [split default = %d, channel default = %d, diploid default = %d]\n",
            POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT, POAFEATURE_CHANNEL_MAX_RUN_LENGTH_DEFAULT, POAFEATURE_DIPLOID_MAX_RUN_LENGTH_DEFAULT);
    fprintf(stderr, "    -u --trueReferenceBam    : true reference aligned to ASSEMBLY_FASTA, for HELEN\n");
    fprintf(stderr, "                               features.  Setting this parameter will include labels\n");
    fprintf(stderr, "                               in output.  If -2/--diploid is set, this parameter must\n");
    fprintf(stderr, "                               contain two comma-separated values\n");
    # endif

    fprintf(stderr, "\nMiscellaneous supplementary output options:\n");
    fprintf(stderr, "    -c --supplementaryChunks : Write supplementary files for each chunk (in additon to writing\n");
    fprintf(stderr, "                               whole genome information)\n");
    fprintf(stderr, "    -C --supplementaryChunksOnly : Only write supplementary files for each chunk (will not write\n");
    fprintf(stderr, "                               whole genome information)\n");
    fprintf(stderr, "    -d --outputPoaDot        : Write out the poa as DOT file (only done per chunk)\n");
    fprintf(stderr, "    -i --outputRepeatCounts  : Write out the repeat counts as CSV file\n");
    fprintf(stderr, "    -j --outputPoaCsv        : Write out the poa as CSV file\n");
    fprintf(stderr, "    -n --outputHaplotypeReads: Write out phased reads and likelihoods as CSV file\n");
    fprintf(stderr, "    -m --outputHaplotypeBAM  : Write out phased BAMs\n");
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
    int numThreads = 1;
    int64_t maxDepth = -1;
    int64_t minPhredScoreForHaplotypePartition = -1;
    bool diploid = FALSE;
    bool inMemory = TRUE;

    // for feature generation
    HelenFeatureType helenFeatureType = HFEAT_NONE;
    bool setDefaultHelenFeature = false;
    char *trueReferenceBam = NULL;
    char *trueReferenceBamHap2 = NULL;
    bool fullFeatureOutput = FALSE;
    int64_t splitWeightMaxRunLength = 0;
    void **helenHDF5Files = NULL;

    // for supplementary output
    bool outputPoaDOT = FALSE;
    bool outputPoaCSV = FALSE;
    bool outputRepeatCounts = FALSE;
    bool outputHaplotypeReads = FALSE;
    bool outputHaplotypeBAM = FALSE;
    bool writeChunkSupplementaryOutput = FALSE;
    bool writeChunkSupplementaryOutputOnly = FALSE;

    if (argc < 4) {
        free(outputBase);
        free(logLevelString);
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    referenceFastaFile = stString_copy(argv[2]);
    paramsFile = stString_copy(argv[3]);

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
                { "depth", required_argument, 0, 'p'},
                { "minPartitionPhred", required_argument, 0, 'P'},
                { "diploid", no_argument, 0, '2'},
                { "produceFeatures", no_argument, 0, 'f'},
                { "featureType", required_argument, 0, 'F'},
                { "trueReferenceBam", required_argument, 0, 'u'},
                { "splitRleWeightMaxRL", required_argument, 0, 'L'},
				{ "supplementaryChunks", no_argument, 0, 'c'},
				{ "supplementaryChunksOnly", no_argument, 0, 'C'},
				{ "outputRepeatCounts", no_argument, 0, 'i'},
				{ "outputPoaCsv", no_argument, 0, 'j'},
				{ "outputPoaDot", no_argument, 0, 'd'},
				{ "outputHaplotypeBAM", no_argument, 0, 'm'},
				{ "outputHaplotypeReads", no_argument, 0, 'n'},
                { "tempFilesToDisk", no_argument, 0, 'k'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "ha:o:v:p:P:2t:r:fF:u:L:cCijdmnk", long_options, &option_index);

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
        case 'p':
            maxDepth = atoi(optarg);
            if (maxDepth < 0) {
                st_errAbort("Invalid maxDepth: %s", optarg);
            }
            break;
        case 'P':
            minPhredScoreForHaplotypePartition = atoi(optarg);
            break;
        case 'F':
            if (stString_eqcase(optarg, "simpleWeight") || stString_eqcase(optarg, "simple")) {
                helenFeatureType = HFEAT_SIMPLE_WEIGHT;
            } else if (stString_eqcase(optarg, "rleWeight") || stString_eqcase(optarg, "splitRleWeight") || stString_eqcase(optarg, "split")) {
                helenFeatureType = HFEAT_SPLIT_RLE_WEIGHT;
            } else if (stString_eqcase(optarg, "channelRleWeight") || stString_eqcase(optarg, "channel")) {
                helenFeatureType = HFEAT_CHANNEL_RLE_WEIGHT;
            } else if (stString_eqcase(optarg, "diploidRleWeight") || stString_eqcase(optarg, "diploid")) {
                helenFeatureType = HFEAT_DIPLOID_RLE_WEIGHT;
            } else {
                fprintf(stderr, "Unrecognized featureType for HELEN: %s\n\n", optarg);
                usage();
                return 1;
            }
            break;
        case 'u':
            trueReferenceBam = stString_copy(optarg);
            break;
        case 'f':
            if (helenFeatureType == HFEAT_NONE) {
                setDefaultHelenFeature = true;
            }
            break;
        case 'L':
            splitWeightMaxRunLength = atoi(optarg);
            if (splitWeightMaxRunLength <= 0) {
                st_errAbort("Invalid splitRleWeightMaxRL: %d", splitWeightMaxRunLength);
            }
            break;
        case 't':
            numThreads = atoi(optarg);
            if (numThreads <= 0) {
                st_errAbort("Invalid thread count: %d", numThreads);
            }
            break;
        case '2':
            diploid = TRUE;
            break;
        case 'k':
            inMemory = FALSE;
            break;
        case 'c':
            writeChunkSupplementaryOutput = TRUE;
            break;
        case 'C':
            writeChunkSupplementaryOutputOnly = TRUE;
            break;
        case 'i':
            outputRepeatCounts = TRUE;
            break;
        case 'j':
            outputPoaCSV = TRUE;
            break;
        case 'd':
            outputPoaDOT = TRUE;
            break;
        case 'm':
            outputHaplotypeBAM = TRUE;
            break;
        case 'n':
            outputHaplotypeReads = TRUE;
            break;
        default:
            usage();
            free(outputBase);
            free(logLevelString);
            free(bamInFile);
            free(referenceFastaFile);
            free(paramsFile);
            if (trueReferenceBam != NULL) free(trueReferenceBam);
            return 0;
        }
    }

    // handle diploid sanity checks
    if (diploid) {
        // true ref
        if (trueReferenceBam != NULL) {
            stList *trueRefParts = stString_splitByString(trueReferenceBam, ",");
            if (stList_length(trueRefParts) != 2) {
                st_errAbort("If --diploid is set, --trueReferenceBam must have two comma-separated values.");
            }
            free(trueReferenceBam);
            trueReferenceBam = stString_copy(stList_get(trueRefParts, 0));
            trueReferenceBamHap2 = stString_copy(stList_get(trueRefParts, 1));
            stList_destruct(trueRefParts);
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
    if (access(paramsFile, R_OK) != 0) {
        st_errAbort("Could not read from params file: %s\n", paramsFile);
    }
    if (trueReferenceBam != NULL) {
        if (access(trueReferenceBam, R_OK) != 0) {
            st_errAbort("Could not read from truth file: %s\n", trueReferenceBam);
        }
        char *idx = stString_print("%s.bai", trueReferenceBam);
        if (access(idx, R_OK) != 0) {
            st_errAbort("BAM does not appear to be indexed: %s\n", trueReferenceBam);
        }
        free(idx);
    }
    if (trueReferenceBamHap2 != NULL) {
        if (access(trueReferenceBamHap2, R_OK ) != 0 ) {
            st_errAbort("Could not read from truth h2 file: %s\n", trueReferenceBamHap2);
        }
        if (stString_eq(trueReferenceBam, trueReferenceBamHap2)) {
            st_errAbort("Truth H1 and H2 files are the same: %", trueReferenceBam);
        }
        char *idx = stString_print("%s.bai", trueReferenceBamHap2);
        if (access(idx, R_OK ) != 0 ) {
            st_errAbort("BAM does not appear to be indexed: %s\n", trueReferenceBamHap2);
        }
        free(idx);
    }

    // Initialization from arguments
    time_t startTime = time(NULL);
    st_setLogLevelFromString(logLevelString);
    free(logLevelString);
# ifdef _OPENMP
    if (numThreads <= 0) {
        numThreads = 1;
    }
    omp_set_num_threads(numThreads);
    st_logCritical("Running OpenMP with %d threads.\n", omp_get_max_threads());
    # endif

    // feature init
    if (helenFeatureType == HFEAT_NONE && setDefaultHelenFeature) {
        helenFeatureType = diploid ? HFEAT_DIPLOID_RLE_WEIGHT : HFEAT_SPLIT_RLE_WEIGHT;
    }
    if (helenFeatureType != HFEAT_NONE && splitWeightMaxRunLength == 0) {
        switch (helenFeatureType) {
            case HFEAT_SPLIT_RLE_WEIGHT:
                splitWeightMaxRunLength = POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT;
                break;
            case HFEAT_CHANNEL_RLE_WEIGHT:
                splitWeightMaxRunLength = POAFEATURE_CHANNEL_MAX_RUN_LENGTH_DEFAULT;
                break;
            case HFEAT_DIPLOID_RLE_WEIGHT:
                splitWeightMaxRunLength = POAFEATURE_DIPLOID_MAX_RUN_LENGTH_DEFAULT;
                break;
            default:
                break;
        }
    }

    // Parse parameters
    st_logCritical("> Parsing model parameters from file: %s\n", paramsFile);
    Params *params = params_readParams(paramsFile);

    // update depth (if set)
    if (maxDepth >= 0) {
        st_logCritical("> Changing maxDepth paramter from %"PRId64" to %"PRId64"\n", params->polishParams->maxDepth,
                       maxDepth);
        params->polishParams->maxDepth = (uint64_t) maxDepth;
    }
    if (minPhredScoreForHaplotypePartition >= 0) {
        st_logCritical("> Changing minPhredScoreForHaplotypePartition paramter from %"PRId64" to %"PRId64"\n",
                       params->phaseParams->minPhredScoreForHaplotypePartition, minPhredScoreForHaplotypePartition);
        params->phaseParams->minPhredScoreForHaplotypePartition = minPhredScoreForHaplotypePartition;
    }

    // Set no RLE if appropriate feature type is set
    if (helenFeatureType == HFEAT_SIMPLE_WEIGHT) {
        if (params->polishParams->useRunLengthEncoding) {
            st_errAbort("Invalid runLengthEncoding parameter because of HELEN feature type.\n");
        }
        // everthing else requires RLE
    } else if (helenFeatureType != HFEAT_NONE) {
        if (!params->polishParams->useRunLengthEncoding) {
            st_errAbort("Invalid runLengthEncoding parameter because of HELEN feature type.\n");
        }
    }

    // Print a report of the parsed parameters
    if (st_getLogLevel() == debug) {
        params_printParameters(params, stderr);
    }

    // get reference sequences (and remove cruft after refName)
    stHash *referenceSequences = parseReferenceSequences(referenceFastaFile);

    // get chunker for bam.  if regionStr is NULL, it will be ignored
    BamChunker *bamChunker = bamChunker_construct2(bamInFile, regionStr, params->polishParams);
    st_logCritical(
            "> Set up bam chunker with chunk size %i and overlap %i (for region=%s), resulting in %i total chunks\n",
            (int) bamChunker->chunkSize, (int) bamChunker->chunkBoundary, regionStr == NULL ? "all" : regionStr,
            bamChunker->chunkCount);
    if (bamChunker->chunkCount == 0) {
        st_errAbort("> Found no valid reads!\n");
    }

    // for feature generation
    #ifdef _HDF5
    if (helenFeatureType != HFEAT_NONE) {
        helenHDF5Files = (void **) openHelenFeatureHDF5FilesByThreadCount(outputBase, numThreads);
    }
    #endif

    // output info
    char *outputSequenceFile = stString_print("%s.fa", outputBase);
    char *outputReadCsvFile = stString_print("%s.reads.csv", outputBase);
    char *outputPoaCsvFile = stString_print("%s.poa.csv", outputBase);
    char *outputRepeatCountFile = stString_print("%s.repeatCount.csv", outputBase);

    // output chunker tracks intermediate output files
    OutputChunkers *outputChunkers = outputChunkers_construct(numThreads, params, outputSequenceFile,
            outputPoaCSV && !writeChunkSupplementaryOutputOnly ? outputPoaCsvFile : NULL,
            outputHaplotypeReads && !writeChunkSupplementaryOutputOnly ? outputReadCsvFile : NULL,
            outputRepeatCounts && !writeChunkSupplementaryOutputOnly ? outputRepeatCountFile : NULL,
            diploid ? ".hap1" : "", diploid ? ".hap2" : NULL, inMemory);

    // (may) need to shuffle chunks
    stList *chunkOrder = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < bamChunker->chunkCount; i++) {
        stList_append(chunkOrder, stIntTuple_construct1(i));
    }
    if (params->polishParams->shuffleChunks) {
        st_logInfo("> Ordering chunks by estimated depth.\n");
        stList_sort2(chunkOrder, compareBamChunkDepthByIndexInList, bamChunker->chunks);
        stList_reverse(chunkOrder);
    }

    // multiproccess the chunks, save to results
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
        char *fullReferenceString = stHash_search(referenceSequences, bamChunk->refSeqName);
        if (fullReferenceString == NULL) {
            st_errAbort(
                    "ERROR: Reference sequence missing from reference map: %s. Perhaps the BAM and REF are mismatched?",
                    bamChunk->refSeqName);
        }
        int64_t fullRefLen = strlen(fullReferenceString);
        if (bamChunk->chunkBoundaryStart > fullRefLen) {
            st_errAbort("ERROR: Reference sequence %s has length %"PRId64", chunk %"PRId64" has start position %"
                        PRId64". Perhaps the BAM and REF are mismatched?",
                        bamChunk->refSeqName, fullRefLen, chunkIdx, bamChunk->chunkBoundaryStart);
        }
        RleString *rleReference = bamChunk_getReferenceSubstring(bamChunk, referenceSequences, params);
        st_logInfo(">%s Going to process a chunk (~%"PRId64"x) for reference sequence: %s, starting at: %i and ending at: %i\n",
                   logIdentifier, bamChunk->estimatedDepth, bamChunk->refSeqName, (int) bamChunk->chunkBoundaryStart,
                   (int) (fullRefLen < bamChunk->chunkBoundaryEnd ? fullRefLen : bamChunk->chunkBoundaryEnd));

        // Convert bam lines into corresponding reads and alignments
        st_logInfo(" %s Parsing input reads from file: %s\n", logIdentifier, bamInFile);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        stList *alignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
        convertToReadsAndAlignments(bamChunk, rleReference, reads, alignments, params->polishParams);

        // do downsampling if appropriate
        if (params->polishParams->maxDepth > 0) {
            // get downsampling structures
            stList *filteredReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
            stList *discardedReads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
            stList *filteredAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);
            stList *discardedAlignments = stList_construct3(0, (void (*)(void *)) stList_destruct);

            bool didDownsample = poorMansDownsample(params->polishParams->maxDepth, bamChunk, reads, alignments,
                                                    filteredReads, filteredAlignments, discardedReads,
                                                    discardedAlignments);

            // we need to destroy the discarded reads and structures
            if (didDownsample) {
                st_logInfo(" %s Downsampled from %"PRId64" to %"PRId64" reads\n", logIdentifier,
                           stList_length(reads), stList_length(filteredReads));
                // free all reads and alignments not used
                stList_destruct(discardedReads);
                stList_destruct(discardedAlignments);
                // still has all the old reads, need to not free these
                stList_setDestructor(reads, NULL);
                stList_setDestructor(alignments, NULL);
                stList_destruct(reads);
                stList_destruct(alignments);
                // and keep the filtered reads
                reads = filteredReads;
                alignments = filteredAlignments;
            }
                // no downsampling, we just need to free the (empty) objects
            else {
                stList_destruct(filteredReads);
                stList_destruct(filteredAlignments);
                stList_destruct(discardedReads);
                stList_destruct(discardedAlignments);
            }
        }

        // prep for ploishing
        Poa *poa = NULL; // The poa alignment
        char *polishedConsensusString = NULL; // The polished reference string

        // Run the polishing method
        int64_t totalNucleotides = 0;
        if (st_getLogLevel() >= info) {
            for (int64_t u = 0; u < stList_length(reads); u++) {
                totalNucleotides += strlen(((BamChunkRead *) stList_get(reads, u))->rleRead->rleString);
            }
            st_logInfo(" %s Running polishing algorithm with %"PRId64" reads and %"PRIu64"K nucleotides\n",
                       logIdentifier, stList_length(reads), totalNucleotides >> 10);
        }

        // Generate partial order alignment (POA) (destroys rleAlignments in the process)
        poa = poa_realignAll(reads, alignments, rleReference, params->polishParams);

        // Log info about the POA
        if (st_getLogLevel() >= info) {
            st_logInfo(" %s Summary stats for POA:\t", logIdentifier);
            poa_printSummaryStats(poa, stderr);
        }
        if (st_getLogLevel() >= debug) {
            poa_print(poa, stderr, reads, 5);
        }

        // Write any optional outputs about repeat count and POA, etc.
        if (writeChunkSupplementaryOutput || writeChunkSupplementaryOutputOnly) {
            poa_writeSupplementalChunkInformation(outputBase, chunkIdx, bamChunk, poa, reads, params,
                                                  outputPoaDOT, outputPoaCSV, outputRepeatCounts);
        }

        // handle diploid case
        if(diploid) {
            // Get the bubble graph representation
            BubbleGraph *bg = bubbleGraph_constructFromPoa2(poa, reads, params->polishParams, TRUE);

            // Now make a POA for each of the haplotypes
            stHash *readsToPSeqs;
            stReference *ref = bubbleGraph_getReference(bg, bamChunk->refSeqName, params);
            stGenomeFragment *gf = bubbleGraph_phaseBubbleGraph(bg, ref, reads, params, &readsToPSeqs);

            stSet *readsBelongingToHap1, *readsBelongingToHap2;
            stGenomeFragment_phaseBamChunkReads(gf, readsToPSeqs, reads, &readsBelongingToHap1, &readsBelongingToHap2,
                    params->phaseParams);
            st_logInfo(" %s After phasing, of %i reads got %i reads partitioned into hap1 and %i reads partitioned "
                       "into hap2 (%i unphased)\n", logIdentifier, (int) stList_length(reads),
                       (int) stSet_size(readsBelongingToHap1), (int) stSet_size(readsBelongingToHap2),
                       (int) (stList_length(reads) - stSet_size(readsBelongingToHap1) -
                              stSet_size(readsBelongingToHap2)));

            // Debug report of hets
            if (st_getLogLevel() <= info) {
                uint64_t totalHets = 0;
                for (uint64_t h = 0; h < gf->length; h++) {
                    Bubble *b = &bg->bubbles[h + gf->refStart];
                    if (gf->haplotypeString1[h] != gf->haplotypeString2[h]) {
                        st_logDebug(" %s Got predicted het at bubble %i %s %s\n", logIdentifier, (int) h + gf->refStart,
                                    b->alleles[gf->haplotypeString1[h]]->rleString,
                                    b->alleles[gf->haplotypeString2[h]]->rleString);
                        totalHets++;
                    } else if (!rleString_eq(b->alleles[gf->haplotypeString1[h]], b->refAllele)) {
                        st_logDebug(" %s Got predicted hom alt at bubble %i %i\n", logIdentifier, (int) h + gf->refStart,
                                    (int) gf->haplotypeString1[h]);
                    }
                }
                st_logInfo(" %s In phasing chunk, got: %i hets from: %i total sites (fraction: %f)\n", logIdentifier,
                        (int) totalHets, (int) gf->length, (float) totalHets / gf->length);
            }

            st_logInfo(" %s Building POA for each haplotype\n", logIdentifier);
            uint64_t *hap1 = getPaddedHaplotypeString(gf->haplotypeString1, gf, bg, params);
            uint64_t *hap2 = getPaddedHaplotypeString(gf->haplotypeString2, gf, bg, params);

            Poa *poa_hap1 = bubbleGraph_getNewPoa(bg, hap1, poa, reads, params);
            Poa *poa_hap2 = bubbleGraph_getNewPoa(bg, hap2, poa, reads, params);

            if(params->polishParams->useRunLengthEncoding) {
                st_logInfo(" %s Using read phasing to reestimate repeat counts in phased manner\n", logIdentifier);
                poa_estimatePhasedRepeatCountsUsingBayesianModel(poa_hap1, reads, params->polishParams->repeatSubMatrix,
                                                                 readsBelongingToHap1, readsBelongingToHap2, params->polishParams);
                poa_estimatePhasedRepeatCountsUsingBayesianModel(poa_hap2, reads, params->polishParams->repeatSubMatrix,
                                                                 readsBelongingToHap2, readsBelongingToHap1, params->polishParams);
            }

            // Output
            outputChunkers_processChunkSequencePhased(outputChunkers, threadIdx, chunkIdx, bamChunk->refSeqName,
                                                      poa_hap1, poa_hap2, reads,
                                                      readsBelongingToHap1, readsBelongingToHap2, gf, params);
            RleString *polishedRleConsensusH1 = rleString_copy(poa_hap1->refString);
            RleString *polishedRleConsensusH2 = rleString_copy(poa_hap2->refString);
            char *polishedConsensusStringH1 = rleString_expand(polishedRleConsensusH1);
            char *polishedConsensusStringH2 = rleString_expand(polishedRleConsensusH2);

            //ancillary files
            if (writeChunkSupplementaryOutput || writeChunkSupplementaryOutputOnly) {
                poa_writeSupplementalChunkInformationDiploid(outputBase, chunkIdx, bamChunk, gf, poa_hap1, poa_hap2,
                        reads, readsBelongingToHap1, readsBelongingToHap2, params, outputPoaDOT, outputPoaCSV,
                        outputRepeatCounts, outputHaplotypeReads, outputHaplotypeBAM, logIdentifier);
            }

            // helen
            #ifdef _HDF5
            if (helenFeatureType != HFEAT_NONE) {
                PoaFeature_handleDiploidHelenFeatures(helenFeatureType,
                                                      splitWeightMaxRunLength,
                                                      helenHDF5Files, fullFeatureOutput, trueReferenceBam,
                                                      trueReferenceBamHap2, params,
                                                      logIdentifier, chunkIdx, bamChunk, reads, poa_hap1, poa_hap2,
                                                      readsBelongingToHap1,
                                                      readsBelongingToHap2, polishedRleConsensusH1,
                                                      polishedRleConsensusH2, rleReference);
            }
            #endif

            // Cleanup
            free(hap1);
            free(hap2);
            stSet_destruct(readsBelongingToHap1);
            stSet_destruct(readsBelongingToHap2);
            rleString_destruct(polishedRleConsensusH1);
            rleString_destruct(polishedRleConsensusH2);
            free(polishedConsensusStringH1);
            free(polishedConsensusStringH2);
            bubbleGraph_destruct(bg);
            stGenomeFragment_destruct(gf);
            stReference_destruct(ref);
            poa_destruct(poa_hap1);
            poa_destruct(poa_hap2);
            stHash_destruct(readsToPSeqs);

        } else {

            // get polished reference string and expand RLE (regardless of whether RLE was applied)
            if (params->polishParams->useRunLengthEncoding) {
                poa_estimateRepeatCountsUsingBayesianModel(poa, reads, params->polishParams->repeatSubMatrix);
            }

            // output
            outputChunkers_processChunkSequence(outputChunkers, threadIdx, chunkIdx, bamChunk->refSeqName, poa, reads);


            //ancillary files
            if (writeChunkSupplementaryOutput || writeChunkSupplementaryOutputOnly) {
                poa_writeSupplementalChunkInformation(outputBase, chunkIdx, bamChunk, poa, reads, params,
                        outputPoaDOT, outputPoaCSV, outputRepeatCounts);
            }

            // HELEN feature outputs
            #ifdef _HDF5
            RleString *polishedRleConsensus = rleString_copy(poa->refString);
            polishedConsensusString = rleString_expand(polishedRleConsensus);
            if (helenFeatureType != HFEAT_NONE) {
                PoaFeature_handleHelenFeatures(helenFeatureType, splitWeightMaxRunLength,
                                               helenHDF5Files, fullFeatureOutput, trueReferenceBam, rleReference, params,
                                               logIdentifier, chunkIdx,
                                               bamChunk, poa, reads, polishedConsensusString, polishedRleConsensus);

            }
            free(polishedConsensusString);
            rleString_destruct(polishedRleConsensus);
            #endif
        }

        // report timing
        if (st_getLogLevel() >= info) {
            st_logInfo(">%s Chunk with %"PRId64" reads and %"PRIu64"K nucleotides processed in %d sec\n",
                       logIdentifier, stList_length(reads), totalNucleotides >> 10,
                       (int) (time(NULL) - chunkStartTime));
        }

        // Cleanup
        rleString_destruct(rleReference);
        poa_destruct(poa);
        stList_destruct(reads);
        stList_destruct(alignments);
        free(logIdentifier);
    }

    // for writing haplotyped chunks
    stList *allReadIdsHap1 = NULL;
    stList *allReadIdsHap2 = NULL;
    if (outputHaplotypeBAM && !writeChunkSupplementaryOutputOnly) {
        // setup
        allReadIdsHap1 = stList_construct();
        allReadIdsHap2 = stList_construct();
    }

    // merge chunks
    time_t mergeStartTime = time(NULL);
    if (params->polishParams->shuffleChunks || inMemory) {
        st_logCritical("> Starting merge\n");
        outputChunkers_stitchAndTrackReadIds(outputChunkers, diploid, bamChunker->chunkCount,
                allReadIdsHap1, allReadIdsHap2);
    } else {
        st_logCritical("> Starting linear merge\n");
        outputChunkers_stitchLinear(outputChunkers, diploid);
    }
    time_t mergeEndTime = time(NULL);
    char *tds = getTimeDescriptorFromSeconds((int) mergeEndTime - mergeStartTime);
    st_logCritical("> Merging took %s\n", tds);
    outputChunkers_destruct(outputChunkers);
    free(tds);
    tds = getTimeDescriptorFromSeconds((int) time(NULL) - mergeEndTime);
    st_logCritical("> Merge cleanup took %s\n", tds);
    free(tds);

    // maybe write final haplotyped bams
    if (allReadIdsHap1 != NULL && allReadIdsHap2 != NULL) {
        time_t hapBamStart = time(NULL);
        st_logInfo("> Writing final haplotyped BAMs\n");

        stSet *allReadIdsForHaplotypingHap1 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
        stSet *allReadIdsForHaplotypingHap2 = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);

        for(int64_t i = 0; i < stList_length(allReadIdsHap1); i++) {
            stSet_insert(allReadIdsForHaplotypingHap1, stList_get(allReadIdsHap1, i));
        }
        for(int64_t i = 0; i < stList_length(allReadIdsHap2); i++) {
            stSet_insert(allReadIdsForHaplotypingHap2, stList_get(allReadIdsHap2, i));
        }

        // write it
        BamChunk *whbBamChunk = NULL;
        if (regionStr != NULL) {
            char regionContig[128] = "";
            int regionStart = 0;
            int regionEnd = 0;
            int scanRet = sscanf(regionStr, "%[^:]:%d-%d", regionContig, &regionStart, &regionEnd);
            whbBamChunk = bamChunk_construct2(regionContig, -1, regionStart, regionStart, regionEnd,
                    regionEnd, 0, bamChunker);
        }
        writeHaplotaggedBam(whbBamChunk, bamChunker->bamFile, outputBase,
                            allReadIdsForHaplotypingHap1, allReadIdsForHaplotypingHap2, params, "");
        if (whbBamChunk != NULL) {
            bamChunk_destruct(whbBamChunk);
        }

        char *hapBamTDS = getTimeDescriptorFromSeconds(time(NULL) - hapBamStart);
        st_logInfo("> Wrote haplotyped bams in %s\n", hapBamTDS);

        // cleanup
        stSet_destruct(allReadIdsForHaplotypingHap1);
        stSet_destruct(allReadIdsForHaplotypingHap2);
        stList_destruct(allReadIdsHap1);
        stList_destruct(allReadIdsHap2);
        free(hapBamTDS);
    }

    // Cleanup
    bamChunker_destruct(bamChunker);
    stHash_destruct(referenceSequences);
    params_destruct(params);
    if (trueReferenceBam != NULL) free(trueReferenceBam);
    if (trueReferenceBamHap2 != NULL) free(trueReferenceBamHap2);
    if (regionStr != NULL) free(regionStr);
#ifdef _HDF5
    if (helenHDF5Files != NULL) {
        for (int64_t i = 0; i < numThreads; i++) {
            HelenFeatureHDF5FileInfo_destruct((HelenFeatureHDF5FileInfo *) helenHDF5Files[i]);
        }
        free(helenHDF5Files);
    }
    #endif
    stList_destruct(chunkOrder);
    free(outputSequenceFile);
    if (outputPoaCsvFile != NULL) free(outputPoaCsvFile);
    if (outputReadCsvFile != NULL) free(outputReadCsvFile);
    if (outputRepeatCountFile != NULL) free(outputRepeatCountFile);
    free(outputBase);
    free(bamInFile);
    free(referenceFastaFile);
    free(paramsFile);

    // log completion
    char *timeDescriptor = getTimeDescriptorFromSeconds(time(NULL) - startTime);
    st_logCritical("> Finished polishing in %s.\n", timeDescriptor);
    free(timeDescriptor);

//    while(1); // Use this for testing for memory leaks

    return 0;
}

