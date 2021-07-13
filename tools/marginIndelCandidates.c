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
#define TOOL_NAME "indelCandidates"

void indel_candidates_usage() {
    fprintf(stderr, "usage: %s <HAPLOTAGGED_BAM> <REFERENCE_FASTA> <INDEL_VARIANT_VCF> <PARAMS> [options]\n", TOOL_NAME);
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Tags reads in ALIGN_BAM using variants in VARIANT_VCF.\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    HAPLOTAGGED_BAM is the alignment of reads to the reference.\n");
    fprintf(stderr, "    REFERENCE_FASTA is the reference sequence BAM file in fasta format.\n");
    fprintf(stderr, "    INDEL_VARIANT_VCF is the set of variants to select INDEL candidates for.\n");
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

    fprintf(stderr, "\n");
}

int stIntTuple_compareIndexZero(const void *a, const void *b) {
    stIntTuple *A = (stIntTuple*) a;
    stIntTuple *B = (stIntTuple*) b;
    if (stIntTuple_get(A, 0) == stIntTuple_get(B, 0)) return 0;
    return stIntTuple_get(A, 0) < stIntTuple_get(B, 0) ? -1 : 1;
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

    if (argc < 4) {
        free(outputBase);
        free(logLevelString);
        indel_candidates_usage();
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
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "ha:o:t:r:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
        case 'a':
            free(logLevelString);
            logLevelString = stString_copy(optarg);
            break;
        case 'h':
            indel_candidates_usage();
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
        default:
            indel_candidates_usage();
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

    // Sanity check for parameters
    if (params->phaseParams->onlyUseSNPVCFEntries) {
        st_errAbort("> %s must have parameter phase.onlyUseSNPVCFEntries = \"false\"", TOOL_NAME);
    } else if (params->phaseParams->onlyUsePassVCFEntries) {
        st_errAbort("> %s must have parameter phase.onlyUsePassVCFEntries = \"false\"", TOOL_NAME);
    } else if (!params->phaseParams->includeHomozygousVCFEntries) {
        st_errAbort("> %s must have parameter phase.includeHomozygousVCFEntries = \"true\"", TOOL_NAME);
    } else if (params->phaseParams->useVariantSelectionAdaptiveSampling) {
        st_errAbort("> %s must have parameter phase.useVariantSelectionAdaptiveSampling = \"false\"", TOOL_NAME);
    } else if (params->phaseParams->minVariantQuality > 0) {
        st_errAbort("> %s must have parameter phase.minVariantQuality = 0", TOOL_NAME);
    } else if (!params->polishParams->useRunLengthEncoding) {
        st_errAbort("> %s must have parameter polish.useRunLengthEncoding = \"true\"", TOOL_NAME);
    } else if (!params->polishParams->useRepeatCountsInAlignment) {
        st_errAbort("> %s must have parameter polish.useRepeatCountsInAlignment = \"true\"", TOOL_NAME);
    } else if (params->polishParams->repeatSubMatrix == NULL) {
        st_errAbort("> %s must have parameter polish.repeatCountSubstitutionMatrix set", TOOL_NAME);
    }

    // used during realignment
    uint64_t maximumRepeatLength = (uint64_t) params->polishParams->repeatSubMatrix->maximumRepeatLength;

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

        // get all VCF entries
        stList *allChunkVcfEntries = getVcfEntriesForRegion(vcfEntries, NULL, bamChunk->refSeqName,
                                                 bamChunk->chunkOverlapStart,  bamChunk->chunkOverlapEnd, params);

        // filter for only INDEL
        int64_t snpVcfEntryCount = 0;
        stList *chunkVcfEntries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);
        for (int64_t v = 0; v < stList_length(allChunkVcfEntries); v++) {
            VcfEntry *vcfEntry = stList_get(allChunkVcfEntries, v);
            RleString *refAllele = stList_get(vcfEntry->alleles, 0);
            bool isIndel = FALSE;
            for (int64_t a = 1; a < stList_length(vcfEntry->alleles); a++) {
                // is this an indel?
                if (refAllele->nonRleLength != ((RleString*) stList_get(vcfEntry->alleles, a))->nonRleLength) {
                    isIndel = TRUE;
                    break;
                }
            }
            if (isIndel) {
                stList_append(chunkVcfEntries, vcfEntry);
            } else {
                snpVcfEntryCount++;
                vcfEntry_destruct(vcfEntry);
            }
        }
        st_logInfo(" %s Removed %"PRId64" SNP and kept %"PRId64" INDEL variants\n", logIdentifier, snpVcfEntryCount, stList_length(chunkVcfEntries));
        stList_setDestructor(allChunkVcfEntries, NULL);
        stList_destruct(allChunkVcfEntries);

        updateVcfEntriesWithSubstringsAndPositionsByRleLength(chunkVcfEntries, chunkReference, strlen(chunkReference),
                FALSE, params);

        // Convert bam lines into corresponding reads and alignments
        st_logInfo(" %s Parsing input reads from file: %s\n", logIdentifier, bamInFile);
        stList *reads = stList_construct3(0, (void (*)(void *)) bamChunkRead_destruct);
        extractReadSubstringsAtVariantPositions(bamChunk, chunkVcfEntries, reads, NULL, params->polishParams);

        // prep for work
        stHash *vcfEntriesToReadSubstrings = buildVcfEntryToReadSubstringsMap(reads, params);
        int64_t vcfEntriesWithoutSubstrings = 0;

        // look at all candidate alleles around vcf entries and pick best alleles
        for (int64_t v = 0; v < stList_length(chunkVcfEntries); v++) {
            VcfEntry *vcfEntry = stList_get(chunkVcfEntries, v);

            // the infrastructure gets these and we want that, but we don't do anything with them
            if (vcfEntry->rawRefPosInformativeOnly < bamChunk->chunkStart || vcfEntry->rawRefPosInformativeOnly >= bamChunk->chunkEnd) {
                continue;
            }

            stList *alleles = vcfEntry->alleleSubstrings;
            assert(stList_length(alleles) >= 2);

            // Get read substrings
            stList *readSubstrings = stHash_search(vcfEntriesToReadSubstrings, vcfEntry);

            // nothing to phase with
            if (stList_length(readSubstrings) == 0) {
                stList_destruct(readSubstrings);
                vcfEntriesWithoutSubstrings++;
                continue;
            }

            // get initial graph
            RleString *reference = stList_get(alleles, 0);
            Poa *poaH1 = poa_getReferenceGraph(reference, params->polishParams->alphabet, maximumRepeatLength);
            Poa *poaH2 = poa_getReferenceGraph(reference, params->polishParams->alphabet, maximumRepeatLength);
            stList *readSubstringsAsBCRs = stList_construct3(0, (void(*)(void*))bamChunkRead_destruct);

            // For each read, update the POAs
            for (int64_t r = 0; r < stList_length(readSubstrings); r++) {
                BamChunkReadSubstring *bcrs = stList_get(readSubstrings, r);

                // Generate set of posterior probabilities for matches, deletes and inserts with respect to reference.
                stList *matches = NULL, *inserts = NULL, *deletes = NULL;

                SymbolString sX = rleString_constructSymbolString(reference, 0, reference->length,
                                                                  params->polishParams->alphabet,
                                                                  params->polishParams->useRepeatCountsInAlignment,
                                                                  poaH1->maxRepeatCount - 1);
                SymbolString sY = rleString_constructSymbolString(bcrs->substring, 0, bcrs->substring->length,
                                                                  params->polishParams->alphabet,
                                                                  params->polishParams->useRepeatCountsInAlignment,
                                                                  poaH1->maxRepeatCount - 1);

                getAlignedPairsWithIndels(bcrs->read->forwardStrand ? params->polishParams->stateMachineForForwardStrandRead :
                                          params->polishParams->stateMachineForReverseStrandRead,
                                          sX, sY, params->polishParams->p, &matches, &deletes, &inserts, FALSE, FALSE);

                char *expandedRead = rleString_expand(bcrs->substring);
                BamChunkRead *readSubstringAsBCR = bamChunkRead_construct3(stString_copy(bcrs->read->readName),
                        expandedRead, NULL, bcrs->read->forwardStrand, 0, TRUE);
                stList_append(readSubstringsAsBCRs, readSubstringAsBCR);

                // Add weights, edges and nodes to the poa
                if (bcrs->read->haplotype == 1) {
                    poa_augment(poaH1, bcrs->substring, bcrs->read->forwardStrand, r, matches, inserts, deletes,
                                params->polishParams);
                } else if (bcrs->read->haplotype == 2) {
                    poa_augment(poaH2, bcrs->substring, bcrs->read->forwardStrand, r, matches, inserts, deletes,
                                params->polishParams);
                } else {
                    poa_augment(poaH1, bcrs->substring, bcrs->read->forwardStrand, r, matches, inserts, deletes,
                                params->polishParams);
                    poa_augment(poaH2, bcrs->substring, bcrs->read->forwardStrand, r, matches, inserts, deletes,
                                params->polishParams);
                }

                // cleanup
                stList_destruct(matches);
                stList_destruct(inserts);
                stList_destruct(deletes);
                free(expandedRead);
                symbolString_destruct(sX);
                symbolString_destruct(sY);
            }

            // sanity
            assert(stList_length(readSubstrings) == stList_length(readSubstringsAsBCRs));
            for (int64_t r = 0; r < stList_length(readSubstrings); r++) {
                BamChunkRead *rs = ((BamChunkReadSubstring*)stList_get(readSubstrings,r))->read;
                BamChunkRead *rsabcr = stList_get(readSubstringsAsBCRs, r);
                assert(stString_eq(rs->readName, rsabcr->readName));
            }

            // debug logging
            char *prevRefH1 = "";
            char *prevRefH2 = "";
            char *postRefH1 = "";
            char *postRefH2 = "";
            if (st_getLogLevel() == info) {
                prevRefH1 = rleString_expand(poaH1->refString);
                prevRefH2 = rleString_expand(poaH2->refString);
            }

            // estimate from observations
            poa_estimateRepeatCountsUsingBayesianModel(poaH1, readSubstringsAsBCRs, params->polishParams->repeatSubMatrix);
            poa_estimateRepeatCountsUsingBayesianModel(poaH2, readSubstringsAsBCRs, params->polishParams->repeatSubMatrix);

            // finish logging
            if (st_getLogLevel() == info) {
                postRefH1 = rleString_expand(poaH1->refString);
                postRefH2 = rleString_expand(poaH2->refString);

                int64_t maxLen = strlen(prevRefH1);
                maxLen = maxLen > strlen(prevRefH2) ? maxLen : strlen(prevRefH2);
                maxLen = maxLen > strlen(postRefH1) ? maxLen : strlen(postRefH1);
                maxLen = maxLen > strlen(postRefH2) ? maxLen : strlen(postRefH2);

                char *logString = stString_print(" %%s   %%8s H1: %%%ds  H2: %%%ds \n", maxLen, maxLen);
                st_logInfo(" %s Estimating run lengths for variant at %s:%"PRId64"\n", logIdentifier,
                        vcfEntry->refSeqName, vcfEntry->refPos);
                st_logInfo(logString, logIdentifier, "Prev:", prevRefH1, prevRefH2);
                st_logInfo(logString, logIdentifier, "Post:", postRefH1, postRefH2);
                free(logString);
                free(prevRefH1);
                free(prevRefH2);
                free(postRefH1);
                free(postRefH2);
            }


            // get best poa reference strings per hap
            char *poaRefStringExpandedH1 = rleString_expand(poaH1->refString);
            RleString *poaRefStringExpandedRawLEH1 = rleString_construct_no_rle(poaRefStringExpandedH1);
            SymbolString sH1 = rleString_constructSymbolString(poaRefStringExpandedRawLEH1, 0, poaRefStringExpandedRawLEH1->length,
                                                               params->polishParams->alphabet,
                                                               params->polishParams->useRepeatCountsInAlignment,
                                                               poaH1->maxRepeatCount - 1);
            char *poaRefStringExpandedH2 = rleString_expand(poaH2->refString);
            RleString *poaRefStringExpandedRawLEH2 = rleString_construct_no_rle(poaRefStringExpandedH2);
            SymbolString sH2 = rleString_constructSymbolString(poaRefStringExpandedRawLEH2, 0, poaRefStringExpandedRawLEH2->length,
                                                               params->polishParams->alphabet,
                                                               params->polishParams->useRepeatCountsInAlignment,
                                                               poaH2->maxRepeatCount - 1);
            stList *alleleRankH1 = stList_construct3(0, (void(*)(void*)) stIntTuple_destruct);
            stList *alleleRankH2 = stList_construct3(0, (void(*)(void*)) stIntTuple_destruct);

            // rank alleles
            for (int64_t a = 0; a < stList_length(vcfEntry->alleles); a++) {

                // prep for alignment ranking
                stList *matchesH1 = NULL, *insertsH1 = NULL, *deletesH1 = NULL;
                stList *matchesH2 = NULL, *insertsH2 = NULL, *deletesH2 = NULL;
                double scoreH1, scoreH2;
                assert(vcfEntry->referenceSuffix != NULL);
                assert(vcfEntry->referencePrefix != NULL);

                // Construct allele string for alignment
                RleString *allele = stList_get(vcfEntry->alleles, a);
                char *alleleRaw = rleString_expand(allele);
                char *alleleSubstringRaw = stString_print("%s%s%s", vcfEntry->referencePrefix, alleleRaw,
                        vcfEntry->referenceSuffix);
                RleString *fullAlleleSubstringRawLE = rleString_construct_no_rle(alleleSubstringRaw);

                SymbolString sAllele = rleString_constructSymbolString(fullAlleleSubstringRawLE, 0, fullAlleleSubstringRawLE->length,
                        params->polishParams->alphabet, params->polishParams->useRepeatCountsInAlignment,
                        poaH1->maxRepeatCount - 1);

                // align and score
                getAlignedPairsWithIndels(params->polishParams->stateMachineForGenomeComparison, sH1, sAllele,
                        params->polishParams->p, &matchesH1, &deletesH1, &insertsH1, FALSE, FALSE);
                getAlignedPairsWithIndels(params->polishParams->stateMachineForGenomeComparison, sH2, sAllele,
                        params->polishParams->p, &matchesH2, &deletesH2, &insertsH2, FALSE, FALSE);
                stList *tmpH1 = getMaximalExpectedAccuracyPairwiseAlignment(matchesH1, deletesH1, insertsH1, sH1.length,
                                                                            sAllele.length, &scoreH1,
                                                                            params->polishParams->p);
                stList *tmpH2 = getMaximalExpectedAccuracyPairwiseAlignment(matchesH2, deletesH2, insertsH2, sH2.length,
                                                                            sAllele.length, &scoreH2,
                                                                            params->polishParams->p);

                // save
                stList_append(alleleRankH1, stIntTuple_construct2((int64_t) (10*scoreH1), a));
                stList_append(alleleRankH2, stIntTuple_construct2((int64_t) (10*scoreH2), a));

                // cleanup
                free(alleleRaw);
                free(alleleSubstringRaw);
                rleString_destruct(fullAlleleSubstringRawLE);
                stList_destruct(matchesH1);
                stList_destruct(matchesH2);
                stList_destruct(insertsH1);
                stList_destruct(insertsH2);
                stList_destruct(deletesH1);
                stList_destruct(deletesH2);
                stList_destruct(tmpH1);
                stList_destruct(tmpH2);
                symbolString_destruct(sAllele);
            }

            // print ranked alleles
            stList_sort(alleleRankH1, stIntTuple_compareIndexZero);
            stList_sort(alleleRankH2, stIntTuple_compareIndexZero);

            // loggit (for now)
            st_logInfo(" %s Hap1 Allele Rank:\n", logIdentifier);
            st_logInfo(" %s   H01: %s\n", logIdentifier, rleString_expand(poaH1->refString));
            for (int64_t a = stList_length(alleleRankH1) - 1; a >= 0; a--) {
                stIntTuple *it = stList_get(alleleRankH1, a);
                int64_t score = stIntTuple_get(it, 0);
                int64_t alleleIdx = stIntTuple_get(it, 1);
                RleString *allele = stList_get(vcfEntry->alleles, alleleIdx);
                char *alleleRaw = rleString_expand(allele);
                char *alleleSubstringRaw = stString_print("%s%s%s", vcfEntry->referencePrefix, alleleRaw,
                                                          vcfEntry->referenceSuffix);
                st_logInfo(" %s   A%2"PRId64": %s (%"PRId64")\n", logIdentifier, alleleIdx, alleleSubstringRaw, score);
                free(alleleRaw);
                free(alleleSubstringRaw);
            }
            st_logInfo(" %s Hap2 Allele Rank:\n", logIdentifier);
            st_logInfo(" %s   H02: %s\n", logIdentifier, rleString_expand(poaH2->refString));
            for (int64_t a = stList_length(alleleRankH2) - 1; a >= 0; a--) {
                stIntTuple *it = stList_get(alleleRankH2, a);
                int64_t score = stIntTuple_get(it, 0);
                int64_t alleleIdx = stIntTuple_get(it, 1);
                RleString *allele = stList_get(vcfEntry->alleles, alleleIdx);
                char *alleleRaw = rleString_expand(allele);
                char *alleleSubstringRaw = stString_print("%s%s%s", vcfEntry->referencePrefix, alleleRaw,
                                                          vcfEntry->referenceSuffix);
                st_logInfo(" %s   A%2"PRId64": %s (%"PRId64")\n", logIdentifier, alleleIdx, alleleSubstringRaw, score);
                free(alleleRaw);
                free(alleleSubstringRaw);
            }

            // update the alleles
            stIntTuple *hap1AlleleTuple = stList_pop(alleleRankH1);
            stIntTuple *hap2AlleleTuple = stList_pop(alleleRankH2);
            vcfEntry->rootVcfEntry->gt1 = stIntTuple_get(hap1AlleleTuple, 1);
            vcfEntry->rootVcfEntry->gt2 = stIntTuple_get(hap2AlleleTuple, 1);
            vcfEntry->rootVcfEntry->wasAnaylzed = TRUE;

            // cleanup
            symbolString_destruct(sH1);
            symbolString_destruct(sH2);
            stList_destruct(readSubstringsAsBCRs);
            poa_destruct(poaH1);
            poa_destruct(poaH2);
        }

        // Cleanup
        if (chunkVcfEntries != NULL) stList_destruct(chunkVcfEntries);
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


    // loggit
    time_t vcfWriteStart = time(NULL);
    char *outputVcfFile = stString_print("%s.candidates.vcf", outputBase);
    st_logCritical("> Writing VCF to %s\n", outputVcfFile);

    // write it
    writeCandidateVcf(vcfFile, regionStr, outputVcfFile, vcfEntries, params);

    // loggit
    char *phasedVcfTDS = getTimeDescriptorFromSeconds(time(NULL) - vcfWriteStart);
    st_logCritical("> Wrote candidate VCF in %s\n", phasedVcfTDS);

    // cleanup
    free(phasedVcfTDS);
    free(outputVcfFile);

    // cleanup
    bamChunker_destruct(bamChunker);
    params_destruct(params);
    if (regionStr != NULL) free(regionStr);
    stList_destruct(chunkOrder);
    free(vcfFile);
    stHash_destruct(vcfEntries);
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

