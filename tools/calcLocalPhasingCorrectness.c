/*
 * Copyright (C) 2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

//#include <stdio.h>
//#include <ctype.h>
//#include <memory.h>
//#include <hashTableC.h>
//#include <unistd.h>
//#include <time.h>
#include <getopt.h>
#include <math.h>

#include "marginVersion.h"
#include "localPhasingCorrectness.h"


void usage() {
    fprintf(stderr, "calcLocalPhasingCorrectness\n");
    fprintf(stderr, "Version: %s\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Generate LPC data for phase sets in both VCFs.\n\n");
    fprintf(stderr, "usage: calcLocalPhasingCorrectness [options] TRUTH_VCF QUERY_VCF > lpc_table.tsv \n\n");
    fprintf(stderr, "options:\n");
    fprintf(stderr, " -n, --grid-num INT     number of length scales to compute LPC for [50]\n");
    fprintf(stderr, " -s, --grid-skew FLOAT  controls evenness of grid between small and large values [0.0]\n");
    fprintf(stderr, " -q, --quiet            do not log progress to stderr\n");
    fprintf(stderr, " -h, --help             print this message and exit\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    
    // for logging
    st_setLogLevel(info);
    
    int64_t numLengthScales = 50;
    // < 1 give more high values, > 1 gives more low values
    double lowValueBias = 1.0;
    
    char* parseEnd = NULL;
    
    int c;
    while (true){
        static struct option long_options[] =
        {
            {"grid-num", required_argument, 0, 'n'},
            {"grid-skew", required_argument, 0, 's'},
            {"quiet", no_argument, 0, 'q'},
            {"help", no_argument, 0, 'h'},
            {0,0,0,0}
        };
        
        int option_index = 0;
        c = getopt_long (argc, argv, "n:s:qh?",
                         long_options, &option_index);
        if (c == -1){
            break;
        }
        
        switch(c){
            case 'n':
                numLengthScales = strtol(optarg, &parseEnd, 10);
                if ((parseEnd - optarg) != strlen(optarg)) {
                    fprintf(stderr, "error: Failed to parse argument %s as an integer\n\n", optarg);
                    usage();
                    exit(1);
                }
                break;
            case 's':
                lowValueBias = exp(-strtod(optarg, &parseEnd));
                if ((parseEnd - optarg) != strlen(optarg)) {
                    fprintf(stderr, "error: Failed to parse argument %s as a float\n\n", optarg);
                    usage();
                    exit(1);
                }
                break;
            case 'q':
                st_setLogLevel(critical);
                break;
            case 'h':
            case '?':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }
    
    if (argc - optind < 2) {
        fprintf(stderr, "error: Missing required positional arguments\n\n");
        usage();
        exit(1);
    }

    char *truthVcfFile = stString_copy(argv[optind++]);
    char *queryVcfFile = stString_copy(argv[optind++]);
    
    if (optind != argc) {
        fprintf(stderr, "error: Unused argument %s\n\n", argv[optind]);
        usage();
        exit(1);
    }
    
    if (numLengthScales < 2) {
        st_errAbort("error: Must have a grid of at least 2 values\n");
    }

    // sanity check (verify files are accessible)
    if (access(truthVcfFile, R_OK) != 0) {
        st_errAbort("error: Could not read from truth vcf file: %s\n", truthVcfFile);
    }
    if (access(queryVcfFile, R_OK) != 0) {
        st_errAbort("error: Could not read from query vcf file: %s\n", queryVcfFile);
    }
        
    st_logDebug("Decay values:\n");
    // use biased bin selection from ARGweaver paper, Rasmussen, et al. (2013)
    double denom = numLengthScales - 1;
    double *decayValues = (double*) malloc(numLengthScales * sizeof(double));
    for (int64_t i = 0; i < numLengthScales; ++i) {
        // handle the boundary conditions separately so we don't have to worry about
        // numerical imprecision
        if (i == 0) {
            decayValues[i] = 0.0;
        }
        else if (i == numLengthScales - 1) {
            decayValues[i] = 1.0;
        }
        else {
            decayValues[i] = (exp((i / denom) * log(1.0 + lowValueBias)) - 1.0) / lowValueBias;
        }
        st_logDebug("\t%f\n", decayValues[i]);
    }
    
    st_logDebug("Length scales:\n");
    double *lengthScales = (double*) malloc(numLengthScales * sizeof(double));
    for (int64_t i = 0; i < numLengthScales; ++i) {
        if (i == 0) {
            lengthScales[i] = 1.0;
        }
        else if (i == numLengthScales - 1) {
            lengthScales[i] = -log(0.0);
        }
        else {
            lengthScales[i] = 1.0 - log(2.0) / log(decayValues[i]);
        }
        st_logDebug("\t%f\n", lengthScales[i]);
    }
    

    // read and parse the VCFs
    st_logInfo("Reading VCF %s...\n", truthVcfFile);
    stHash *truthVariants = getPhasedVariants(truthVcfFile);
    st_logInfo("Reading VCF %s...\n", queryVcfFile);
    stHash *queryVariants = getPhasedVariants(queryVcfFile);
    
    free(truthVcfFile);
    free(queryVcfFile);

    stList *sharedContigs = getSharedContigs(truthVariants, queryVariants);
    st_logInfo("Found %"PRId64" shared contigs (truth %"PRId64", query %"PRId64")\n", stList_length(sharedContigs),
            stHash_size(truthVariants), stHash_size(queryVariants));
    
    
    double *correctnessValues = (double*) malloc(sizeof(double) * numLengthScales * stList_length(sharedContigs));
    double *meanCorrectnessValues = (double*) malloc(sizeof(double) * numLengthScales);
    
    int64_t reportIterations[5];
    for (int64_t i = 0; i < 5; ++i) {
        reportIterations[i] = ((i + 1) * numLengthScales) / 5;
    }
    int64_t nextReportIdx = 0;
    while (reportIterations[nextReportIdx] == 0 && nextReportIdx < 5) {
        ++nextReportIdx;
    }
    
    // Compute the correctness values
    for (int64_t i = 0; i < numLengthScales; ++i) {
        
        double weightedMeanNumer = 0.0;
        double weightedMeanDenom = 0.0;
        
        for (int64_t j = 0; j < stList_length(sharedContigs); ++j) {
            stList *contigTruthVariants = stHash_search(truthVariants, stList_get(sharedContigs, j));
            stList *contigQueryVariants = stHash_search(queryVariants, stList_get(sharedContigs, j));
            
            st_logDebug("\tComputing correctness for contig %s\n", stList_get(sharedContigs, j));
            
            int64_t phasedLength = 0;
            double correctness = phasingCorrectness(contigTruthVariants, contigQueryVariants, decayValues[i],
                                                    &phasedLength);
            
            correctnessValues[i * stList_length(sharedContigs) + j] = correctness;
            
            // TODO: would it be more consistent to weight by number of pairs?
            weightedMeanDenom += phasedLength;
            weightedMeanNumer += phasedLength * correctness;
        }
                
        meanCorrectnessValues[i] = weightedMeanNumer / weightedMeanDenom;
        
        if (i + 1 == reportIterations[nextReportIdx]) {
            st_logInfo("Finished computing correctness for %"PRId64" of %"PRId64" length scales\n", i + 1, numLengthScales);
            while (i + 1 == reportIterations[nextReportIdx] && nextReportIdx < 5) {
                ++nextReportIdx;
            }
        }
    }
    
    // print out the results in a table
    printf("decay\tlength_scale");
    for (int64_t i = 0; i < stList_length(sharedContigs); ++i) {
        printf("\t%s", (char*) stList_get(sharedContigs, i));
    }
    printf("\tweighted_mean\n");
    for (int64_t i = 0; i < numLengthScales; ++i) {
        printf("%f\t%f", decayValues[i], lengthScales[i]);
        for (int64_t j = 0; j < stList_length(sharedContigs); ++j) {
            printf("\t%f", correctnessValues[i * stList_length(sharedContigs) + j]);
        }
        printf("\t%f\n", meanCorrectnessValues[i]);
    }
    
    free(correctnessValues);
    free(meanCorrectnessValues);
    free(decayValues);
    free(lengthScales);
    stList_destruct(sharedContigs);
    stHash_destruct(truthVariants);
    stHash_destruct(queryVariants);
    
    return 0;
}

