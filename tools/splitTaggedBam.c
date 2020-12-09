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
#include <htslib/thread_pool.h>
#include <htslib/faidx.h>
#include <bedidx.h>


/*
 * Main functions
 */

void stb_usage() {
    fprintf(stderr, "usage: splitTaggedBam <TAGGED_BAM> <OUT_BASE>\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Split reads with HP tag in TAGGED_BAM into OUT_H1_BAM and OUT_H2_BAM (optional OUT_H0_BAM).\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    TAGGED_BAM input file with HP present\n");
    fprintf(stderr, "    OUT_BASE basename string for output files (OUT_BASE_H1.bam, OUT_BASE_H2.bam, OUT_BASE_H0.bam\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "\nDefault options:\n");
    fprintf(stderr, "    -h --help                : Print this help screen\n");
    fprintf(stderr, "    -a --logLevel            : Set the log level [default = info]\n");
# ifdef _OPENMP
    fprintf(stderr, "    -t --threads             : Set number of concurrent threads [default = 1]\n");
#endif
    fprintf(stderr, "    -r --region              : If set, will only compute for given chromosomal region\n");
    fprintf(stderr, "                                 Format: chr:start_pos-end_pos (chr3:2000-3000)\n");
    fprintf(stderr, "\n");
}

typedef struct samview_settings {
    void *bed;
} samview_settings_t;

int main(int argc, char *argv[]) {

    // params
    char *logString = stString_copy("critical");
    char *regionStr = NULL;
    int numThreads = 1;

    if (argc < 2) {
        stb_usage();
        return 0;
    }

    // Parse the options
    while (1) {
        static struct option long_options[] = {
                { "help", no_argument, 0, 'h' },
                { "logLevel", required_argument, 0, 'a' },
# ifdef _OPENMP
                { "threads", required_argument, 0, 't'},
#endif
                { "region", required_argument, 0, 'r'},
                { 0, 0, 0, 0 } };

        int option_index = 0;
        int key = getopt_long(argc-2, &argv[2], "ha:t:r:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                free(logString);
                logString = stString_copy(optarg);
                break;
            case 'h':
                stb_usage();
                return 0;
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
                stb_usage();
                return 1;
        }
    }

    // Parameters / arguments
    char *bamInFile = stString_copy(argv[1]);
    char *bamH1File = stString_print("%s_H1.bam", argv[2]);
    char *bamH2File = stString_print("%s_H2.bam", argv[2]);
    char *bamH0File = stString_print("%s_H0.bam", argv[2]);

    // for logging
    st_setLogLevelFromString(logString);

    // sanity check (verify files are accessible)
    if (access(bamInFile, R_OK) != 0) {
        st_errAbort("Could not read from input bam file: %s\n", bamInFile);
    }

    // counting
    int64_t h1Count = 0;
    int64_t h2Count = 0;
    int64_t h0Count = 0;
    int64_t hUnCount = 0;

    // input file management
    samFile *in = hts_open(bamInFile, "r");
    st_logCritical("Reading from BAM: %s \n", bamInFile);
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }

    // bam index
    hts_idx_t *idx = NULL;
    if ((idx = sam_index_load(in, bamInFile)) == 0) {
        st_errAbort("ERROR: Cannot open index for bam file %s\n", bamInFile);
    }

    // reading data structures
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    int r;

    // thread stuff
    htsThreadPool threadPool = {NULL, 0};
    if (!(threadPool.pool = hts_tpool_init(numThreads))) {
        fprintf(stderr, "Error creating thread pool\n");
    }
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &threadPool);

    // prep for index (not entirely sure what all this does.  see samtools/sam_view.c
    hts_itr_multi_t *iter = NULL;
    int filter_state = ALL, filter_op = 0;
    samview_settings_t settings = {.bed = NULL};
    char *region[1] = {};
    int regcount = 0;
    hts_reglist_t *reglist = NULL;
    if (regionStr != NULL) {
        region[0] = stString_copy(regionStr);
        settings.bed = bed_hash_regions(settings.bed, region, 0, 1,
                                        &filter_op); //insert(1) or filter out(0) the regions from the command line in the same hash table as the bed file
        if (!filter_op) filter_state = FILTERED;
        reglist = bed_reglist(settings.bed, filter_state, &regcount);
        if (!reglist) {
            st_errAbort("ERROR: Could not create list of regions for read conversion");
        }
        if ((iter = sam_itr_regions(idx, bamHdr, reglist, regcount)) == 0) {
            st_errAbort("ERROR: Cannot open iterator for region %s for bam file %s\n", region[0], bamInFile);
        }
        st_logCritical("Reading from region %s\n", regionStr);
    }

    // open bam to write to
    st_logCritical("Writing to BAMs: H1:%s, H2:%s, H0:%s \n", bamH1File, bamH2File, bamH0File == NULL ? "" : bamH0File);
    samFile *outH1 = hts_open(bamH1File, "wb");
    samFile *outH2 = hts_open(bamH2File, "wb");
    samFile *outH0 = hts_open(bamH0File, "wb");
    r = sam_hdr_write(outH1, bamHdr);
    r = sam_hdr_write(outH2, bamHdr);
    r = sam_hdr_write(outH0, bamHdr);

    // fetch alignments
    while ((regionStr == NULL ? sam_read1(in, bamHdr, aln) : sam_itr_multi_next(in, iter, aln)) > 0) {
        uint8_t *tag = bam_aux_get(aln, "HP");
        if (tag == NULL) {
            hUnCount++;
            continue;
        }
        int64_t hp = bam_aux2i(tag);

        switch (hp) {
            case 1:
                h1Count++;
                r = sam_write1(outH1, bamHdr, aln);
                break;
            case 2:
                h2Count++;
                r = sam_write1(outH2, bamHdr, aln);
                break;
            case 0:
                h0Count++;
                r = sam_write1(outH0, bamHdr, aln);
                break;
            default:
                st_logInfo("Unrecognized HP tag value: %"PRId64"\n", hp);
        }
    }
    st_logCritical("Wrote reads with divisions: H1 %"PRId64", H2 %"PRId64", H0 %"PRId64", untagged %"PRId64"\n",
               h1Count, h2Count, h0Count, hUnCount);

    // Cleanup

    if (regionStr != NULL) {
        hts_itr_multi_destroy(iter);
        free(region[0]);
        bed_destroy(settings.bed);
    }
    hts_idx_destroy(idx);
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    hts_tpool_destroy(threadPool.pool);
    sam_close(in);
    sam_close(outH1);
    sam_close(outH2);
    sam_close(outH0);

    st_logInfo("Fin.");

    return 0;
}

