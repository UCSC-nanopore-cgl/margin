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
    fprintf(stderr, "usage: marginTagBam <TAGGED_BAM> <OUT_H1_BAM> <OUT_H2_BAM> [OUT_H0_BAM]\n");
    fprintf(stderr, "Version: %s \n\n", MARGIN_POLISH_VERSION_H);
    fprintf(stderr, "Split reads with HP tag in TAGGED_BAM into OUT_H1_BAM and OUT_H2_BAM (optional OUT_H0_BAM).\n");

    fprintf(stderr, "\nRequired arguments:\n");
    fprintf(stderr, "    TAGGED_BAM input file with HP present\n");
    fprintf(stderr, "    OUT_H1_BAM, OUT_H2_BAM filenames for output bams with split reads");
    fprintf(stderr, "    OUT_H0_BAM will write reads with H0 tag and untagged reads\n");
    fprintf(stderr, "\n");
}


int main(int argc, char *argv[]) {

    // Parameters / arguments
    char *bamInFile = NULL;
    char *bamH1File = NULL;
    char *bamH2File = NULL;
    char *bamH0File = NULL;

    if (argc < 4) {
        usage();
        return 0;
    }

    bamInFile = stString_copy(argv[1]);
    bamH1File = stString_copy(argv[2]);
    bamH2File = stString_copy(argv[3]);
    if (stString_eq(bamH1File, bamH2File)) {
        st_errAbort("Output BAMs must have different names");
    }
    if (argc > 4) {
        bamH0File = stString_copy(argv[4]);
        if (stString_eq(bamH0File, bamH1File) || stString_eq(bamH0File, bamH2File)) {
            st_errAbort("Output BAMs must have different names (including H0 bam)");
        }
    }

    // for logging
    st_setLogLevel(critical);

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
    hts_itr_multi_t *iter = NULL;
    if (in == NULL) {
        st_errAbort("ERROR: Cannot open bam file %s\n", bamInFile);
    }
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    int r;

    // open bam to write to
    st_logCritical("Writing to BAM: %s \n", bamOutFile);
    samFile *out = hts_open(bamOutFile, "wb");
    r = sam_hdr_write(out, bamHdr);

    // fetch alignments
    while (sam_read1(in,bamHdr,aln) >= 0) {
        // get haplotag value
        char *readName = bam_get_qname(aln);
        int32_t hp = 0;
        if (stHash_search(readInfo, readName) == NULL) {
            st_logInfo("Did not find %s in read info\n", readName);
            hUnCount++;
        } else {
            int64_t foundTag = (int64_t) stHash_search(readInfo, readName);
            if (foundTag == -1) {
                h0Count++;
            } else if (foundTag == 1) {
                hp = 1;
                h1Count++;
            } else if (foundTag == 2) {
                hp = 2;
                h2Count++;
            } else {
                st_errAbort("Got haplotag value %"PRId64" for read %s!", foundTag, readName);
            }
        }

        if (bam_aux_get(aln, "HP") != NULL) {
            bam_aux_update_int(aln, "HP", hp);
        } else {
            bam_aux_append(aln, "HP", 'i', sizeof(hp), (uint8_t*) &hp);
        }
        r = sam_write1(out, bamHdr, aln);
    }
    st_logCritical("Wrote reads with divisions: H1 %"PRId64", H2 %"PRId64", and H0 %"PRId64"\n",
               h1Count, h2Count, h0Count);
    st_logCritical("Found %"PRId64" reads which were not annotated in info file (tagged as H0, but not counted above).\n", hUnCount);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
    sam_close(out);
    stHash_destruct(readInfo);

    st_logInfo("Fin.");

    return 0;
}

