/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"
#include <sys/stat.h>
#include <htslib/faidx.h>
#include "htsIntegration.h"
#include "bedidx.h"
#include <htslib/hts.h>
#include <htslib/vcf.h>

static char *POLISH_PARAMS_FILE = "../params/polish/ont/r9.4/allParams.np.human.r94-g360.json";
static char *PHASE_PARAMS_FILE = "../params/phase/allParams.phase_vcf.ont.json";
static char *REF_FILE = "../tests/data/realData/hg38.chr20_59M_100k.fa";
static char *BAM_FILE = "../tests/data/realData/HG002.r94g360.chr20_59M_100k.bam";
static char *VCF_FILE = "../tests/data/realData/HG002.r94g360.chr20_59M_100k.vcf";
static bool verbose = FALSE;

void checkCSV(CuTest *testCase, char *csvFile, char *sequence); // Defined in stitching test
char *getSequence(CuTest *testCase, char *outputSequenceFile, char *sequenceName);  // Defined in stitching test

int64_t marginPolishIntegrationTest(char *bamFile, char *referenceFile, char *paramsFile, char *region, char *base,
                                    bool verbose,
                                    bool diploid, bool outputRepeatCounts, bool outputPoaCsv, bool outputReadPhasingCsv,
                                    bool inMemoryBuffers) {
    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *regionStr = region == NULL ? stString_print("") : stString_print("--region %s", region);
    char *diploidString = diploid ? "--diploid" : "";
    char *outputRepeatCountsString = outputRepeatCounts ? "--outputRepeatCounts" : "";
    char *outputPoaCsvString = outputPoaCsv ? "--outputPoaCsv" : "";
    char *outputReadPhasingCsvString = outputReadPhasingCsv ? "--outputHaplotypeReads" : "";
    char *inMemoryBuffersString = inMemoryBuffers ? "" : "--tempFilesToDisk";
    char *command = stString_print("./margin polish %s %s %s %s %s %s %s %s %s %s --outputBase %s",
                                   bamFile, referenceFile, paramsFile, regionStr, logString, diploidString,
                                   outputRepeatCountsString, outputPoaCsvString, outputReadPhasingCsvString,
                                   inMemoryBuffersString, base);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(regionStr);
    free(command);
    return i;
}

static stSet *getReadNamesFromPartitionFile(CuTest *testCase, char *readPartitionFile) {
    /*
     * Parse the names of the reads from the lines of output representing the relative read phasing and return as a set
     * of strings.
     */
    stList *readLines = stFile_getLinesFromFile(readPartitionFile);
    stSet *readNames = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    for (int64_t i = 1; i < stList_length(readLines); i++) { // Start from first line after the CSV header line
        char *line = stList_get(readLines, i);
        stList *tokens = stString_splitByString(line, ",");
        CuAssertTrue(testCase, stSet_search(readNames, stList_get(tokens, 0)) ==
                               NULL); // Sanity check that read name is not present twice
        stSet_insert(readNames, stList_removeFirst(tokens)); // First field is the read name
        stList_destruct(tokens);
    }
    return readNames;
}

void test_marginPolishIntegration2(CuTest *testCase, bool inMemoryBuffers) {
    char *base = "temp_output_polish";

    // Make a temporary params file with smaller default chunk sizes
    char *tempParamsFile = "params_polish.temp";
    FILE *fh = fopen(tempParamsFile, "w");
    fprintf(fh, "{ \"include\" : \"%s\", \"polish\": { \"filterAlignmentsWithMapQBelowThisThreshold\": 0,\"chunkSize\": 20000,\"chunkBoundary\": 500 } }", POLISH_PARAMS_FILE);
    fclose(fh);

    // Run in diploid mode and get all the file outputs
    st_logInfo("\tTesting diploid polishing on %s\n", BAM_FILE);
    int64_t i = marginPolishIntegrationTest(BAM_FILE, REF_FILE, tempParamsFile, NULL, base, verbose,
                                            1, 1, 1, 1, inMemoryBuffers);
    CuAssertTrue(testCase, i == 0);

    // outputs
    char *outputHap1File = "temp_output_polish.fa.hap1";
    char *outputHap2File = "temp_output_polish.fa.hap2";
    char *outputHap1PoaFile = "temp_output_polish.poa.csv.hap1";
    char *outputHap2PoaFile = "temp_output_polish.poa.csv.hap2";
    char *outputHap1RepeatCountFile = "temp_output_polish.repeatCount.csv.hap1";
    char *outputHap2RepeatCountFile = "temp_output_polish.repeatCount.csv.hap2";
    char *outputHap1ReadPhasingFile = "temp_output_polish.reads.csv.hap1";
    char *outputHap2ReadPhasingFile = "temp_output_polish.reads.csv.hap2";

    // Parse the sequences
    char *sequenceName = "chr20";
    char *seq1 = getSequence(testCase, outputHap1File, sequenceName);
    char *seq2 = getSequence(testCase, outputHap2File, sequenceName);
    RleString *seq1Rle = rleString_construct(seq1);
    RleString *seq2Rle = rleString_construct(seq2);

    // Check the repeat counts
    checkCSV(testCase, outputHap1RepeatCountFile, seq1Rle->rleString);
    checkCSV(testCase, outputHap2RepeatCountFile, seq2Rle->rleString);

    // Check the poa outputs
    checkCSV(testCase, outputHap1PoaFile, seq1Rle->rleString);
    checkCSV(testCase, outputHap2PoaFile, seq2Rle->rleString);

    // Check the read outputs
    stSet *readsHap1 = getReadNamesFromPartitionFile(testCase, outputHap1ReadPhasingFile);
    stSet *readsHap2 = getReadNamesFromPartitionFile(testCase, outputHap2ReadPhasingFile);
    CuAssertIntEquals(testCase, 0, stSet_sizeOfIntersection(readsHap1, readsHap2));

    // Cleanup
    stFile_rmrf(outputHap1File);
    stFile_rmrf(outputHap2File);
    stFile_rmrf(outputHap1PoaFile);
    stFile_rmrf(outputHap2PoaFile);
    stFile_rmrf(outputHap1RepeatCountFile);
    stFile_rmrf(outputHap2RepeatCountFile);
    stFile_rmrf(outputHap1ReadPhasingFile);
    stFile_rmrf(outputHap2ReadPhasingFile);
    stFile_rmrf(tempParamsFile);
}

void test_marginPolishIntegration(CuTest *testCase) {
    test_marginPolishIntegration2(testCase, 0);
}

void test_marginIntegrationInMemory(CuTest *testCase) {
    test_marginPolishIntegration2(testCase, 1);
}

int64_t marginPhaseIntegrationTest(char *bamFile, char *referenceFile, char *vcfFile, char *paramsFile, char *region,
                                   char *base, bool verbose) {
    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *regionStr = region == NULL ? stString_print("") : stString_print("--region %s", region);
    char *command = stString_print("./margin phase %s %s %s %s %s %s --outputBase %s",
                                   bamFile, referenceFile, vcfFile, paramsFile, regionStr, logString, base);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(regionStr);
    free(command);
    return i;
}

void verifyHaplotaggedReads(CuTest *testCase, char *bamLocation) {
    stSet *hap1Reads = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);
    stSet *hap2Reads = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, free);

    // input file management
    samFile *in = hts_open(bamLocation, "r");
    bam_hdr_t *bamHdr = sam_hdr_read(in);
    bam1_t *aln = bam_init1();
    int r;

    // fetch alignments (either all reads or without reads)
    while ((sam_read1(in,bamHdr,aln)) >= 0) {
        char *readName = bam_get_qname(aln);
        uint8_t *tag = bam_aux_get(aln, "HP");
        if (tag == NULL)
            continue;
        int64_t hp = bam_aux2i(tag);
        if (hp == 1) {
            stSet_insert(hap1Reads, stString_copy(readName));
        } else if (hp == 2) {
            stSet_insert(hap2Reads, stString_copy(readName));
        }
    }

    // verify
    int64_t hap1Count = stSet_size(hap1Reads);
    int64_t hap2Count = stSet_size(hap2Reads);
    CuAssertTrue(testCase, hap1Count > 0);
    CuAssertTrue(testCase, hap2Count > 0);
    CuAssertTrue(testCase, hap1Count > hap2Count * 2 / 3);
    CuAssertTrue(testCase, hap2Count > hap1Count * 2 / 3);
    stSet_destruct(hap1Reads);
    stSet_destruct(hap2Reads);

    // Cleanup
    bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
    sam_close(in);
}


void verifyVcfGenotypes(CuTest *testCase, char *vcfFile1, char *vcfFile2) {
    //open vcf file
    htsFile *fp1 = hts_open(vcfFile1,"rb");
    htsFile *fp2 = hts_open(vcfFile2,"rb");

    //read header
    bcf_hdr_t *hdr1 = bcf_hdr_read(fp1);
    bcf_hdr_t *hdr2 = bcf_hdr_read(fp2);

    bcf1_t *rec1    = bcf_init();
    bcf1_t *rec2    = bcf_init();

    //save for each vcf record
    while ( true ) {
        int res1 = bcf_read(fp1, hdr1, rec1);
        int res2 = bcf_read(fp2, hdr2, rec2);
        CuAssertTrue(testCase, res1 == res2);
        if (res1 == 0) break;

        //unpack for read REF,ALT,INFO,etc
        bcf_unpack(rec1, BCF_UN_ALL);
        bcf_unpack(rec2, BCF_UN_ALL);

        // location data
        const char *chrom1 = bcf_hdr_id2name(hdr1, rec1->rid);
        const char *chrom2 = bcf_hdr_id2name(hdr2, rec2->rid);
        int64_t pos1 = rec1->pos;
        int64_t pos2 = rec2->pos;
        CuAssertTrue(testCase, stString_eq(chrom1, chrom2));
        CuAssertTrue(testCase, pos1 == pos2);

        // genotype
        int vcf1_gt1 = -1;
        int vcf1_gt2 = -1;
        int vcf2_gt1 = -1;
        int vcf2_gt2 = -1;
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(hdr1, rec1, &gt_arr, &ngt_arr);
        if (ngt>0 && !bcf_gt_is_missing(gt_arr[0])  && gt_arr[1] != bcf_int32_vector_end) {
            vcf1_gt1 = bcf_gt_allele(gt_arr[0]);
            vcf1_gt2 = bcf_gt_allele(gt_arr[1]);
        }
        ngt = bcf_get_genotypes(hdr2, rec2, &gt_arr, &ngt_arr);
        if (ngt>0 && !bcf_gt_is_missing(gt_arr[0])  && gt_arr[1] != bcf_int32_vector_end) {
            vcf2_gt1 = bcf_gt_allele(gt_arr[0]);
            vcf2_gt2 = bcf_gt_allele(gt_arr[1]);
        }
        free(gt_arr);

        // cis or trans
        bool cis = vcf1_gt1 == vcf2_gt1 && vcf1_gt2 == vcf2_gt2;
        bool trans = vcf1_gt1 == vcf2_gt2 && vcf1_gt2 == vcf2_gt1;
        CuAssertTrue(testCase, cis || trans);
    }

    // cleanup
    bcf_destroy(rec1);
    bcf_hdr_destroy(hdr1);
    bcf_destroy(rec2);
    bcf_hdr_destroy(hdr2);
    hts_close(fp1);
    hts_close(fp2);
}

void test_marginPhaseIntegration(CuTest *testCase) {
    struct stat st;
    char *base = "temp_output_phase";

    // Make a temporary params file with smaller default chunk sizes
    char *tempParamsFile = "params_phase.temp";
    FILE *fh = fopen(tempParamsFile, "w");
    fprintf(fh, "{ \"include\" : \"%s\", \"polish\": { \"chunkSize\": 20000,\"chunkBoundary\": 500 } }", PHASE_PARAMS_FILE);
    fclose(fh);

    // Run in diploid mode and get all the file outputs
    int64_t i = marginPhaseIntegrationTest(BAM_FILE, REF_FILE, VCF_FILE, tempParamsFile, NULL, base, verbose);
    CuAssertTrue(testCase, i == 0);

    // outputs
    char *outputBamFile = "temp_output_phase.haplotagged.bam";
    char *outputVcfFile = "temp_output_phase.phased.vcf";
    char *outputPhasesetFile = "temp_output_phase.phaseset.bed";

    // test bam
    CuAssertTrue(testCase, access(outputBamFile, F_OK) == 0);
    stat(outputBamFile, &st);
    CuAssertTrue(testCase, st.st_size > 0);
    verifyHaplotaggedReads(testCase, outputBamFile);

    CuAssertTrue(testCase, access(outputVcfFile, F_OK) == 0);
    stat(outputVcfFile, &st);
    CuAssertTrue(testCase, st.st_size > 0);
    verifyVcfGenotypes(testCase, outputVcfFile, VCF_FILE);

    CuAssertTrue(testCase, access(outputPhasesetFile, F_OK) == 0);
    stat(outputPhasesetFile, &st);
    CuAssertTrue(testCase, st.st_size > 0);

    // Parse the sequences
    stFile_rmrf(tempParamsFile);
    stFile_rmrf(outputVcfFile);
    stFile_rmrf(outputBamFile);
    stFile_rmrf(outputPhasesetFile);
}

CuSuite *marginIntegrationTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_marginPolishIntegration);
    SUITE_ADD_TEST(suite, test_marginIntegrationInMemory);
    SUITE_ADD_TEST(suite, test_marginPhaseIntegration);

    return suite;
}
