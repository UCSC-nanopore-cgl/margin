/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"

static char *polishParamsFile = "../params/allParams.np.json";

void checkCSV(CuTest *testCase, char *csvFile, char *sequence); // Defined in stitching test
char *getSequence(CuTest *testCase, char *outputSequenceFile, char *sequenceName);  // Defined in stitching test

int64_t
marginIntegrationTest(char *bamFile, char *referenceFile, char *paramsFile, char *region, char *base, bool verbose,
                      bool diploid,
                      bool outputRepeatCounts, bool outputPoaCsv, bool outputReadPhasingCsv) {
    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *regionStr = region == NULL ? stString_print("") : stString_print("--region %s", region);
    char *diploidString = diploid ? "--diploid" : "";
    char *outputRepeatCountsString = outputRepeatCounts ? "--outputRepeatCounts" : "";
    char *outputPoaCsvString = outputPoaCsv ? "--outputPoaCsv" : "";
    char *outputReadPhasingCsvString = outputReadPhasingCsv ? "--outputReadPhasingCsv" : "";
    char *command = stString_print("./margin %s %s %s %s %s %s %s %s %s --output %s",
                                   bamFile, referenceFile, paramsFile, regionStr, logString, diploidString,
                                   outputRepeatCountsString, outputPoaCsvString, outputReadPhasingCsvString, base);
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

void test_marginIntegration(CuTest *testCase) {
    char *referenceFile = "../tests/data/diploidTestExamples/AVG-chr7/HG002.shasta.g305.122-10980000-11086000.fasta";
    bool verbose = false;
    char *bamFile = "../tests/data/diploidTestExamples/AVG-chr7/HG002.shasta.g305.122-10980000-11086000.bam ";
    char *region = NULL;
    char *base = "temp_output";

    // Make a temporary params file with smaller default chunk sizes
    char *tempParamsFile = "params.temp";
    st_system(
            "cat %s | sed 's/\"chunkSize\": 100000/\"chunkSize\": 1000/' | sed 's/\"chunkBoundary\": 1000/\"chunkBoundary\": 500/' > %s",
            polishParamsFile, tempParamsFile);

    // Run in diploid mode and get all the file outputs
    st_logInfo("\tTesting diploid polishing on %s\n", bamFile);
    int i = marginIntegrationTest(bamFile, referenceFile, tempParamsFile, region, base, verbose,
                                  1, 1, 1, 1);
    CuAssertTrue(testCase, i == 0);

    // outputs
    char *outputHap1File = "temp_output.fa.hap1";
    char *outputHap2File = "temp_output.fa.hap2";
    char *outputHap1PoaFile = "temp_output_poa.csv.hap1";
    char *outputHap2PoaFile = "temp_output_poa.csv.hap2";
    char *outputHap1RepeatCountFile = "temp_output_repeat_counts.csv.hap1";
    char *outputHap2RepeatCountFile = "temp_output_repeat_counts.csv.hap2";
    char *outputHap1ReadPhasingFile = "temp_output_reads.csv.hap1";
    char *outputHap2ReadPhasingFile = "temp_output_reads.csv.hap2";

    // Parse the sequences
    char *sequenceName = "122:10980000-11086000";
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

CuSuite *marginIntegrationTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_marginIntegration);

    return suite;
}
