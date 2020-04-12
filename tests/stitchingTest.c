//
// Created by Benedict Paten on 3/26/20.
//

#include "CuTest.h"
#include "margin.h"

static char *paramsFile = "../params/allParams.np.json";
static char *outputSequenceFile = "./testStitchingSequenceFile.fa";
static char *outputPoaFile = "./testStitchingPoaFile.csv";
static char *outputRepeatCountFile = "./testStitchingRepeatCountFile.csv";

void checkCSV(CuTest *testCase, char *csvFile, char *sequence) {
    /*
     * Checks the CSV file we output is organized as expected
     */
    stList *csvLines = stFile_getLinesFromFile(csvFile);
    for (int64_t i = 1; i < stList_length(csvLines); i++) {
        stList *tokens = stString_splitByString(stList_get(csvLines, i), ",");
        // Check the sequence of positions is contiguous
        int64_t j = strtol(stList_get(tokens, 0), NULL, 10);
        CuAssertIntEquals(testCase, i - 1, j);
        // Check we have the expected reference base
        if (i > 1) {
            char *refBase = stList_get(tokens, 1);
            CuAssertTrue(testCase, strlen(refBase) == 1);
            CuAssertTrue(testCase, refBase[0] == sequence[i - 2]);
        }
    }
    CuAssertIntEquals(testCase, strlen(sequence) + 2, stList_length(csvLines));
    stList_destruct(csvLines);
}

char *getSequence(CuTest *testCase, char *outputSequenceFile, char *sequenceName) {
    FILE *fh = fopen(outputSequenceFile, "r");
    stHash *seqMap = fastaReadToMap(fh);
    fclose(fh);
    stList *seqNames = stHash_getKeys(seqMap);
    stList *seqs = stHash_getValues(seqMap);
    CuAssertTrue(testCase, stList_length(seqNames) == 1);
    CuAssertStrEquals(testCase, sequenceName, stList_peek(seqNames));
    char *sequence = stString_copy(stList_pop(seqs));
    // Cleanup
    stList_destruct(seqs);
    stList_destruct(seqNames);
    stHash_destruct(seqMap);

    return sequence;
}

void test_stitching(CuTest *testCase) {
    /*
     * Runs the stitcher with a set of chunks and checks we get back the original sequence
     */

    for (int64_t test = 0; test < 10; test++) {
        // Get params
        Params *params = params_readParams(paramsFile);
        params->polishParams->useRunLengthEncoding = 0; // Turn off RLE for this test to work

        // Sequences
        char *sequence = "AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGG";
        char *sequenceName = "seq1";
        char *chunk0 = "AAAA";
        char *chunk1 = "AAAAAAAAAAT";
        char *chunk2 = "AAATTT";
        char *chunk3 = "AAATTTTTTTTTTCCCCC";
        char *chunk4 = "TTTCCCCCCCCCCG";
        char *chunk5 = "CGGGGGGGGGG";
        char *chunk6 = ""; // Test a 0 length chunk :)
        stList *chunks = stList_construct();
        stList_append(chunks, chunk0);
        stList_append(chunks, chunk1);
        stList_append(chunks, chunk2);
        stList_append(chunks, chunk3);
        stList_append(chunks, chunk4);
        stList_append(chunks, chunk5);
        stList_append(chunks, chunk6);
        stList *randomizedChunks = stList_copy(chunks, NULL);
        stList_shuffle(randomizedChunks); // This randomizes the order of the chunks
        int64_t noOfOutputChunkers = st_randomInt(1, 5);

        // Get chunker
        OutputChunkers *outputChunkers = outputChunkers_construct(noOfOutputChunkers, params,
                                                                  outputSequenceFile, outputPoaFile, NULL,
                                                                  outputRepeatCountFile,
                                                                  NULL, NULL, 0);

        // Now process the chunks
        for (int64_t i = 0; i < stList_length(randomizedChunks); i++) {
            char *chunk = stList_get(randomizedChunks, i);
            int64_t chunkOrdinal = stList_find(chunks, chunk);
            // Build the inputs for the chunker
            RleString *rle_chunk = rleString_construct_no_rle(chunk);
            Poa *poa = poa_getReferenceGraph(rle_chunk, params->polishParams->alphabet,
                                             params->polishParams->repeatSubMatrix->maximumRepeatLength);
            stList *reads = stList_construct(); // No reads, currently
            // Output the chunk
            outputChunkers_processChunkSequence(outputChunkers, st_randomInt(0, noOfOutputChunkers), chunkOrdinal,
                                                sequenceName, poa, reads);
            // Cleanup
            poa_destruct(poa);
            stList_destruct(reads);
            rleString_destruct(rle_chunk);
        }

        // Do stitching
        outputChunkers_stitch(outputChunkers, 0, stList_length(randomizedChunks));

        // Destroy chunkers
        outputChunkers_destruct(outputChunkers);

        // Read output sequence and check it agrees with what we expect
        char *seqFromFile = getSequence(testCase, outputSequenceFile, sequenceName);
        CuAssertStrEquals(testCase, sequence, seqFromFile);

        // Check poa
        checkCSV(testCase, outputPoaFile, sequence);

        // Check repeat count file
        checkCSV(testCase, outputRepeatCountFile, sequence);

        // Cleanup
        stList_destruct(randomizedChunks);
        stFile_rmrf(outputSequenceFile);
        stFile_rmrf(outputPoaFile);
        stFile_rmrf(outputRepeatCountFile);
        params_destruct(params);
        stList_destruct(chunks);
    }
}

CuSuite *stitchingTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_stitching);
    return suite;
}
