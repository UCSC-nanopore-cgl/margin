/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <htsIntegration.h>
#include <sys/stat.h>
#include "CuTest.h"
#include "margin.h"
#include "helenFeatures.h"

static char *FEATURE_TEST_PARAMS = "../params/allParams.np.json";
static char *FEATURE_TEST_NO_RLE_PARAMS = "../params/allParams.np.no_rle.json";
static char *FEATURE_TEST_BAM = "../tests/data/featureTest/featureTest.bam";
static char *FEATURE_TEST_FA = "../tests/data/featureTest/featureTest.fa";
static char *FEATURE_TEST_TRUTH_BAM = "../tests/data/featureTest/featureTestTruth.bam";
static char *FEATURE_TEST_TRUTH_SEQ = "ACGATAACCGGTTAAACCCCGGGTTTCAAACCCCGGGGTTGATTACACAT";


int64_t polishingFeatureTest(char *bamFile, char *referenceFile, char *paramsFile, char *featureType,
                             char *featureTruth, char *outputName, bool verbose) {

    // Run margin phase
    char *logString = verbose ? "--logLevel DEBUG" : "--logLevel INFO";
    char *featureTruthCmd =
            featureTruth == NULL ? stString_copy("") : stString_print("--trueReferenceBam %s", featureTruth);
    char *featureTypeCmd = featureType == NULL ? stString_copy("-f") : stString_print("-F %s", featureType);
    char *command = stString_print("./marginPolish %s %s %s --outputBase %s %s %s %s",
                                   bamFile, referenceFile, paramsFile, outputName, logString, featureTypeCmd,
                                   featureTruthCmd);
    st_logInfo("> Running command: %s\n", command);
    fprintf(stderr, "> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(featureTruthCmd);
    free(command);
    return i;
}


void test_defaultFeatureGeneration(CuTest *testCase) {

    struct stat st;
    char *outputName = "test.default";
    char *expectedOutputFa = stString_print("%s.fa", outputName);
    char *expectedOutputFeature = stString_print("%s.T00.h5", outputName);

    CuAssertTrue(testCase, access(expectedOutputFa, R_OK) == 0 || remove(expectedOutputFa) != 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, R_OK) == 0 || remove(expectedOutputFeature) != 0);

    st_logInfo("\n\nTesting default feature polishing on %s and %s\n", FEATURE_TEST_BAM, FEATURE_TEST_FA);
    int64_t retCode = polishingFeatureTest(FEATURE_TEST_BAM, FEATURE_TEST_FA, FEATURE_TEST_PARAMS, NULL,
                                           NULL, outputName, FALSE);

    CuAssertTrue(testCase, retCode == 0);
    CuAssertTrue(testCase, access(expectedOutputFa, F_OK) == 0);
    stat(expectedOutputFa, &st);
    CuAssertTrue(testCase, st.st_size > 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, F_OK) == 0);
    stat(expectedOutputFeature, &st);
    CuAssertTrue(testCase, st.st_size > 0);
}

void test_simpleWeightFeatureGeneration(CuTest *testCase) {

    struct stat st;
    char *featureType = "simpleWeight";
    char *outputName = "test.simpleWeightFeature";
    char *expectedOutputFa = stString_print("%s.fa", outputName);
    char *expectedOutputFeature = stString_print("%s.T00.h5", outputName);

    CuAssertTrue(testCase, access(expectedOutputFa, R_OK) == 0 || remove(expectedOutputFa) != 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, R_OK) == 0 || remove(expectedOutputFeature) != 0);

    st_logInfo("\n\nTesting %s feature polishing on %s and %s\n", featureType, FEATURE_TEST_BAM, FEATURE_TEST_FA);
    int64_t retCode = polishingFeatureTest(FEATURE_TEST_BAM, FEATURE_TEST_FA, FEATURE_TEST_NO_RLE_PARAMS,
                                           featureType, NULL, outputName, FALSE);

    CuAssertTrue(testCase, retCode == 0);
    CuAssertTrue(testCase, access(expectedOutputFa, F_OK) == 0);
    stat(expectedOutputFa, &st);
    CuAssertTrue(testCase, st.st_size > 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, F_OK) == 0);
    stat(expectedOutputFeature, &st);
    CuAssertTrue(testCase, st.st_size > 0);
}


void test_splitRleWeightFeatureGeneration(CuTest *testCase) {

    struct stat st;
    char *featureType = "splitRleWeight";
    char *outputName = "test.splitRleWeightFeature";
    char *expectedOutputFa = stString_print("%s.fa", outputName);
    char *expectedOutputFeature = stString_print("%s.T00.h5", outputName);

    CuAssertTrue(testCase, access(expectedOutputFa, R_OK) == 0 || remove(expectedOutputFa) != 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, R_OK) == 0 || remove(expectedOutputFeature) != 0);

    st_logInfo("\n\nTesting %s feature polishing on %s and %s\n", featureType, FEATURE_TEST_BAM, FEATURE_TEST_FA);
    int64_t retCode = polishingFeatureTest(FEATURE_TEST_BAM, FEATURE_TEST_FA, FEATURE_TEST_PARAMS, featureType,
                                           NULL, outputName, FALSE);

    CuAssertTrue(testCase, retCode == 0);
    CuAssertTrue(testCase, access(expectedOutputFa, F_OK) == 0);
    stat(expectedOutputFa, &st);
    CuAssertTrue(testCase, st.st_size > 0);
    CuAssertTrue(testCase, access(expectedOutputFeature, F_OK) == 0);
    stat(expectedOutputFeature, &st);
    CuAssertTrue(testCase, st.st_size > 0);
}

void test_simpleWeightIndex(CuTest *testCase) {

    int idx;
    PoaFeatureSimpleWeight *feature = PoaFeature_SimpleWeight_construct(0, 0);

    for (int64_t c = 0; c < SYMBOL_NUMBER_NO_N; c++) {
        idx = PoaFeature_SimpleWeight_charIndex((Symbol) c, TRUE);
        CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
        feature->weights[idx] += 1;

        idx = PoaFeature_SimpleWeight_charIndex((Symbol) c, FALSE);
        CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
        feature->weights[idx] += 1;
    }

    idx = PoaFeature_SimpleWeight_gapIndex(TRUE);
    CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    feature->weights[idx] += 1;

    idx = PoaFeature_SimpleWeight_gapIndex(FALSE);
    CuAssertTrue(testCase, idx < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE);
    feature->weights[idx] += 1;

    for (int64_t i = 0; i < POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE; i++) {
        CuAssertTrue(testCase, feature->weights[i] == 1);
    }

    PoaFeature_SimpleWeight_destruct(feature);

}


void test_splitRleWeightIndex(CuTest *testCase) {

    int idx;
    int maxRunLength = POAFEATURE_SPLIT_MAX_RUN_LENGTH_DEFAULT;
    int maxIndex = ((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2;
    PoaFeatureSplitRleWeight *feature = PoaFeature_SplitRleWeight_construct(0, 0, 0, maxRunLength);

    for (int64_t c = 0; c < SYMBOL_NUMBER - 1; c++) {
        for (int64_t l = 0; l <= maxRunLength; l++) {
            idx = PoaFeature_SplitRleWeight_charIndex(maxRunLength, (Symbol) c, l, TRUE);
            CuAssertTrue(testCase, idx < maxIndex);
            feature->weights[idx] += 1;

            idx = PoaFeature_SplitRleWeight_charIndex(maxRunLength, (Symbol) c, l, FALSE);
            CuAssertTrue(testCase, idx < maxIndex);
            feature->weights[idx] += 1;
        }
    }

    idx = PoaFeature_SplitRleWeight_gapIndex(maxRunLength, TRUE);
    CuAssertTrue(testCase, idx < maxIndex);
    feature->weights[idx] += 1;

    idx = PoaFeature_SplitRleWeight_gapIndex(maxRunLength, FALSE);
    CuAssertTrue(testCase, idx < maxIndex);
    feature->weights[idx] += 1;

    for (int64_t i = 0; i < maxIndex; i++) {
        if (feature->weights[i] != 1) {
            CuAssertTrue(testCase, feature->weights[i] == 1);
        }
    }

    PoaFeature_SplitRleWeight_destruct(feature);
}


void test_channelRleWeightIndex(CuTest *testCase) {

//    int PoaFeature_ChannelRleWeight_charNuclIndex(Symbol character, bool forward);
//    int PoaFeature_ChannelRleWeight_gapNuclIndex(bool forward);
//    int PoaFeature_ChannelRleWeight_charRLIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward);

    int idx;
    int maxRunLength = POAFEATURE_CHANNEL_MAX_RUN_LENGTH_DEFAULT;
    int maxNuclIndex = (SYMBOL_NUMBER) * 2;
    int maxRLIndex = (SYMBOL_NUMBER - 1) * (1 + maxRunLength) * 2;
    PoaFeatureChannelRleWeight *feature = PoaFeature_ChannelRleWeight_construct(0, 0, 0, maxRunLength);

    for (int64_t c = 0; c < SYMBOL_NUMBER - 1; c++) {
        idx = PoaFeature_ChannelRleWeight_charNuclIndex(c, TRUE);
        CuAssertTrue(testCase, idx < maxNuclIndex);
        feature->nucleotideWeights[idx] += 1;

        idx = PoaFeature_ChannelRleWeight_charNuclIndex(c, FALSE);
        CuAssertTrue(testCase, idx < maxNuclIndex);
        feature->nucleotideWeights[idx] += 1;

        for (int64_t l = 0; l <= maxRunLength; l++) {
            idx = PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, (Symbol) c, l, TRUE);
            CuAssertTrue(testCase, idx < maxRLIndex);
            feature->runLengthWeights[idx] += 1;

            idx = PoaFeature_ChannelRleWeight_charRLIndex(maxRunLength, (Symbol) c, l, FALSE);
            CuAssertTrue(testCase, idx < maxRLIndex);
            feature->runLengthWeights[idx] += 1;
        }
    }

    idx = PoaFeature_ChannelRleWeight_gapNuclIndex(TRUE);
    CuAssertTrue(testCase, idx < maxNuclIndex);
    feature->nucleotideWeights[idx] += 1;

    idx = PoaFeature_ChannelRleWeight_gapNuclIndex(FALSE);
    CuAssertTrue(testCase, idx < maxNuclIndex);
    feature->nucleotideWeights[idx] += 1;

    for (int64_t i = 0; i < maxNuclIndex; i++) {
        if (feature->nucleotideWeights[i] != 1) {
            CuAssertTrue(testCase, feature->nucleotideWeights[i] == 1);
        }
    }

    for (int64_t i = 0; i < maxRLIndex; i++) {
        if (feature->runLengthWeights[i] != 1) {
            CuAssertTrue(testCase, feature->runLengthWeights[i] == 1);
        }
    }

    PoaFeature_ChannelRleWeight_destruct(feature);
}


void test_diploidRleWeightIndex(CuTest *testCase) {

    int idx;
    int maxRunLength = POAFEATURE_DIPLOID_MAX_RUN_LENGTH_DEFAULT;
    int maxIndex = ((SYMBOL_NUMBER - 1) * (1 + maxRunLength) + 1) * 2;
    PoaFeatureDiploidRleWeight *feature = PoaFeature_DiploidRleWeight_construct(0, 0, 0, maxRunLength);

    for (int64_t c = 0; c < SYMBOL_NUMBER - 1; c++) {
        for (int64_t l = 0; l <= maxRunLength; l++) {
            idx = PoaFeature_DiploidRleWeight_charIndex(maxRunLength, (Symbol) c, l, TRUE);
            CuAssertTrue(testCase, idx < maxIndex);
            feature->weightsHOn[idx] += 1;

            idx = PoaFeature_DiploidRleWeight_charIndex(maxRunLength, (Symbol) c, l, FALSE);
            CuAssertTrue(testCase, idx < maxIndex);
            feature->weightsHOn[idx] += 1;
        }
    }

    idx = PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, TRUE);
    CuAssertTrue(testCase, idx < maxIndex);
    feature->weightsHOn[idx] += 1;

    idx = PoaFeature_DiploidRleWeight_gapIndex(maxRunLength, FALSE);
    CuAssertTrue(testCase, idx < maxIndex);
    feature->weightsHOn[idx] += 1;

    for (int64_t i = 0; i < maxIndex; i++) {
        if (feature->weightsHOn[i] != 1) {
            CuAssertTrue(testCase, feature->weightsHOn[i] == 1);
        }
    }

    PoaFeature_DiploidRleWeight_destruct(feature);
}

CuSuite* featureTestSuite(void) {
    CuSuite* suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_simpleWeightIndex);
    SUITE_ADD_TEST(suite, test_splitRleWeightIndex);
    SUITE_ADD_TEST(suite, test_channelRleWeightIndex);
    SUITE_ADD_TEST(suite, test_diploidRleWeightIndex);
    SUITE_ADD_TEST(suite, test_defaultFeatureGeneration);
    SUITE_ADD_TEST(suite, test_simpleWeightFeatureGeneration);
    SUITE_ADD_TEST(suite, test_splitRleWeightFeatureGeneration);

    return suite;
}
