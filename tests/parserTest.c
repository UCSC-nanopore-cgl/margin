/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"
#include "htsIntegration.h"

void test_jsmnParsing(CuTest *testCase) {

    char *paramsFile = "../params/ont/r10.3/allParams.np.human.r103-g3210.json";

    //Params *allParams = params_readParams(paramsFile);
    Params *allParams = params_readParams(paramsFile);
    stRPHmmParameters *params = allParams->phaseParams;

    // Check stRPHmmParameters
    // Check a few phase parameters parsed correctly
    CuAssertTrue(testCase, params->maxNotSumTransitions);
    CuAssertIntEquals(testCase, params->maxPartitionsInAColumn, 100);
    CuAssertIntEquals(testCase, params->maxCoverageDepth, 64);
    CuAssertIntEquals(testCase, params->minReadCoverageToSupportPhasingBetweenHeterozygousSites, 2);

    // TODO: Check more parameters

    // cleanup
    stRPHmmParameters_destruct(params);
}

CuSuite *parserTestSuite(void) {
    st_setLogLevelFromString("debug");
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_jsmnParsing);

    return suite;
}
