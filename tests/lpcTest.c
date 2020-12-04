/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "margin.h"
#include "localPhasingCorrectness.h"



int64_t lpcIntegrationTest(char *truthVcfFile, char *queryVcfFile) {
    // Run localPhasingCorrectness
    char *command = stString_print("./localPhasingCorrectness %s %s", truthVcfFile, queryVcfFile);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(command);
    return i;
}

void test_executableExecutes(CuTest *testCase) {
    int64_t ret = lpcIntegrationTest("", "");
    CuAssertTrue(testCase, ret == 0);
}

CuSuite *lpcTestSuite(void) {
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_executableExecutes);

    return suite;
}