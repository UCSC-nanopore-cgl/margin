/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <math.h>

#include "CuTest.h"
#include "margin.h"
#include "localPhasingCorrectness.h"

static char *TEST_VCF = "../tests/data/localPhasingCorrectness/smallPhased.vcf";

double directLPC(stList *queryPhasedVariants, stList *truthPhasedVariants, double decay) {
    
    double lpcNumer = 0.0;
    double lpcDenom = 0.0;
    
    // i'll assume that the variants are identical between the two lists for
    // this dumber version of the algorithm, and also assume that ref and alt are the
    // only variants, and that every variant is heterozygous (which should also
    // be true in the full algorithm)
    assert(stList_length(queryPhasedVariants) == stList_length(truthPhasedVariants));
    
    for (int64_t i = 0; i < stList_length(queryPhasedVariants); ++i) {
        for (int64_t j = 0; j < stList_length(queryPhasedVariants); ++j) {
            if (i == j) {
                continue;
            }
            
            PhasedVariant *qpvi = stList_get(queryPhasedVariants, i);
            PhasedVariant *qpvj = stList_get(queryPhasedVariants, j);
            PhasedVariant *tpvi = stList_get(truthPhasedVariants, i);
            PhasedVariant *tpvj = stList_get(truthPhasedVariants, j);
            
            double summand = pow(decay, fabs((double) (i - j)));
            
            lpcDenom += summand;
            
            if (!stString_eq(qpvi->phaseSet, qpvj->phaseSet) ||
                !stString_eq(tpvi->phaseSet, tpvj->phaseSet) ||
                (qpvi->gt1 == tpvi->gt1) == (qpvj->gt1 == tpvj->gt1)) {
                // these are either separated by a phase break or have a consistent
                // relative phasing, so we count this pair as correct;
                lpcNumer += summand;
            }
        }
    }
    
    return lpcNumer / lpcDenom;
}




int64_t lpcIntegrationTest(char *truthVcfFile, char *queryVcfFile) {
    // Run localPhasingCorrectness
    char *command = stString_print("./calcLocalPhasingCorrectness %s %s", truthVcfFile, queryVcfFile);
    st_logInfo("> Running command: %s\n", command);

    int64_t i = st_system(command);
    free(command);
    return i;
}

void test_correctValue(CuTest *testCase) {

    
    int64_t numSites = 5;
    stList *alleles[2];
    for (int64_t i = 0; i < 2; ++i) {
        alleles[i] = stList_construct();
        for (int64_t j = 0; j < numSites; ++j) {
            stList *a = stList_construct3(0, free);
            stList_append(a, stString_copy("A"));
            stList_append(a, stString_copy("C"));
            stList_append(alleles[i], a);
        }
    }
    
    stList *variants[2];
    for (int64_t i = 0; i < 2; ++i) {
        variants[i] = stList_construct3(0, (void (*)(void *)) partialPhaseSums_destruct);
        for (int64_t j = 0; j < numSites; ++i) {
            stList_append(variants[i],
                          phasedVariant_construct("ref", 1, 60.0,
                                                  stList_get(alleles[i], j),
                                                  0, 1, stString_copy("ps")));
        }
    }
    
    double decayValues[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    
    for (int64_t i = 0; i < 11; ++i) {
        int64_t numPhased;
        double correctness = phasingCorrectness(variants[0], variants[1], decayValues[i],
                                                &numPhased);
        double direct = directLPC(variants[0], variants[1], decayValues[i]);
        CuAssertTrue(testCase, numPhased == numSites);
        
        CuAssertDblEquals(testCase, correctness, direct, 0.000001);
        // everything in this test is perfectly phased
        CuAssertDblEquals(testCase, correctness, 1.0, 0.000001);
    }
    
    for (int64_t i = 0; i < 2; ++i) {
        stList_destruct(variants[i]);
        stList_destruct(alleles[i]);
    }
}

void test_executableExecutes(CuTest *testCase) {
    int64_t ret = lpcIntegrationTest(TEST_VCF, TEST_VCF);
    CuAssertTrue(testCase, ret == 0);
}



CuSuite *lpcTestSuite(void) {
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_executableExecutes);
    SUITE_ADD_TEST(suite, test_correctValue);

    return suite;
}
