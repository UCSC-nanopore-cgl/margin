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

// in lieu of a separate algorithm for switch correctness, i'll just have an independent
// implementation
double redundantSwitchCorrectness(stList *queryPhasedVariants, stList *truthPhasedVariants,
                                  bool crossBlockCorrect) {
    
    
    int64_t numUnswitched = 0;
    int64_t numPossibleUnswitched = 0;
    
    for (int64_t i = 1; i < stList_length(queryPhasedVariants); ++i) {
        
        PhasedVariant *qpv_prev = stList_get(queryPhasedVariants, i - 1);
        PhasedVariant *tpv_prev = stList_get(truthPhasedVariants, i - 1);
        PhasedVariant *qpv = stList_get(queryPhasedVariants, i);
        PhasedVariant *tpv = stList_get(truthPhasedVariants, i);
        
        if (stString_eq(qpv_prev->phaseSet, qpv->phaseSet) &&
            stString_eq(tpv_prev->phaseSet, tpv->phaseSet)) {
            if ((qpv_prev->gt1 == tpv_prev->gt1) == (qpv->gt1 == tpv->gt1)) {
                ++numUnswitched;
            }
            ++numPossibleUnswitched;
        }
        else if (crossBlockCorrect) {
            ++numUnswitched;
            ++numPossibleUnswitched;
        }
    }
    
    return ((double) numUnswitched) / numPossibleUnswitched;
}

double redundantBaseSwitchCorrectness(stList *queryPhasedVariants, stList *truthPhasedVariants,
                                      bool crossBlockCorrect) {
    
    
    int64_t numUnswitched = 0;
    
    int64_t minDist = INT_MAX;
    for (int64_t i = 1; i < stList_length(queryPhasedVariants); ++i) {
        
        PhasedVariant *qpv_prev = stList_get(queryPhasedVariants, i - 1);
        PhasedVariant *tpv_prev = stList_get(truthPhasedVariants, i - 1);
        PhasedVariant *qpv = stList_get(queryPhasedVariants, i);
        PhasedVariant *tpv = stList_get(truthPhasedVariants, i);
        
        if (qpv->refPos - qpv_prev->refPos < minDist &&
            ((stString_eq(qpv_prev->phaseSet, qpv->phaseSet) &&
              stString_eq(tpv_prev->phaseSet, tpv->phaseSet))
             || crossBlockCorrect)) {
            minDist = qpv->refPos - qpv_prev->refPos;
        }
    }
    
    st_logDebug("calc min dist at %"PRId64"\n", minDist);
    
    int64_t numAtMinDist = 0;
    for (int64_t i = 1; i < stList_length(queryPhasedVariants); ++i) {
        
        st_logDebug("iter i = %"PRId64"\n", i);
        
        PhasedVariant *qpv_prev = stList_get(queryPhasedVariants, i - 1);
        PhasedVariant *tpv_prev = stList_get(truthPhasedVariants, i - 1);
        PhasedVariant *qpv = stList_get(queryPhasedVariants, i);
        PhasedVariant *tpv = stList_get(truthPhasedVariants, i);
        
        if (qpv->refPos - qpv_prev->refPos == minDist) {
            if (stString_eq(qpv_prev->phaseSet, qpv->phaseSet) &&
                stString_eq(tpv_prev->phaseSet, tpv->phaseSet)) {
                if ((qpv_prev->gt1 == tpv_prev->gt1) == (qpv->gt1 == tpv->gt1)) {
                    
                    ++numUnswitched;
                }
                ++numAtMinDist;
            }
            else if (crossBlockCorrect) {
                ++numUnswitched;
                ++numAtMinDist;
            }
            
            st_logDebug("num unswitched %"PRId64", num min dist %"PRId64"\n", numUnswitched, numAtMinDist);
        }
    }
    
    return ((double) numUnswitched) / numAtMinDist;
}

// a less efficient but simpler algorithm for LPC
double directLPC(stList *queryPhasedVariants, stList *truthPhasedVariants, double decay,
                 bool bySeqDist, bool crossBlockCorrect) {
    
    // i'll assume that the variants are identical between the two lists for
    // this dumber version of the algorithm, and also assume that ref and alt are the
    // only variants, and that every variant is heterozygous (which should also
    // be true in the full algorithm)
    
    assert(stList_length(queryPhasedVariants) == stList_length(truthPhasedVariants));
    
    if (decay == 0.0) {
        if (bySeqDist) {
            return redundantBaseSwitchCorrectness(queryPhasedVariants, truthPhasedVariants,
                                                  crossBlockCorrect);
        }
        else {
            return redundantSwitchCorrectness(queryPhasedVariants, truthPhasedVariants,
                                              crossBlockCorrect);
        }
    }
    
    double lpcNumer = 0.0;
    double lpcDenom = 0.0;
    
    for (int64_t i = 0; i < stList_length(queryPhasedVariants); ++i) {
        
        for (int64_t j = 0; j < stList_length(queryPhasedVariants); ++j) {
            if (i == j) {
                continue;
            }
            
            //st_logDebug("i = %"PRId64", j = %"PRId64"\n", i, j);
            
            PhasedVariant *qpvi = stList_get(queryPhasedVariants, i);
            PhasedVariant *qpvj = stList_get(queryPhasedVariants, j);
            PhasedVariant *tpvi = stList_get(truthPhasedVariants, i);
            PhasedVariant *tpvj = stList_get(truthPhasedVariants, j);
            
            if ((!stString_eq(qpvi->phaseSet, qpvj->phaseSet)
                 || !stString_eq(tpvi->phaseSet, tpvj->phaseSet))
                && !crossBlockCorrect) {
                continue;
            }
            
            double summand;
            if (bySeqDist) {
                summand = pow(decay, fabs((double) (qpvi->refPos - qpvj->refPos)));
            }
            else {
                summand = pow(decay, fabs((double) (i - j)));
            }
            
            lpcDenom += summand;
            
            if (!stString_eq(qpvi->phaseSet, qpvj->phaseSet) ||
                !stString_eq(tpvi->phaseSet, tpvj->phaseSet) ||
                (qpvi->gt1 == tpvi->gt1) == (qpvj->gt1 == tpvj->gt1)) {
                // these are either separated by a phase break or have a consistent
                // relative phasing, so we count this pair as correct;
                //st_logDebug("\tphased consistently\n");
                lpcNumer += summand;
            }
            //st_logDebug("\tterm %f, numer %f, denom %f numer\n", summand, lpcNumer, lpcDenom);
        }
    }
    
    return lpcNumer / lpcDenom;
}


void test_correctValueSimple(CuTest *testCase) {
    
    int64_t numDecayValues = 11;
    double decayValues[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    
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
    
    // all perfectly phased variants
    stList *variants[2];
    for (int64_t i = 0; i < 2; ++i) {
        variants[i] = stList_construct3(0, (void (*)(void *)) phasedVariant_destruct);
        for (int64_t j = 0; j < numSites; ++j) {
            stList_append(variants[i],
                          phasedVariant_construct("ref", j * j + 1, 60.0,
                                                  stList_get(alleles[i], j),
                                                  0, 1, "ps"));
        }
    }
    
    st_logDebug("TEST BLOCK: perfect phasing\n");
    for (int64_t i = 0; i < numDecayValues; ++i) {
        for (int64_t byDist = 0; byDist < 2; ++byDist) {
            for (int64_t crossBlockCorrect = 0; crossBlockCorrect < 2; ++crossBlockCorrect) {
                bool doByDist = byDist;
                bool doCrossBlockCorrect = crossBlockCorrect;
                
                int64_t length;
                double correctness = phasingCorrectness(variants[0], variants[1], decayValues[i],
                                                        doByDist, doCrossBlockCorrect, &length);
                double direct = directLPC(variants[0], variants[1], decayValues[i], doByDist,
                                          doCrossBlockCorrect);
                st_logDebug("decay %f, by dist %d, cross block correct %d, algorithm %f, direct %f\n", decayValues[i], doByDist, doCrossBlockCorrect, correctness, direct);
                CuAssertDblEquals(testCase, correctness, direct, 0.000001);
                // everything in this test is perfectly phased
                CuAssertDblEquals(testCase, correctness, 1.0, 0.000001);
                
                if (doByDist) {
                    // the distance between
                    CuAssertTrue(testCase, length == (numSites - 1) * (numSites - 1));
                }
                else {
                    CuAssertTrue(testCase, length == numSites);
                }
            }
        }
    }
    
    // set up new variants with errors
    for (int64_t i = 0; i < 2; ++i) {
        for (int64_t j = 0; j < numSites; ++j) {
            int64_t gt1, gt2;
            if (i == 0 || j % 2 == 0) {
                gt1 = 0;
                gt2 = 1;
            }
            else {
                gt1 = 1;
                gt2 = 0;
            }
            PhasedVariant *var = stList_get(variants[i], j);
            var->gt1 = gt1;
            var->gt2 = gt2;
        }
    }
    
    st_logDebug("TEST BLOCK: phasing with errors\n");
    for (int64_t i = 0; i < numDecayValues; ++i) {
        for (int64_t byDist = 0; byDist < 2; ++byDist) {
            for (int64_t crossBlockCorrect = 0; crossBlockCorrect < 2; ++crossBlockCorrect) {
                bool doByDist = byDist;
                bool doCrossBlockCorrect = crossBlockCorrect;
                
                int64_t length;
                double correctness = phasingCorrectness(variants[0], variants[1], decayValues[i],
                                                        doByDist, doCrossBlockCorrect, &length);
                double direct = directLPC(variants[0], variants[1], decayValues[i], doByDist,
                                          doCrossBlockCorrect);
                st_logDebug("decay %f, by dist %d, cross block correct %d, algorithm %f, direct %f\n", decayValues[i], doByDist, doCrossBlockCorrect, correctness, direct);
                CuAssertDblEquals(testCase, correctness, direct, 0.000001);
                // everything in this test is perfectly phased
                CuAssertTrue(testCase, correctness >= 0.0);
                CuAssertTrue(testCase, correctness < 1.00000001); // a little room for numerical slop
                
                if (doByDist) {
                    // the distance between
                    CuAssertTrue(testCase, length == (numSites - 1) * (numSites - 1));
                }
                else {
                    CuAssertTrue(testCase, length == numSites);
                }
            }
        }
    }
    
    
    // clean up
    for (int64_t i = 0; i < 2; ++i) {
        stList_destruct(variants[i]);
        stList_destruct(alleles[i]);
    }
}

void test_correctValueWithPhaseSets(CuTest *testCase) {
    
    int64_t numDecayValues = 11;
    double decayValues[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    
    int64_t numSites = 10;
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
    
    // all perfectly phased variants
    stList *variants[2];
    for (int64_t i = 0; i < 2; ++i) {
        variants[i] = stList_construct3(0, (void (*)(void *)) phasedVariant_destruct);
        for (int64_t j = 0; j < numSites; ++j) {
            int64_t gt1, gt2;
            if (i == 0 || j < 4) {
                gt1 = 0;
                gt2 = 1;
            }
            else {
                gt1 = 1;
                gt2 = 0;
            }
            const char *phaseSet = (i == 0 || j < 4) ? "ps1" : (j < 8 ? "ps2" :"ps3");
            stList_append(variants[i],
                          phasedVariant_construct("ref", j * j + 1, 60.0,
                                                  stList_get(alleles[i], j),
                                                  gt1, gt2, phaseSet));
        }
    }
    
    st_logDebug("TEST BLOCK: perfect phasing in two blocks\n");
    for (int64_t i = 0; i < numDecayValues; ++i) {
        for (int64_t byDist = 0; byDist < 2; ++byDist) {
            for (int64_t crossBlockCorrect = 0; crossBlockCorrect < 2; ++crossBlockCorrect) {
                bool doByDist = byDist;
                bool doCrossBlockCorrect = crossBlockCorrect;
                
                int64_t length;
                double correctness = phasingCorrectness(variants[0], variants[1], decayValues[i],
                                                        doByDist, doCrossBlockCorrect, &length);
                double direct = directLPC(variants[0], variants[1], decayValues[i], doByDist,
                                          doCrossBlockCorrect);
                st_logDebug("decay %f, by dist %d, cross block correct %d, algorithm %f, direct %f\n", decayValues[i], doByDist, doCrossBlockCorrect, correctness, direct);
                CuAssertDblEquals(testCase, correctness, direct, 0.000001);
                // everything in this test is perfectly phased
                CuAssertDblEquals(testCase, correctness, 1.0, 0.000001);
                
                if (doByDist) {
                    // the distance between
                    CuAssertTrue(testCase, length == (numSites - 1) * (numSites - 1));
                }
                else {
                    CuAssertTrue(testCase, length == numSites);
                }
            }
        }
    }
    
    PhasedVariant* var1 = stList_get(variants[1], 1);
    PhasedVariant* var2 = stList_get(variants[1], 6);
    var1->phaseSet[2] = '2'; // "ps1" -> "ps2"
    var2->phaseSet[2] = '1'; // "ps2" -> "ps1"
    
    st_logDebug("TEST BLOCK: imperfect phasing across interlocking sets\n");
    for (int64_t i = 0; i < numDecayValues; ++i) {
        for (int64_t byDist = 0; byDist < 2; ++byDist) {
            for (int64_t crossBlockCorrect = 0; crossBlockCorrect < 2; ++crossBlockCorrect) {
                bool doByDist = byDist;
                bool doCrossBlockCorrect = crossBlockCorrect;
                
                int64_t length;
                double correctness = phasingCorrectness(variants[0], variants[1], decayValues[i],
                                                        doByDist, doCrossBlockCorrect, &length);
                double direct = directLPC(variants[0], variants[1], decayValues[i], doByDist,
                                          doCrossBlockCorrect);
                st_logDebug("decay %f, by dist %d, cross block correct %d, algorithm %f, direct %f\n", decayValues[i], doByDist, doCrossBlockCorrect, correctness, direct);
                CuAssertDblEquals(testCase, correctness, direct, 0.000001);
                // everything in this test is perfectly phased
                CuAssertTrue(testCase, correctness >= 0.0);
                CuAssertTrue(testCase, correctness < 1.00000001); // a little room for numerical slop
                
                if (doByDist) {
                    // the distance between
                    CuAssertTrue(testCase, length == (numSites - 1) * (numSites - 1));
                }
                else {
                    CuAssertTrue(testCase, length == numSites);
                }
            }
        }
    }
    
    // clean up
    for (int64_t i = 0; i < 2; ++i) {
        stList_destruct(variants[i]);
        stList_destruct(alleles[i]);
    }
}

int64_t lpcIntegrationTest(char *truthVcfFile, char *queryVcfFile) {
    // Run localPhasingCorrectness
    char *command = stString_print("./calcLocalPhasingCorrectness %s %s > /dev/null", truthVcfFile, queryVcfFile);
    st_logInfo("> Running command: %s\n", command);
    
    int64_t i = st_system(command);
    free(command);
    return i;
}

void test_executableExecutes(CuTest *testCase) {
    int64_t ret = lpcIntegrationTest(TEST_VCF, TEST_VCF);
    CuAssertTrue(testCase, ret == 0);
}

CuSuite *lpcTestSuite(void) {
    
    //st_setLogLevel(debug);
    
    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_executableExecutes);
    SUITE_ADD_TEST(suite, test_correctValueSimple);
    SUITE_ADD_TEST(suite, test_correctValueWithPhaseSets);

    return suite;
}
