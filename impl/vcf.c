//
// Created by tpesout on 5/28/20.
//

#include "margin.h"

VcfEntry *vcfEntry_construct(char *refSeqName, int64_t refPos, RleString *allele1, RleString *allele2) {
    VcfEntry *vcfEntry = st_calloc(1, sizeof(VcfEntry));
    vcfEntry->refSeqName = stString_copy(refSeqName);
    vcfEntry->refPos = refPos;
    vcfEntry->allele1 = allele1;
    vcfEntry->allele2 = allele2;
}

void vcfEntry_destruct(VcfEntry *vcfEntry) {
    if (vcfEntry->allele1 != NULL) rleString_destruct(vcfEntry->allele1);
    if (vcfEntry->allele2 != NULL) rleString_destruct(vcfEntry->allele2);
    free(vcfEntry->refSeqName);
    free(vcfEntry);
}

stList *parseVcf2(char *vcfFile, bool hetOnly, PolishParams *params) {
    st_logCritical("Parsing VCF for %s entries: %s\n", hetOnly ? "HET" : "all", vcfFile);
    stList *entries = stList_construct3(0, (void(*)(void*))vcfEntry_destruct);
    FILE *fp = fopen(vcfFile, "r");
    if (fp == NULL) {
        st_errAbort("Could not open VCF %s\n", vcfFile);
    }

    char *line = NULL;
    while ((line = stFile_getLineFromFile(fp)) != NULL) {
        if (strlen(line) == 0 || line[0] == '#') {
            free(line);
            continue;
        }
        stList *elements = stString_split(line);
        assert(stList_length(elements) >= 10);

        // get genotype idx
        stList *format = stString_splitByString(stList_get(elements, 8), ":");
        int64_t gtIdx = -1;
        for (int64_t i = 0; i < stList_length(format); i++) {
            if (stString_eq("GT", stList_get(format, i))) {
                gtIdx = i;
                break;
            }
        }
        if (gtIdx == -1) {
            st_logInfo("  Missing GT record for format string '%s' in VCF line:\n\t\t%s\n", stList_get(elements, 8), line);
            continue;
        }

        // get genotype
        stList *sample = stString_splitByString(stList_get(elements, 9), ":");
        char *genotypeStr = stList_get(sample, gtIdx);
        if (strlen(genotypeStr) != 3) {
            st_logInfo("  Unexpected genotype str '%s' in VCF line:\n\t\t%s\n", genotypeStr, line);
            continue;
        }
        int64_t gt1 = genotypeStr[0] - '0';
        int64_t gt2 = genotypeStr[2] - '0';

        // get alleles
        stList *alleles = stList_construct(); //ref is freed in elements, alts are freed in altAlleles
        stList_append(alleles, stString_copy(stList_get(elements, 3))); //ref
        stList *altAlleles = stString_splitByString(stList_get(elements, 4), ",");
        stList_appendAll(alleles, altAlleles);
        assert(stList_length(alleles) >= (gt1 > gt2 ? gt1 : gt2));
        char *allele1 = stList_get(alleles, gt1);
        char *allele2 = stList_get(alleles, gt2);

        // location data
        char *chrom = stList_get(elements, 0);
        int64_t pos = atol(stList_get(elements, 1));

        // save it
        if (gt1 != gt2 || !hetOnly) {
            RleString *rleAllele1 = params->useRunLengthEncoding ? rleString_construct(allele1)
                                                                 : rleString_construct_no_rle(allele1);
            RleString *rleAllele2 = params->useRunLengthEncoding ? rleString_construct(allele2)
                                                                 : rleString_construct_no_rle(allele2);
            VcfEntry *entry = vcfEntry_construct(chrom, pos, rleAllele1, rleAllele2);
            stList_append(entries, entry);
        }

        // cleanup
        stList_destruct(format);
        stList_destruct(sample);
        stList_destruct(elements);
        free(line);

    }

    // cleanup
    fclose(fp);

    st_logCritical("Read %"PRId64" VCF entries from %s.\n", stList_length(entries), vcfFile);
    return entries;
}
stList *parseVcf(char *vcfFile, PolishParams *params) {
    return parseVcf2(vcfFile, TRUE, params);
}