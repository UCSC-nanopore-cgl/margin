//
// Created by tpesout on 3/29/19.
//

#ifndef MARGINPHASE_HELENFEATURES_H
#define MARGINPHASE_HELENFEATURES_H

#ifdef _HDF5

#include "margin.h"
#include <hdf5.h>

#define SYMBOL_NUMBER 5
#define SYMBOL_NUMBER_NO_N 4
#define POAFEATURE_SYMBOL_GAP_POS SYMBOL_NUMBER_NO_N
#define POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE ((SYMBOL_NUMBER) * 2) // {A,C,G,T,gap} x {fwd,bkwd}

typedef struct _poaFeatureSimpleWeight PoaFeatureSimpleWeight;
struct _poaFeatureSimpleWeight {
    int64_t refPosition;
    int64_t insertPosition;
    char label;
    double weights[POAFEATURE_SIMPLE_WEIGHT_TOTAL_SIZE];
    PoaFeatureSimpleWeight *nextInsert; //so we can model all inserts after a position
};

typedef struct _poaFeatureSplitRleWeight PoaFeatureSplitRleWeight;
struct _poaFeatureSplitRleWeight {
    int64_t refPosition;
    int64_t insertPosition;
    int64_t runLengthPosition;
    char labelChar;
    int64_t labelRunLength;
    PoaFeatureSplitRleWeight *nextRunLength; //so we can model all inserts after a position
    PoaFeatureSplitRleWeight *nextInsert; //so we can model all inserts after a position
    double *weights;
    int64_t maxRunLength;
};

typedef struct _poaFeatureChannelRleWeight PoaFeatureChannelRleWeight;
struct _poaFeatureChannelRleWeight {
    int64_t refPosition;
    int64_t insertPosition;
    int64_t runLengthPosition;
    char labelChar;
    int64_t labelRunLength;
    PoaFeatureChannelRleWeight *nextRunLength; //so we can model all inserts after a position
    PoaFeatureChannelRleWeight *nextInsert; //so we can model all inserts after a position
    double *nucleotideWeights;
    double *runLengthWeights;
    int64_t maxRunLength;
};

typedef struct _poaFeatureDiploidRleWeight PoaFeatureDiploidRleWeight;
struct _poaFeatureDiploidRleWeight {
    int64_t refPosition;
    int64_t insertPosition;
    int64_t runLengthPosition;
    int64_t referencePosition;
    char labelChar;
    int64_t labelRunLength;
    PoaFeatureDiploidRleWeight* nextRunLength;
    PoaFeatureDiploidRleWeight* nextInsert;
    double* weightsHOn;
    double* weightsHOff;
    int64_t maxRunLength;
};

typedef struct _HelenFeatureHDF5FileInfo HelenFeatureHDF5FileInfo;
struct _HelenFeatureHDF5FileInfo {
    char *filename;
    hid_t file;
    hid_t int64Type;
    hid_t uint32Type;
    hid_t uint8Type;
    hid_t floatType;
    hid_t groupPropertyList;
};

HelenFeatureHDF5FileInfo *HelenFeatureHDF5FileInfo_construct(char *filename);

void HelenFeatureHDF5FileInfo_destruct(HelenFeatureHDF5FileInfo *fileInfo);

HelenFeatureHDF5FileInfo **openHelenFeatureHDF5FilesByThreadCount(char *filenameBase, int64_t threadCount);

int PoaFeature_SimpleWeight_charIndex(Symbol character, bool forward);

int PoaFeature_SimpleWeight_gapIndex(bool forward);
int PoaFeature_SplitRleWeight_charIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward);

int PoaFeature_SplitRleWeight_gapIndex(int64_t maxRunLength, bool forward);

int PoaFeature_ChannelRleWeight_charNuclIndex(Symbol character, bool forward);

int PoaFeature_ChannelRleWeight_gapNuclIndex(bool forward);

int PoaFeature_ChannelRleWeight_charRLIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward);
int PoaFeature_DiploidRleWeight_charIndex(int64_t maxRunLength, Symbol character, int64_t runLength, bool forward);
int PoaFeature_DiploidRleWeight_gapIndex(int64_t maxRunLength, bool forward);

PoaFeatureSimpleWeight *PoaFeature_SimpleWeight_construct(int64_t refPos, int64_t insPos);

void PoaFeature_SimpleWeight_destruct(PoaFeatureSimpleWeight *feature);

PoaFeatureSplitRleWeight *
PoaFeature_SplitRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos, int64_t maxRunLength);

void PoaFeature_SplitRleWeight_destruct(PoaFeatureSplitRleWeight *feature);

PoaFeatureChannelRleWeight *
PoaFeature_ChannelRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos, int64_t maxRunLength);

void PoaFeature_ChannelRleWeight_destruct(PoaFeatureChannelRleWeight *feature);

PoaFeatureDiploidRleWeight *PoaFeature_DiploidRleWeight_construct(int64_t refPos, int64_t insPos, int64_t rlPos, int64_t maxRunLength);
void PoaFeature_DiploidRleWeight_destruct(PoaFeatureDiploidRleWeight *feature);

stList *poa_getSimpleWeightFeatures(Poa *poa, stList *bamChunkReads);

stList *poa_getSplitRleWeightFeatures(Poa *poa, stList *bamChunkReads, int64_t maxRunLength);

stList *poa_getChannelRleWeightFeatures(Poa *poa, stList *bamChunkReads, int64_t maxRunLength);
stList *poa_getDiploidRleWeightFeatures(Poa *poa, stList *bamChunkReads, stSet *onHapReads, const int64_t maxRunLength);

void handleHelenFeatures(HelenFeatureType helenFeatureType, BamChunker *trueReferenceBamChunker,
                         int64_t splitWeightMaxRunLength, void **helenHDF5Files, bool fullFeatureOutput,
                         char *trueReferenceBam,
                         Params *params, char *logIdentifier, int64_t chunkIdx, BamChunk *bamChunk, Poa *poa,
                         stList *bamChunkReads,
                         char *polishedConsensusString, RleString *polishedRleConsensus);

void handleDiploidHelenFeatures(HelenFeatureType helenFeatureType, BamChunker *trueReferenceBamChunker,
        int64_t splitWeightMaxRunLength, void **helenHDF5Files, bool fullFeatureOutput,
        char *trueReferenceBamA, char *trueReferenceBamB, Params *params,
        char *logIdentifier, int64_t chunkIdx, BamChunk *bamChunk, stList *bamChunkReads, Poa *poaH1, Poa *poaH2,
        stSet *readsInH1, stSet *readsInH2, RleString *polishedRleConsensusH1, RleString *polishedRleConsensusH2,
        RleString *originalReference);

void poa_writeHelenFeatures(HelenFeatureType type, Poa *poa, stList *bamChunkReads,
                            char *outputFileBase, BamChunk *bamChunk, stList *trueRefAlignment,
                            RleString *consensusRleString,
                            RleString *trueRefRleString, bool fullFeatureOutput, int64_t splitWeightMaxRunLength,
                            HelenFeatureHDF5FileInfo **helenHDF5Files);

void poa_writeDiploidHelenFeatures(HelenFeatureType type, stList *bamChunkReads, char *outputFileBase,
        BamChunk *bamChunk, Poa *poaH1, Poa *poaH2, stSet *readsInH1, stSet *readsInH2,
        stList *trueRefAlignmentToH1, stList *trueRefAlignmentToH2,
        RleString *trueRefRleStringToH1, RleString *trueRefRleStringToH2,
        int64_t maxRunLength, HelenFeatureHDF5FileInfo** helenHDF5Files);

stList *alignConsensusAndTruthSSW(char *consensusStr, char *truthStr, uint16_t *score);
stList *alignConsensusAndTruthCPECAN(char *consensusStr, char *truthStr, double *score, PolishParams *polishParams);
stList *alignConsensusAndTruthRLE(RleString *consensusStr, RleString *truthStr, double *score, PolishParams *polishParams);
stList *alignConsensusAndTruthRLEWithKmerAnchors(RleString *consensusStr, RleString *truthStr, double *score,
                                                 PolishParams *polishParams);
stList *alignConsensusAndTruthRLEWithSSWAnchors(RleString *consensusStr, RleString *truthStr, double *score,
                                                PolishParams *polishParams);
void poa_annotateHelenFeaturesWithTruth(stList *features, HelenFeatureType featureType, stList *trueRefAlignment,
                                        RleString *trueRefRleString, int64_t *firstMatchedFeaure,
                                        int64_t *lastMatchedFeature);

void printMEAAlignment(char *X, char *Y, int64_t lX, int64_t lY, stList *alignedPairs, uint64_t *Xrl, uint64_t *Yrl);
void printMEAAlignment2(RleString *X, RleString *Y, stList *alignedPairs);

void
writeSimpleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                   BamChunk *bamChunk, bool outputLabels, stList *features,
                                   int64_t featureStartIdx, int64_t featureEndIdxInclusive);

void
writeSplitRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                     BamChunk *bamChunk, bool outputLabels, stList *features,
                                     int64_t featureStartIdx, int64_t featureEndIdxInclusive, int64_t maxRunLength);

void
writeChannelRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo *hdf5FileInfo, char *outputFileBase,
                                       BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx,
                                       int64_t featureEndIdxInclusive, const int64_t maxRunLength);

void writeChannelRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo* hdf5FileInfo, char *outputFileBase,
        BamChunk *bamChunk, bool outputLabels, stList *features, int64_t featureStartIdx,
        int64_t featureEndIdxInclusive, int64_t maxRunLength);

void writeDiploidRleWeightHelenFeaturesHDF5(Alphabet *alphabet, HelenFeatureHDF5FileInfo* hdf5FileInfo,
                                            char *outputFileBase, BamChunk *bamChunk, bool outputLabels,
                                            stList *featuresH1, int64_t featureStartIdxH1, int64_t featureEndIdxInclusiveH1,
                                            stList *featuresH2, int64_t featureStartIdxH2, int64_t featureEndIdxInclusiveH2,
                                            int64_t maxRunLength);
#endif //_HDF5
#endif //MARGINPHASE_HELENFEATURES_H
