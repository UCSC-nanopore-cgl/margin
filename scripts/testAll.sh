makePlots=$1

pathToData=/Users/benedictpaten/CLionProjects/MarginPolish/tests/data/externalData

# Phased
time ../scripts/phaseTestMultipleRegions.sh ${pathToData} ont r10.3 allParams.np.human.r103-g3210.json $makePlots
time ../scripts/phaseTestMultipleRegions.sh ${pathToData} pacbio hifi allParams.hifi.json $makePlots
time ../scripts/phaseTestMultipleRegions.sh ${pathToData} ont r9.4 allParams.np.human.r94-g344.json $makePlots

# Unphased
time ../scripts/unphaseTestMultipleRegions.sh ${pathToData} ont r9.4 allParams.np.microbial.r94-g344.json
