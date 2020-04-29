makePlots=$1

pathToData=/Users/benedictpaten/data/margin/externalData

# Phased
time ../scripts/phaseTestMultipleRegions.sh ${pathToData} ont r10.3 allParams.np.human.r103-g3210.json TRUE $makePlots &
time ../scripts/phaseTestMultipleRegions.sh ${pathToData} pacbio hifi allParams.hifi.json TRUE $makePlots &
time ../scripts/phaseTestMultipleRegions.sh ${pathToData} ont r9.4 allParams.np.human.r94-g344.json FALSE $makePlots &

# Unphased
time ../scripts/unphaseTestMultipleRegions.sh ${pathToData} ont r9.4 allParams.np.microbial.r94-g344.json &

# Wait until these things are finished
wait

# Cleanup the mess
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT
