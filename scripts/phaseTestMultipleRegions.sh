pathToData=$1
company=$2
seq=$3
params=$4
reverseComplement=$5
makePlots=$6

for region in LOW-chr2 AVG-chr20 AVG-chr7 AVG-chr8; do #HIGH-chr12 ## Removing HIGH region cos broken for hifi
  time ../scripts/phaseTest.sh ${pathToData} ${region} ${company} ${seq} $params $reverseComplement $makePlots >&out_${company}_${seq}_${params}_${region}_phased.txt &
done

wait
trap 'kill $(jobs -p)' EXIT
