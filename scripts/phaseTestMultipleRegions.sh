pathToData=$1
company=$2
seq=$3
params=$4
makePlots=$5

for region in LOW-chr2 AVG-chr20 AVG-chr7 AVG-chr8 HIGH-chr12; do
  time ../scripts/phaseTest.sh ${pathToData} ${region} ${company} ${seq} $params $makePlots >&out_${company}_${seq}_${params}_${region}_phased.txt
done
