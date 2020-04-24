pathToData=$1
company=$2
seq=$3
params=$4

for region in MICROBE_100k1 MICROBE_100k2; do
  for coverage in 32 64 128; do
    time ../scripts/unphaseTest.sh ${pathToData} ${region} ${company} ${seq} $params ${coverage} >& out_${company}_${seq}_${params}_${region}_${coverage}_unphased.txt
  done
done
