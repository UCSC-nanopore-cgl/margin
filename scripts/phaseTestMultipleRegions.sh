options=${1}

for i in LOW-chr2 AVG-chr20 AVG-chr7 AVG-chr8 HIGH-chr12; do
  time ../scripts/phaseTest.sh ${i} ${options}
done
