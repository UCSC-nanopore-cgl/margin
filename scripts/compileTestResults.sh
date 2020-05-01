#!/bin/zsh

## Input is path containing output files
tempFile=temp.txt # Temporary file used

function get_qv {
  echo `grep total-xLen $i | head -n4 | cut -f${j} -d' ' | sort | tail -n2 | cut -f1 -d',' | python3 -c "import sys; l=[float(l) for l in sys.stdin]; print(sum(l)/len(l))"` ${i} >> ${tempFile};
}

function print_qvs {
  cat ${tempFile} | cut -f1 -d' ' | python3 -c "import sys; l=[float(l) for l in sys.stdin]; print(sum(l)/len(l))"
  cat ${tempFile}
  rm -f ${tempFile}
}

rm -f ${tempFile}

j=10; echo "Average qvs"; for i in ${1}/out*.txt; do; get_qv; done;
print_qvs

j=12; echo "Average qv-matches"; for i in ${1}/out*.txt; do; get_qv; done;
print_qvs

j=14; echo "Average qv-inserts"; for i in ${1}/out*.txt; do; get_qv; done;
print_qvs

j=16; echo "Average qv-deletes"; for i in ${1}/out*.txt; do; get_qv; done;
print_qvs

