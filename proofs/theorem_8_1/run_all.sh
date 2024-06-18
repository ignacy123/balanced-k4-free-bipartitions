#!/bin/bash
PROBLEM_FILE="problem.txt"

if [ ! -f flag ]; then
    g++ flag_small.cpp -o flag
fi

cat objective.txt > $PROBLEM_FILE
echo "1 1 0" >> $PROBLEM_FILE

add_assumption () {
  echo "
# assumption $1 from paper
0 0 " >> $PROBLEM_FILE
  cat "assump_$1.txt" >> $PROBLEM_FILE
}

add_assumption "20"
add_assumption "22"
add_assumption "24"

./flag -n 7 -ub -fp -obj $PROBLEM_FILE
