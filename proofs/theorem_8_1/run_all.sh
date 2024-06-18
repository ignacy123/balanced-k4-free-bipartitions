#!/bin/bash
PROBLEM_FILE="problem.txt"

if [ ! -f flag ]; then
    g++ flag_small.cpp -o flag
fi

cat objective.txt > $PROBLEM_FILE
echo "1 1 0" >> $PROBLEM_FILE

echo "" >> $PROBLEM_FILE
cat "assump_20.txt" >> $PROBLEM_FILE
echo "" >> $PROBLEM_FILE
cat "assump_22.txt" >> $PROBLEM_FILE
echo "" >> $PROBLEM_FILE
cat "assump_24.txt" >> $PROBLEM_FILE

./flag -n 7 -ub -fp -obj $PROBLEM_FILE

sage round.sage
