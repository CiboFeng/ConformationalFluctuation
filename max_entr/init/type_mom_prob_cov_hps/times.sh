#!/bin/bash

touch times.dat
echo "init:" >> times.dat

lines=$(grep 'init' bash.log)

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $3}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | awk 'NR%3==1' >> times.dat
    fi
done

echo " " >> times.dat
echo "file:" >> times.dat

lines=$(grep 'file' bash.log)

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $3}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | awk 'NR%3==1' >> times.dat
    fi
done

echo " " >> times.dat
echo "sim:" >> times.dat

lines=$(grep 'sim' bash.log)

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $3}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | awk 'NR%4==1' >> times.dat
    fi
done

echo " " >> times.dat
echo "pc:" >> times.dat

lines=$(grep 'pc' bash.log)

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $3}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | awk 'NR%3==1' >> times.dat
    fi
done

echo " " >> times.dat
echo "coeff:" >> times.dat

lines=$(grep 'coeff' bash.log)

numbers=()
while IFS= read -r line; do
    number=$(echo "$line" | awk '{print $3}')
    if [[ $number =~ ^[0-9]+$ ]]; then
        numbers+=($number)
    fi
done <<< "$lines"

for x in "${numbers[@]}"; do
    result=$(sacct -j $x -n -o Elapsed 2>/dev/null)
    if [ -n "$result" ]; then
        echo "$x: $result" | awk 'NR%3==1' >> times.dat
    fi
done

echo " " >> times.dat
echo "---------------------------------------------" >> times.dat
echo " " >> times.dat

