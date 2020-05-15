#!/bin/bash

executable="../src/solver"

if [ ! -x "$executable" ]; then
    echo "Error: expecting executable -> $executable"
    exit 1
fi

err=0.000001

$executable ./input/input_12g.dat
./regression output.dat ./ref/ref_12g.dat > result.dat
res="$(cat result.dat)"
if [ $(echo "$res <= $err" | bc -l) = 1 ]; then
	echo "the difference between analytical and numerical solutions are too large"
        exit 1
fi

$executable ./input/input_12j.dat
./regression output.dat ./ref/ref_12j.dat > result.dat
res="$(cat result.dat)"
if [ $(echo "$res <= $err" | bc -l) = 1 ]; then
        echo "the difference between analytical and numerical solutions are too large"
	exit 1
fi

$executable ./input/input_14g.dat
./regression output.dat ./ref/ref_14g.dat > result.dat
res="$(cat result.dat)"
if [ $(echo "$res <= $err" | bc -l) = 1 ]; then
        echo "the difference between analytical and numerical solutions are too large"
        exit 1
fi

$executable ./input/input_22g.dat
./regression output.dat ./ref/ref_22g.dat > result.dat
res="$(cat result.dat)"
if [ $(echo "$res <= $err" | bc -l) = 1 ]; then
        echo "the difference between analytical and numerical solutions are too large"
        exit 1
fi

$executable ./input/input_22j.dat
./regression output.dat ./ref/ref_22j.dat > result.dat
res="$(cat result.dat)"
if [ $(echo "$res <= $err" | bc -l) = 1 ]; then
        echo "the difference between analytical and numerical solutions are too large"
        exit 1
fi

$executable ./input/input_24g.dat
./regression output.dat ./ref/ref_24g.dat > result.dat
res="$(cat result.dat)"
if [ $(echo "$res <= $err" | bc -l) = 1 ]; then
        echo "the difference between analytical and numerical solutions are too large"
        exit 1
fi
