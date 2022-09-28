#!/bin/bash

./compile.sh

./PiDecimals.x GMP 0 30000 1
./PiDecimals.x GMP 1 30000 1
./PiDecimals.x GMP 2 30000 1
./PiDecimals.x GMP 3 30000 1
./PiDecimals.x GMP 4 30000 1
./PiDecimals.x GMP 5 30000 1

./PiDecimals.x MPFR 0 30000 1
./PiDecimals.x MPFR 1 30000 1
./PiDecimals.x MPFR 2 30000 1

./PiDecimals.x GMP 0 30000 16
./PiDecimals.x GMP 1 30000 16
./PiDecimals.x GMP 2 30000 16
./PiDecimals.x GMP 3 30000 16
./PiDecimals.x GMP 4 30000 16
./PiDecimals.x GMP 5 30000 16

./PiDecimals.x MPFR 0 30000 16
./PiDecimals.x MPFR 1 30000 16
./PiDecimals.x MPFR 2 30000 16


