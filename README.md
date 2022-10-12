# PiDecimals

## Introduction

This program computes the Pi number using different Spigot algorithms. 
It can be used as a benchmark to compare and test the CPU performance in a single thread way or even in a multiple threads way.
To perform the operations it is used a floating point precision library. This code allows you to use GMP (The GNU Multiple Precision Arithmetic Library) or MPFR (The GNU Multiple Precision Floating-Point Reliable Library) and comnbine this libraries with the algorithms supported.

### Spigot Algorithms

Currently, PiDecimals allows you to compute Pi using three different algorithms.

* Bailey-Borwein-Plouffe

* Bellard

* Chudnovsky

### Multiple Precision Floating Point Libraries

Currently, PiDecimals allows yoy to compute Pi using two different 

* GMP 

* MPFR 

## Compilation and Installation

To compile the code succesfully it is necessary to have installed OpenMP, GMP and MPFR library. 
If you are using a Linux distro it is very likely that you already have these dependencies installed.

To compile the source code use the "compile.sh" located at the root project directory. 
In the future it is expected to replace the bash script with a make file.   

## Launch

When the source code is compiled you are ready to launch. 

```console
./PiDecimals.out library algorithm precision num_threads [-csv]
```

* Library can be 'GMP' of 'MPFR'
* algorithm is a number between 0-5 if you choose GMP or is a number between 0-5 if you chose de MPFR library
* precision is the precision you want to use to perform the operations. 
* num_threads are the number of threads that you want to use to perform the operations.

