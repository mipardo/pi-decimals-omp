# PiDecimals

## Introduction

This program computes the Pi number using different Spigot algorithms. 
It can be used as a benchmark to compare and test the CPU performance in a single thread way or even with multiple threads.
To perform the operations it is used a floating point precision arithmetic library. PiDecimals allows you to use GMP (The GNU Multiple Precision Arithmetic Library) or MPFR (The GNU Multiple Precision Floating-Point Reliable Library) and combine this libraries with the algorithms supported.

### Spigot Algorithms

Currently, PiDecimals allows you to compute Pi using three different algorithms.

* Bailey-Borwein-Plouffe

<img src="https://user-images.githubusercontent.com/60443339/195336253-bf6aeeea-c255-458c-9f16-7fcc91d5b2c7.png" alt="drawing" width="200"/>

![image](https://user-images.githubusercontent.com/60443339/195336253-bf6aeeea-c255-458c-9f16-7fcc91d5b2c7.png)

* Bellard

![image](https://user-images.githubusercontent.com/60443339/195336107-7465da26-237c-4a67-8d18-00bc4136e8ca.png)

* Chudnovsky

![image](https://user-images.githubusercontent.com/60443339/195336414-27422fd3-4884-4cf4-a7b8-47bf49f5b67a.png)

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

