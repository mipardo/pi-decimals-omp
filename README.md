<img src="https://user-images.githubusercontent.com/60443339/195345025-2f54a691-7918-4a3a-89de-2d488625f4b1.png" alt="drawing"/>

This program computes the Pi number using different Spigot algorithms. 
To perform the operations it is used a floating point precision arithmetic library. PiDecimals allows you to use GMP (The GNU Multiple Precision Arithmetic Library) or MPFR (The GNU Multiple Precision Floating-Point Reliable Library) to compute the Pi number with the algorithms supported.

It can be used as a benchmark to compare and test the CPU performance in a single thread way or even with multiple threads. Threads are managed with OpenMP.

<img src="https://user-images.githubusercontent.com/60443339/195342306-1eb14b7d-ce25-41a8-87b5-545011edf172.png" alt="drawing" height="100"/>

#### Spigot Algorithms

Currently, PiDecimals allows you to compute Pi using three different algorithms:

* Bailey-Borwein-Plouffe. The expression is presented below:

<img src="https://user-images.githubusercontent.com/60443339/195336253-bf6aeeea-c255-458c-9f16-7fcc91d5b2c7.png" alt="drawing" height="85" />

$$64\pi = \sum_{n=0}^{\infty}\frac{-1^{n}}{1024^{n}} (-\frac{32}{4n + 1} - \frac{1}{4n + 3} + \frac{256}{10n + 1} - \frac{64}{10n + 3} - \frac{4}{10n + 5} - \frac{4}{10n + 7} + \frac{1}{10n + 9})
    \label{equation:Bellard}
$$

* Bellard. The expression is presented as follows:

<img src="https://user-images.githubusercontent.com/60443339/195340916-7508ee10-2209-413a-b24a-92cede2aea44.png" alt="drawing" />

* Chudnovsky. The expression is shown below:

<img src="https://user-images.githubusercontent.com/60443339/195336414-27422fd3-4884-4cf4-a7b8-47bf49f5b67a.png" alt="drawing" height="85" />

#### Multiple Precision Floating Point Libraries

Currently, PiDecimals allows yoy to compute Pi using two different floating point arithmetic libraries: 

* GMP (https://gmplib.org/)

* MPFR (https://www.mpfr.org/)

## Compilation and Installation

To compile the code succesfully it is necessary to have installed OpenMP, GMP and MPFR library. 
If you are using a Linux distro it is very likely that you already have these dependencies installed.

To compile the source code use the "compile.sh" located at the root project directory. 
In the future it is expected to replace the bash script with a Make file.   

## Launch

When the source code is compiled you are ready to launch: 

```console
./PiDecimals.x library algorithm precision num_threads [-csv]
```

* library can be 'GMP' or 'MPFR'
* algorithm is a value between 0 and X. The X value may depend on the library used.
* precision param is the value of precision you want to use to perform the operations. 
* num_threads param is the number of threads that you want to use to perform the operations.
* -csv param is optional. If this param is used the program will show the results in csv format.

En example of use could be:
```console
./PiDecimals.x MPFR 1 50000 4 
```
And the output could be:

<img src="https://user-images.githubusercontent.com/60443339/195339706-35d6c8d7-a0e1-4732-87dd-57f886090160.png" alt="drawing" />



