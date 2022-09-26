#!/bin/bash

program="$1"
error=""
RED_OUTPUT="tput setaf 1"
GREEN_OUTPUT="tput setaf 2"
RESET_OUTPUT="tput sgr0"

errors(){
    echo "No params are needed"
    exit 1
}

#CHECK PARAMS
if [ "$#" -ne 0 ]; then
   errors
fi

error=$(gcc -fopenmp -o PiDecimals.x Sources/Common/*.c Sources/GMP/*.c Sources/GMP/Algorithms/*.c Sources/MPFR/*.c Sources/MPFR/Algorithms/*.c -lmpfr -lgmp -lm 2>&1 1>/dev/null)


#GIVE FEEDBACK ABOUT COMPILATION
if [[ -z "$error" ]]; then
    echo -n "COMPILATION "
    ${GREEN_OUTPUT}
    echo -n "DONE"
    ${RESET_OUTPUT}
    echo " SUCCESFULLY "

else
    ${RED_OUTPUT}
    echo -n "ERROR: "
    ${RESET_OUTPUT}
    echo "$error"
fi

