#!/bin/csh

## Compilation 
ifort  -c fft.f bb86.f 
ifort  -o bb86 -i8 -r8  fft.o bb86.o

\rm *.o
echo 'Compilation done !'

