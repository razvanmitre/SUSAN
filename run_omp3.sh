#!/bin/bash

time ./susan_serial large3.pgm serial3.pgm
time ./susan_omp 2 large3.pgm omp3_2.pgm
time ./susan_omp 4 large3.pgm omp3_4.pgm
time ./susan_omp 8 large3.pgm omp3_8.pgm
time ./susan_omp 16 large3.pgm omp3_16.pgm

