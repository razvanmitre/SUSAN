#!/bin/bash

time ./susan_serial large.pgm serial.pgm
time ./susan_omp 2 large.pgm omp_2.pgm
time ./susan_omp 4 large.pgm omp_4.pgm
time ./susan_omp 8 large.pgm omp_8.pgm
time ./susan_omp 16 large.pgm omp_16.pgm

