#!/bin/bash

time ./susan_serial large3.pgm serial3.pgm
time ./susan_pth 2 large3.pgm omp3_2.pgm
time ./susan_pth 4 large3.pgm omp3_4.pgm
time ./susan_pth 8 large3.pgm omp3_8.pgm
time ./susan_pth 16 large3.pgm omp3_16.pgm

