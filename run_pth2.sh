#!/bin/bash

time ./susan_serial large2.pgm serial2.pgm
time ./susan_pth 2 large2.pgm omp2_2.pgm
time ./susan_pth 4 large2.pgm omp2_4.pgm
time ./susan_pth 8 large2.pgm omp2_8.pgm
time ./susan_pth 16 large2.pgm omp2_16.pgm

