#!/bin/sh
NUM_PROCS=1
mpirun -np $NUM_PROCS ../build/prog input.txt