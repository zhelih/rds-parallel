#!/bin/bash
set -x
mpirun -np 8 -perhost 1 -genv OMP_NUM_THREADS=16 $@
