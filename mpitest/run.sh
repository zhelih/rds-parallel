#!/bin/bash
set -x
mpirun -np 8 OMP_NUM_THREADS=16 $@
