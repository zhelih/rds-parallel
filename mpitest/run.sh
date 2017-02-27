#!/bin/bash
set -x
mpirun -np 8 -perhost 1 $@
