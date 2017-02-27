#!/bin/bash
set -x
mpirun -np 2 -perhost 1 $@
