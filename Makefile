all:
	mpicxx graph.cpp verifier.cpp rds.cpp main.cpp -O3 -lm -o rds -Wall -Wextra -std=c++11 -fopenmp
