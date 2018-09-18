all:
	mpicxx graph.cpp rds.cpp main.cpp -O3 -lm -o rds -Wall -Wextra -std=c++11 -fopenmp
all-nw:
	mpicxx graph.cpp rds.cpp main.cpp -O3 -lm -o rds -Wall -Wextra -std=c++11 -fopenmp -Wno-unused-parameter
