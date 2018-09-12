#ifndef _RDS_H
#define _RDS_H

#include <vector>
#include "graph.h"
#include "verifiers/verifier.hpp"

typedef unsigned int uint;

//Russian-Doll-Search
// input: verifier v, graph g, vector for storing the result res
// time limit in seconds
// output: maximum subgraph size in g satisfying v
// res contains the solution, use graph::restore_order to get original nodes
uint rds(verifier* v, graph*g, std::vector<uint>& res, uint time_lim = 0, bool slave_out=false);

void print_lb_atomic(int signal);

#endif
