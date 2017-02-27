#include "rds.h"
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <atomic>
#include <chrono>
#include <mpi.h>

#include <csignal> // Display best result for SIGINT before exit

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

using namespace std;

typedef unsigned int uint;

static uint lb;

atomic_uint lb_a;

#define READYTAG 1
#define WORKTAG  2
#define DONETAG  3
#define EXITTAG  4

// for debug
void print_cont(const vector<uint>& c)
{
  for(auto it = c.begin(); it != c.end(); ++it)
    printf("%u ", *it);
  printf("\n");
}

void print_lb_atomic(int signal)
{
  printf("\nReceived SIGINT\n");
  printf("Best lower bound found: %u\n", lb_a.load());
  exit(0);
}

static uint iter = 0;
static bool should_exit = false;

uint find_max(vector<vector <uint> >& c, vector<uint>& weight_c, vector<uint>& p, uint& weight_p, const uint* mu, verifier *v, graph* g, vector<uint>& res, int level, const chrono::time_point<chrono::steady_clock> start, const uint time_lim)
{
  if(should_exit)
    return lb;
  if(c[level].size() == 0)
  {
    if(weight_p > lb)
    {
      res = p; //copy
      return weight_p;
    }
    else
      return lb;
  }

  while(c[level].size() > 0)
  {
    iter++;
    if(iter % 1000 == 0 && time_lim > 0)
    {
      chrono::duration<double> d = chrono::steady_clock::now() - start;
      if(d.count() >= (double)time_lim)
      {
        should_exit = true;
        return lb;
      }
    }

    if(weight_c[level] + weight_p <= lb) // Prune 1
      return lb;
    uint i = c[level][c[level].size()-1];
    if(mu[i] + weight_p <= lb) // Prune 2
      return lb;
    c[level].pop_back(); weight_c[level] -= g->weight(i);
//    NB: exploit that we adding only 1 vertex to p
//    thus verifier can prepare some info using prev calculations
    v->prepare_aux(g, p, i, c[level]);
    p.push_back(i); weight_p += g->weight(i);
    c[level+1].resize(0); weight_c[level+1] = 0;
    for(uint it2 = 0; it2 < c[level].size(); ++it2)
    {
      if(c[level][it2] != i && v->check(g, p, c[level][it2]))
      {
        c[level+1].push_back(c[level][it2]);
        weight_c[level+1] += g->weight(c[level][it2]);
      }
    }
    lb = find_max(c, weight_c, p, weight_p, mu, v, g, res, level+1, start, time_lim);
    lb_a = lb;
    p.pop_back(); weight_p -= g->weight(i);
    v->undo_aux(g, p, i, c[level]);
  }
  return lb;
}

uint rds(verifier* v, graph* g, vector<uint>& res, uint time_lim)
{
  MPI_Init(NULL, NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  chrono::time_point<chrono::steady_clock> start = chrono::steady_clock::now(); // C++11 only
  should_exit = false;
  uint n = g->nr_nodes;
  // order V
  lb = 0; // best solution size found so far
  lb_a = 0;
  uint *mu = new uint[n];

  if(world_rank == 0) // master
  {
    // synchronize here
    for(int i = 1; i < world_size; ++i)
      MPI_Recv(NULL, 0, MPI_BYTE, i, READYTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Master: all slaves ready!\n");
    int cur_node = n-1;
    int rec_left = cur_node;
    for(int i = 1; cur_node >= 0 && i < world_size; ++i)
    {
      MPI_Send(&cur_node, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
      cur_node--;
    } // tasks send
    while(cur_node >= 0 || rec_left >= 0)
    {
      int buf[3];
      MPI_Status st;
      MPI_Recv(&buf, 3, MPI_INT, MPI_ANY_SOURCE, DONETAG, MPI_COMM_WORLD, &st);
      rec_left--;
      printf("Master : received result from slave %d : mu[%d] = %d, lb = %d\n", st.MPI_SOURCE, buf[0], buf[1], buf[2]);
      //TODO update mu[i] and lb here
      if(cur_node >= 0)
      {
        MPI_Send(&cur_node, 1, MPI_INT, st.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD); // send new job
        //TODO Send slaves new updates
        cur_node--;
      }
    }
   for(int i = 1; i < world_size; ++i)
      MPI_Send(NULL, 0, MPI_BYTE, i, EXITTAG, MPI_COMM_WORLD);
    printf("Master done, lb = %d, mu[0] = %d\n", 0, 0); // TODO out

  } else { // slave
    MPI_Send(NULL, 0, MPI_BYTE, 0, READYTAG, MPI_COMM_WORLD); // report master ready
    while(!should_exit)
    {
      int cur_i;
      MPI_Status st;
      MPI_Recv(&cur_i, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
      switch(st.MPI_TAG)
      {
        case WORKTAG:
          printf("Slave %d : received work %d\n", world_rank, cur_i);
          //TODO start working
          int buf[3]; // TODO assing the result
          buf[0] = cur_i;
          MPI_Send(&buf, 3, MPI_INT, 0, DONETAG, MPI_COMM_WORLD);
          break;
        case EXITTAG:
          should_exit = true;
          printf("Slave %d : exiting...\n", world_rank);
          break;
        default:
          printf("Slave %d : wrong message tag %d, ignoring...\n", world_rank, st.MPI_TAG);
      }
    }
  }
/*
  int i;
  for(i = n-1; i >= 0; --i)
  {
    // form candidate set
    // take vertices from v \in {i+1, n} for which pair (i,v) satisfies \Pi
    // first iteration c is empty, that must set bound to 1
    vector<vector<uint> > c; vector<uint> weight_c(g->nr_nodes);
    c.resize(g->nr_nodes);
    for(uint j = 0; j < g->nr_nodes; ++j)
    {
      c[j].reserve(g->nr_nodes);
      c[j].resize(0);
      weight_c[j] = 0;
    }
    for(uint j = i+1; j < n; ++j)
    {
      if(v->check_pair(g, i, j))
      {
        // add to C
        c[0].push_back(j);
        weight_c[0] += g->weight(j);
      }
    }
    reverse(c[0].begin(), c[0].end()); // for efficient deletion of min element
    vector<uint> p; uint weight_p = 0;
    p.push_back(i); weight_p += g->weight(i);
    v->init_aux(g, i, c[0]);
    printf("i = %u, c.size = %lu, ", i, c[0].size());
    mu[i] = find_max(c, weight_c, p, weight_p, mu, v, g, res, 0, start, time_lim);
    printf("mu[%d] = %d\n", i, mu[i]);
    v->free_aux();
    if(time_lim > 0)
    {
      chrono::duration<double> d = chrono::steady_clock::now() - start;
      if((uint)(d.count()) >= time_lim)
        break;
    }
  }
  printf("RDS done\n");
  uint fres = mu[i+1]; // last
  delete [] mu;*/
  uint fres = 0;
  chrono::duration<double> d = chrono::steady_clock::now() - start;
  printf("rds: time elapsed = %.8lf secs\n", d.count());
  return fres;
}
