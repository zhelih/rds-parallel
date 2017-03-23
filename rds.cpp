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

#define READYTAG  1
#define WORKTAG   2
#define DONETAG   3
#define EXITTAG   4
#define UPDATETAG 5

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

uint find_max(vector<vector <uint> >& c, vector<uint>& weight_c, vector<uint>& p, uint& weight_p, uint* mu, verifier *v, graph* g, vector<uint>& res, int level, const uint time_lim)
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
    if(iter % 1000 == 0) //&& time_lim > 0)
    {
      int flag = false; MPI_Status st;
      MPI_Iprobe(0, UPDATETAG, MPI_COMM_WORLD, &flag, &st); // TODO EXITTAG
      if(flag) // there is a message from master // FIXME copypaste
      {
        int buf[3];
        MPI_Recv(&buf, 3, MPI_INT, 0, UPDATETAG, MPI_COMM_WORLD, &st);
        int cur_i = buf[0];
        for(int i = cur_i; i >= 0; --i)
          mu[i] = max(mu[i], buf[1]);
        lb = max(lb, buf[2]);
        printf("Slave %%d : Received update from master, mu[%d] = %d, lb = %d\n", cur_i, mu[cur_i], lb);
      }
/*      chrono::duration<double> d = chrono::steady_clock::now() - start;
      if(d.count() >= (double)time_lim)
      {
        should_exit = true;
        return lb;
      }*/
      ;
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
    lb = find_max(c, weight_c, p, weight_p, mu, v, g, res, level+1, time_lim);
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
  should_exit = false;
  uint n = g->nr_nodes;
  // order V
  lb = 0; // best solution size found so far
  lb_a = 0;
  uint *mu = new uint[n];
  for(int i = 0; i < n; ++i)
    mu[i] = 1;

  if(world_rank == 0) // master
  {
    chrono::time_point<chrono::steady_clock> start = chrono::steady_clock::now(); // C++11 only
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
      chrono::duration<double> d = chrono::steady_clock::now() - start;
      printf("Master : time elapsed %.8lf secs\n", d.count());
      lb = max(lb, buf[2]);
      lb_a = lb;
      mu[buf[0]] = buf[1];
      if(cur_node >= 0)
      {
        MPI_Send(&cur_node, 1, MPI_INT, st.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD); // send new job
        //TODO Send slaves new updates
        for(int i = 1; i < world_size; ++i)
        {
//          MPI_Request req;
          MPI_Send(&buf, 3, MPI_INT, i, UPDATETAG, MPI_COMM_WORLD);
          //MPI_Request_free(req); FIXME memory leak
        }
        cur_node--;
      }
    }
   for(int i = 1; i < world_size; ++i)
      MPI_Send(NULL, 0, MPI_BYTE, i, EXITTAG, MPI_COMM_WORLD);
    printf("Master done, lb = %d, mu[0] = %d\n", lb, mu[0]); // TODO when CTRL-C
    chrono::duration<double> d = chrono::steady_clock::now() - start;
    printf("rds: time elapsed = %.8lf secs\n", d.count());
  } else { // slave
    MPI_Send(NULL, 0, MPI_BYTE, 0, READYTAG, MPI_COMM_WORLD); // report master ready
    while(!should_exit)
    {
      int cur_i;
      int buf[3];
      MPI_Status st;
      MPI_Recv(&buf, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
      vector<vector<uint> > c; vector<uint> weight_c(g->nr_nodes);
      vector<uint> p; uint weight_p = 0;
      switch(st.MPI_TAG)
      {
        case WORKTAG:
          cur_i = buf[0];
          printf("Slave %d : received work %d\n", world_rank, cur_i);
          // form candidate set
          // take vertices from v \in {i+1, n} for which pair (i,v) satisfies \Pi
          // first iteration c is empty, that must set bound to 1
          c.resize(g->nr_nodes);
          for(uint j = 0; j < g->nr_nodes; ++j)
          {
            c[j].reserve(g->nr_nodes);
            c[j].resize(0);
            weight_c[j] = 0;
          }
          for(uint j = cur_i+1; j < n; ++j)
          {
            if(v->check_pair(g, cur_i, j))
            {
              // add to C
              c[0].push_back(j);
              weight_c[0] += g->weight(j);
            }
          }
          reverse(c[0].begin(), c[0].end()); // for efficient deletion of min element
          p.push_back(cur_i); weight_p += g->weight(cur_i);
          v->init_aux(g, cur_i, c[0]);
          printf("i = %u, c.size = %lu\n", cur_i, c[0].size());
          mu[cur_i] = find_max(c, weight_c, p, weight_p, mu, v, g, res, 0, time_lim);
          printf("mu[%d] = %d\n", cur_i, mu[cur_i]);
          v->free_aux();
          int buf[3];
          buf[0] = cur_i;
          buf[1] = mu[cur_i];
          buf[2] = lb;
          MPI_Send(&buf, 3, MPI_INT, 0, DONETAG, MPI_COMM_WORLD);
          break;
        case UPDATETAG:
          int cur_i = buf[0];
          for(int i = cur_i; i >= 0; --i)
            mu[i] = max(mu[i], buf[1]);
          lb = max(lb, buf[2]);
          printf("Slave %%d : Received update from master, mu[%d] = %d, lb = %d\n", cur_i, mu[cur_i], lb);
          break;
        case EXITTAG:
          should_exit = true;
          printf("Slave %d : exiting...\n", world_rank);
          break;
        default: // FIXME UPDATETAG ?
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
  return fres;
}
