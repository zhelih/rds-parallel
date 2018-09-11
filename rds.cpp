#include "rds.h"
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <omp.h>

#include <mpi.h>

#include <csignal> // Display best result for SIGINT before exit

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

using namespace std;

typedef unsigned int uint;

atomic_uint lb;

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
  printf("Best lower bound found: %u\n", lb.load());
  exit(0);
}

atomic_uint iter (0);
atomic_bool should_exit (false);

// MPI TAGS
#define READYTAG  1
#define WORKTAG   2
#define DONETAG   3
#define EXITTAG   4
#define UPDATETAG 5

void find_max(vector<vector <uint> >& c, vector<uint>& weight_c, vector<uint>& p, uint weight_p, const int* mu, verifier *v, graph* g, vector<uint>& res, int level, const chrono::time_point<chrono::steady_clock> start, const uint time_lim)
{
  if(should_exit)
    return;
  if(c[level].size() == 0)
  {
    if(weight_p > lb)
    {
      #pragma omp critical (lbupdate)
      {
//      res = p; //copy
      lb.store(max(lb.load(), weight_p));
      }
      return;
    }
    else
      return;
  }

  for(uint c_i = 0; c_i < c[level].size(); ++c_i)
  {
    iter++;
    if(iter % 1000 == 0 && time_lim > 0)
    {
      chrono::duration<double> d = chrono::steady_clock::now() - start;
      if(d.count() >= (double)time_lim)
      {
        should_exit = true;
        return;
      }
    }

    if(weight_c[level] + weight_p <= lb) // Prune 1
    {
      return;
    }
    uint i = c[level][c_i];
    if(mu[i] > 0 && mu[i] + weight_p <= lb) // Prune 2
    {
      return;
    }
    weight_c[level] -= g->weight(i);
//    NB: exploit that we adding only 1 vertex to p
//    thus verifier can prepare some info using prev calculations
    v->prepare_aux(p, i, c[level]);
    p.push_back(i); weight_p += g->weight(i);
    c[level+1].resize(0); weight_c[level+1] = 0;
    for(uint it2 = c_i; it2 < c[level].size(); ++it2)
    {
      if(c[level][it2] != i && v->check(p, c[level][it2]))
      {
        c[level+1].push_back(c[level][it2]);
        weight_c[level+1] += g->weight(c[level][it2]);
      }
    }
    find_max(c, weight_c, p, weight_p, mu, v, g, res, level+1, start, time_lim);
    p.pop_back(); weight_p -= g->weight(i);
    v->undo_aux(p, i, c[level]);
  }
  return;
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
  int* mu = new int[n];
  for(uint i = 0; i < n; ++i)
    mu[i] = 0;

  int i;
  if(world_rank == 0) // master
  {
    // wait for slave to load and report readiness
    for(int i = 1; i < world_size; ++i)
      MPI_Recv(NULL, 0, MPI_BYTE, i, READYTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Using %d slaves\n", world_size);
    printf("All slaves are ready!\n");
    int cur_node = n-1;
    int rec_left = cur_node;
    // sending initial tasks
    for(int i = 1; cur_node >= 0 && i < world_size; ++i)
    {
      MPI_Send(&cur_node, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
      cur_node--;
    }
    // main loop, dynamic scheduling
    while(cur_node >= 0 || rec_left >= 0)
    {
      int buf[2];
      MPI_Status st;
      MPI_Recv(&buf, 2, MPI_INT, MPI_ANY_SOURCE, DONETAG, MPI_COMM_WORLD, &st);
      rec_left--;
      printf("Master : received result from slave %d : mu[%d] = %d\n", st.MPI_SOURCE, buf[0], buf[1]);
      mu[buf[0]] = buf[1];
      // now send this data to everybody
      //FIXME non-blockingly
      for(int i = 1; i < world_size; ++i)
      {
        if(i == st.MPI_SOURCE)
          continue;
        MPI_Send(&buf, 2, MPI_INT, i, UPDATETAG, MPI_COMM_WORLD);
      }
      // if there are more jobs, send to this slave
      if(cur_node >= 0)
      {
        MPI_Send(&cur_node, 1, MPI_INT, st.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        cur_node--;
      }
    }

    // exit
    for(int i = 1; i < world_size; ++i)
      MPI_Send(NULL, 0, MPI_BYTE, i, EXITTAG, MPI_COMM_WORLD);
    chrono::duration<double> d = chrono::steady_clock::now() - start;
    printf("Master : rds: time elapsed = %.8lf secs\n", d.count());
    // dumping mu to a file
    if(mu[g->nr_nodes-1] < 0)
      mu[g->nr_nodes-1] = -mu[g->nr_nodes-1];
    for(int j = g->nr_nodes-2; j >= 0; --j)
    {
      if(mu[j] > 0 && mu[j+1] > 0)
        continue;
      if(mu[j] < 0 && mu[j+1] > 0)
        mu[j] = max(mu[j+1], -mu[j]);
      else
        break;
     }
     FILE* f = fopen("final_mu.txt", "w");
     for(int j = g->nr_nodes-1; j >= 0; --j)
       fprintf(f, "mu[%d] = %d\n", j, mu[j]);
     fclose(f);
     printf("Master done, result is %u\n", mu[0]);
     MPI_Finalize();
     return mu[0];
  } else { // slave
    MPI_Send(NULL, 0, MPI_BYTE, 0, READYTAG, MPI_COMM_WORLD);

    int current_work = -1;
    lb = 0;

    while(!should_exit)
    {
      int buf[2];
      MPI_Status st;
      MPI_Recv(&buf, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

      vector<vector<uint> > c(g->nr_nodes); vector<uint> weight_c(g->nr_nodes, 0);
      vector<uint> p; uint weight_p = 0;
      int i;
      switch(st.MPI_TAG)
      {
      case EXITTAG:
        should_exit = true;
        printf("Slave %d : exiting...\n", world_rank);
        break;
      case UPDATETAG:
        // update mu
        mu[buf[0]] = buf[1];
        printf("Slave %d : bound update, now mu[%d] = %d\n", world_rank, buf[0], buf[1]);
        if(mu[g->nr_nodes-1] < 0)
          mu[g->nr_nodes-1] = -mu[g->nr_nodes-1];
        // update all -1
        for(int j = g->nr_nodes-2; j >= 0; --j)
        {
          if(mu[j] > 0 && mu[j+1] > 0)
            continue;
          if(mu[j] < 0 && mu[j+1] > 0)
            mu[j] = max(mu[j+1], -mu[j]);
          else
            break;
        }
        break;
      case WORKTAG:
        printf("Slave %d : received work %d\n", world_rank, buf[0]);
        current_work = buf[0];
        #pragma omp parallel
        {
          #pragma omp single
          printf("Slave %d : Using up to %d threads (OMP)\n", world_rank, omp_get_num_threads()); //FIXME might change in the future
        }

        //NB: SMART recomputing of LB
        lb = 0;
        for(int j = g->nr_nodes-1; j > current_work; --j)
        {
          lb = max(lb.load(), abs(mu[j]));
        }

        //NOW RDS-SERIAL CODE FOLLOWS
        //FIXME FACTOR
        i = current_work;
        // form candidate set
        // take vertices from v \in {i+1, n} for which pair (i,v) satisfies \Pi
        // first iteration c is empty, that must set bound to 1
        for(uint j = 0; j < g->nr_nodes; ++j)
        {
          c[j].reserve(g->nr_nodes);
          c[j].resize(0);
        }
        for(uint j = i+1; j < n; ++j)
        {
          if(v->check_pair(i, j))
          {
            // add to C
            c[0].push_back(j);
            weight_c[0] += g->weight(j);
          }
        }
        p.push_back(i); weight_p += g->weight(i);
        printf("Slave %d : i = %u, c.size = %lu\n", world_rank, i, c[0].size());
        // run for level = 0 manually with respect to the thread number
        if(c[0].size() == 0)
        {
          if(weight_p > lb)
          {
            mu[i] = weight_p;
            lb.store(max(lb.load(), weight_p));
//            res = p;
          }
          else
            mu[i] = lb.load();
        } else {
        #pragma omp parallel
        {
          // clone for separate threads
          verifier* v_ = v->clone();
          v_->init_aux(i, c[0]);
          vector<vector<uint> > c_(c); vector<uint> weight_c_(weight_c);
          vector<uint> p_(p); uint weight_p_ = weight_p;

          uint thread_i = omp_get_thread_num();
          uint num_threads = omp_get_num_threads();

          uint mu_i = 0;

          for(uint c_i = thread_i; c_i < c_[0].size() && !should_exit; c_i += num_threads) // split by threads
          {
            // adjust weight_c
            // we remove nodes [0; c_i), adjust weight accordingly
            for(uint j = 0; j < num_threads; ++j)
              if(c_i >= j+1)
                weight_c_[0] -= g->weight(c_i - j - 1);
            if(weight_c_[0] + weight_p_ <= lb) // Prune 1
            {
              mu_i = lb.load();
              break;
            }
            uint i_ = c_[0][c_i];
            if(mu[i_] > 0 && mu[i_] + weight_p_ <= lb) // Prune 2
            {
              mu_i = lb.load();
              break;
            } else {
              v_->prepare_aux(p_, i_, c_[0]);
              p_.push_back(i_); weight_p_ += g->weight(i_);
              c_[1].resize(0); weight_c_[1] = 0;
              for(uint it2 = c_i; it2 < c_[0].size(); ++it2)
              {
                if(c_[0][it2] != i_ && v_->check(p_, c_[0][it2])) //TODO only swap check?
                {
                  c_[1].push_back(c_[0][it2]);
                  weight_c_[1] += g->weight(c_[0][it2]);
                }
              }
              find_max(c_, weight_c_, p_, weight_p_, mu, v_, g, res, 1, start, time_lim);
              p_.pop_back(); weight_p_ -= g->weight(i_);
              v_->undo_aux(p_, i_, c_[0]);
            }
          }
          mu_i = lb.load();
          #pragma omp critical
          {
            mu[i] = max(mu_i, mu[i]);
          }
          v_->free_aux();
          delete v_;
        } // pragma omp parallel
        } // else
        mu[i] = -mu[i];
        printf("Slave %d : calculated mu[%d] = %d\n", world_rank, i, mu[i]);
        int buf[2];
        buf[0] = i;
        buf[1] = mu[i];
        MPI_Send(&buf, 2, MPI_INT, 0, DONETAG, MPI_COMM_WORLD);
        break;
      default:
        printf("Slave %d : wrong mesage tag %d, ingoring...\n", world_rank, st.MPI_TAG);
        break;
      } // switch
    } // while
    MPI_Finalize();
  } // world_rank
  delete [] mu;
  return 0;
}
