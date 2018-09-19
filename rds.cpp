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

inline void update_mu(graph* g, int* mu, int buf0, int buf1)
{
  mu[buf0] = buf1;
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
}

void find_max(vector<vertex_set>& c, vertex_set& p, int* mu, verifier *v, graph* g, vector<uint>& res, int level, const chrono::time_point<chrono::steady_clock> start, const uint time_lim)
{
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  vertex_set& curC = c[level];
  if(should_exit)
    return;
  if(curC.empty())
  {
    if(p.weight > lb)
    {
      #pragma omp critical (lbupdate)
      {
      res = p; //copy
      lb.store(max(lb.load(), p.weight));
      }
      return;
    }
    else
      return;
  }

  vertex_set& nextC = c[level+1];
  for(uint c_i = 0; c_i < curC.size(); ++c_i)
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
      int flag = 0;
      MPI_Status st;
      auto error = MPI_Iprobe(0, UPDATETAG, MPI_COMM_WORLD, &flag, &st); // TODO EXITTAG
      if(flag)
      {
        fprintf(stderr, "Slave %d probe: Flag value was %d, error code was %d (%d), MPI_status is: SOURCE = %d, TAG = %d, ERROR = %d\n", world_rank, flag, error, MPI_SUCCESS, st.MPI_SOURCE, st.MPI_TAG, st.MPI_ERROR);
        // process UPDATETAG
        int32_t buf[2];
        MPI_Recv(buf, 2, MPI_INT32_T, 0, UPDATETAG, MPI_COMM_WORLD, &st);
        fprintf(stderr, "Slave %d : Got a fancy mu update mu[%d] = %d!\nSlave %d: MPI_status is: SOURCE = %d, TAG = %d, ERROR = %d\n", world_rank, buf[0], buf[1], world_rank, st.MPI_SOURCE, st.MPI_TAG, st.MPI_ERROR);
        #pragma omp critical (muupdate)
        {
          update_mu(g, mu, buf[0], buf[1]);
        }
      }
    }

    if(curC.weight + p.weight <= lb) // Prune 1
    {
      return;
    }

    uint i = curC[c_i];
    if(mu[i] > 0 && mu[i] + p.weight <= lb) // Prune 2
    {
      return;
    }

    curC.weight -= g->weight(i);
//    NB: exploit that we adding only 1 vertex to p
//    thus verifier can prepare some info using prev calculations
    v->prepare_aux(p, i, curC);
    p.add_vertex(i, g->weight(i));
    nextC.clear();

    for(uint it2 = c_i; it2 < curC.size(); ++it2)
    {
      uint u = curC[it2];
      if(u != i && v->check(p, u))
      {
        nextC.add_vertex(u, g->weight(u));
      }
    }
    find_max(c, p, mu, v, g, res, level+1, start, time_lim);
    p.pop_vertex(g->weight(i));
    v->undo_aux(p, i, curC);
  }
  return;
}



uint rds(verifier* v, graph* g, vector<uint>& res, uint time_lim, bool slave_out)
{
  MPI_Init(NULL, NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if(world_size < 2)
  {
    fprintf(stderr, "MPI program requires at least 2 nodes (got %d), bailing out...\n", world_size);
    return 0;
  }
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
    printf("Using %d slaves\n", world_size);
    // wait for slave to load and report readiness
    for(int i = 1; i < world_size; ++i)
      MPI_Recv(NULL, 0, MPI_INT32_T, i, READYTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("All slaves are ready!\n");
    int32_t cur_node = n-1;
    int rec_left = cur_node;
    // sending initial tasks
    for(int i = 1; cur_node >= 0 && i < world_size; ++i)
    {
      MPI_Send(&cur_node, 1, MPI_INT32_T, i, WORKTAG, MPI_COMM_WORLD);
      cur_node--;
    }
    // main loop, dynamic scheduling
    while(cur_node >= 0 || rec_left >= 0)
    {
      int32_t buf[2];
      MPI_Status st;
      MPI_Recv(buf, 2, MPI_INT32_T, MPI_ANY_SOURCE, DONETAG, MPI_COMM_WORLD, &st);
      rec_left--;
      printf("Master : received result from slave %d : mu[%d] = %d\n", st.MPI_SOURCE, buf[0], buf[1]);
      mu[buf[0]] = buf[1];
      // now send this data to everybody
      //FIXME non-blockingly
      for(int i = 1; i < world_size; ++i)
      {
/*        if(i == st.MPI_SOURCE)
          continue;*/
        MPI_Send(buf, 2, MPI_INT32_T, i, UPDATETAG, MPI_COMM_WORLD);
      }
      // if there are more jobs, send to this slave
      if(cur_node >= 0)
      {
        MPI_Send(&cur_node, 1, MPI_INT32_T, st.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        cur_node--;
      }
    }

    // exit
    for(int i = 1; i < world_size; ++i)
      MPI_Send(NULL, 0, MPI_INT32_T, i, EXITTAG, MPI_COMM_WORLD);
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
     int ret = mu[0];
     delete [] mu;
     return ret;
  } else { // slave
    MPI_Send(NULL, 0, MPI_INT32_T, 0, READYTAG, MPI_COMM_WORLD);

    int current_work = -1;
    lb = 0;

    if(slave_out) {
    #pragma omp parallel
    {
      #pragma omp single
      printf("Slave %d : Using up to %d threads (OMP)\n", world_rank, omp_get_num_threads()); //FIXME might change in the future
    }
    }

    while(!should_exit)
    {
      int32_t buf[2];
      MPI_Status st;
      MPI_Recv(buf, 2, MPI_INT32_T, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

      vector<vertex_set> c(g->nr_nodes);
      vertex_set p;
      int i;
      switch(st.MPI_TAG)
      {
      case EXITTAG:
        should_exit = true;
        if(slave_out)
          printf("Slave %d : exiting...\n", world_rank);
        break;
      case UPDATETAG:
        // update mu
        update_mu(g, mu, buf[0], buf[1]);
        if(slave_out)
          printf("Slave %d : bound update, now mu[%d] = %d\n", world_rank, buf[0], buf[1]);
        break;
      case WORKTAG:
        if(slave_out)
          printf("Slave %d : received work %d\n", world_rank, buf[0]);
        current_work = buf[0];
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
        }
        for(uint j = i+1; j < n; ++j)
        {
          if(v->check_pair(i, j))
          {
            // add to C
            c[0].add_vertex(j, g->weight(j));
          }
        }
        p.add_vertex(i, g->weight(i));
        if(slave_out)
          printf("Slave %d : i = %u, c.size = %lu\n", world_rank, i, c[0].size());
        // run for level = 0 manually with respect to the thread number
        if(c[0].empty())
        {
          if(p.weight > lb)
          {
            mu[i] = p.weight;
            lb.store(max(lb.load(), p.weight));
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
          vector<vertex_set> c_(c);
          vertex_set p_(p);

          uint thread_i = omp_get_thread_num();
          uint num_threads = omp_get_num_threads();

          uint mu_i = 0;

          for(uint c_i = thread_i; c_i < c_[0].size() && !should_exit; c_i += num_threads) // split by threads
          {
            // adjust weight_c
            // we remove nodes [0; c_i), adjust weight accordingly
            for(uint j = 0; j < num_threads; ++j)
              if(c_i >= j+1)
                c_[0].weight -= g->weight(c_i - j - 1);
            if(c_[0].weight + p_.weight <= lb) // Prune 1
            {
              mu_i = lb.load();
              break;
            }
            uint i_ = c_[0][c_i];
            if(mu[i_] > 0 && mu[i_] + p_.weight <= lb) // Prune 2
            {
              mu_i = lb.load();
              break;
            } else {
              v_->prepare_aux(p_, i_, c_[0]);
              p_.add_vertex(i_, g->weight(i_));
              c_[1].clear();
              for(uint it2 = c_i; it2 < c_[0].size(); ++it2)
              {
                uint u = c_[0][it2];
                if(u != i_ && v_->check(p_, u)) //TODO only swap check?
                {
                  c_[1].add_vertex(u, g->weight(u));
                }
              }
              find_max(c_, p_, mu, v_, g, res, 1, start, time_lim);
              p_.pop_vertex(g->weight(i_));
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
        if(slave_out)
          printf("Slave %d : calculated mu[%d] = %d\n", world_rank, i, mu[i]);
        int32_t buf[2];
        buf[0] = i;
        buf[1] = mu[i];
        MPI_Send(buf, 2, MPI_INT32_T, 0, DONETAG, MPI_COMM_WORLD);
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
