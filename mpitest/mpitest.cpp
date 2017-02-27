// master-slaves structure

#include <mpi.h>
#include <stdio.h>

#define READYTAG 1
#define WORKTAG  2
#define DONETAG  3
#define EXITTAG  4

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    if(world_rank == 0) // I'm master
    {
      // synchronize here, ensure everybody loaded the data
      for(int i = 1; i < world_size; ++i)
      {
        MPI_Recv(NULL, 0, MPI_BYTE, i, READYTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      printf("Master: all slaves ready!\n");

      // run main
      int nr_nodes = 200;
      int nodes_left = nr_nodes;

      for(int i = 1; nodes_left > 0 && i < world_size; ++i)
      {
        MPI_Send(&nodes_left, 1, MPI_INT, i, WORKTAG, MPI_COMM_WORLD);
        nodes_left--;
      }

      // when slave is ready then send more work if have any
      while(nodes_left > 0)
      {
        int buf[2];
        MPI_Status st;
        MPI_Recv(&buf, 2, MPI_INT, MPI_ANY_SOURCE, DONETAG, MPI_COMM_WORLD, &st);
        MPI_Send(&nodes_left, 1, MPI_INT, st.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
        nodes_left--;
      }

      // send exit to to slaves
      for(int i = 1; i < world_size; ++i)
      {
        MPI_Send(NULL, 0, MPI_BYTE, i, EXITTAG, MPI_COMM_WORLD);
      }
    } else { // I'm slave
      // load data
      ;
      // report read to master
      MPI_Send(NULL, 0, MPI_BYTE, 0, READYTAG, MPI_COMM_WORLD);

      // wait for work to be done
      int buf[1];
      MPI_Status st;
      bool should_exit = false;
      while(!should_exit)
      {
        MPI_Recv(&buf, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        switch(st.MPI_TAG)
        {
          case WORKTAG:
            printf("Slave %d : recevied work %d\n", world_rank, buf[0]);
            // do here
            ;
            // send result to master
            int res_buf[2];
            MPI_Send(&res_buf, 2, MPI_INT, 0, DONETAG, MPI_COMM_WORLD);
            break;
          case EXITTAG:
            should_exit = true;
            printf("Slave %d : exiting...\n", world_rank);
            break;
          default:
            printf("Slave %d : unknown message tag\n", world_rank);
            should_exit = true;
         }
      }
    }
    // Print off a hello world message
    /*printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);
    */
    // Finalize the MPI environment.
    MPI_Finalize();
}
