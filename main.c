#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// test


typedef struct context {

} context_t ;

int main(int argc, char **argv) {

    int rc, my_rank, n_processors, name_len;
    char name[MPI_MAX_PROCESSOR_NAME];

    rc = MPI_Init(&argc, &argv);

    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI test program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processors);

    int *p_buffer;

    if (my_rank == 0){
        p_buffer = malloc(sizeof(int) * 2);
        //p_buffer[0] = 1;
        //p_buffer[1] = 2;
    }

   // MPI_Bcast(p_buffer, 2, MPI_INT, 0, MPI_COMM_WORLD);

    //printf("%d: first item is %d\n", my_rank, p_buffer[1]);


    name_len = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &name_len);


    MPI_Finalize();

    return 0;
}
