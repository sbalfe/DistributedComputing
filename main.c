#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

// test
typedef struct context context;

struct context {
    context *p_next;
};

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

    int *send_buffer = calloc(64, sizeof(int));
    int *receive_value = calloc(4,sizeof(int));

    MPI_Scatter(send_buffer, 4 , MPI_INT, receive_value, 4, MPI_INT, 1, MPI_COMM_WORLD);

    printf("%d: hey the second value is %d\n", my_rank, receive_value[0]);

    name_len = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &name_len);

    MPI_Finalize();

    return 0;
}
