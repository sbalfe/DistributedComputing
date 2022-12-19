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

    int *p_send_buffer = calloc(n_processors, sizeof(int));

    int *p_receive_value = NULL;
    p_receive_value = malloc(sizeof(int));

    MPI_Scatter(p_send_buffer, 1 , MPI_INT, p_receive_value, 1, MPI_INT, 1, MPI_COMM_WORLD);

    printf("%d: hey the second value is %d\n", my_rank, p_receive_value[1]);

    name_len = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &name_len);

    MPI_Finalize();

    return 0;
}
