#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

struct Context {
    uint array_size;
    int n_processors;
    int rank;
    double precision;
    int *block_size;
    int *displacements;
    double *local_buffer;
} typedef context_t;


void print_array(void* arr, uint arrSize) {

    double *array = (double *) arr;
    for (uint y = 0; y < arrSize; ++y) {

        for (int x = 0; x < arrSize; ++x) {
            printf("%f ", array[arrSize * y + x]);
        }
        printf("\n");
    }
}

int main(int argc, char **argv) {

    int rc = MPI_Init(&argc, &argv);
    int name_len;

    context_t context;
    context.array_size = (int) strtol(argv[1], 0,10);
    context.precision = strtof(argv[2], 0);

    MPI_Comm_rank(MPI_COMM_WORLD, &context.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &context.n_processors);

    char name[MPI_MAX_PROCESSOR_NAME];

    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI test program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    double *input_buffer = malloc(sizeof(double) * ((size_t) pow(context.array_size,2)));
    context.block_size = malloc(sizeof(double) * context.n_processors);
    context.displacements = malloc(sizeof(double) * context.n_processors);

    uint remainder = (uint) (pow(context.array_size,2)) % context.n_processors;

    int sum = 0;
    if (context.rank == 0) {
        for (uint i = 0; i < context.n_processors; ++i) {
            context.block_size[i] = (int) (pow(context.array_size,2)) / context.n_processors;
            if (remainder > 0) {
                context.block_size[i]++;
                remainder--;
            }
            context.displacements[i] = sum;
            sum += context.block_size[i];
        }
    }

    MPI_Bcast(context.block_size, context.n_processors, MPI_INT, 0, MPI_COMM_WORLD);

    context.local_buffer = malloc(sizeof(double) * context.block_size[context.rank]);

    // make one processor allocate the array.
    if (context.rank == 0) {
        for (int y = 0; y < context.array_size; ++y) {
            if (y == 0) {
                for (int x = 0; x < context.array_size; ++x) {
                    input_buffer[context.array_size * y + x] = 1;
                }
            } else {
                for (int x = 0; x < context.array_size; ++x) {
                    if (x == 0) {
                        input_buffer[context.array_size * y + x] = 1;
                    } else {
                        input_buffer[context.array_size * y + x] = 0;
                    }
                }
            }
        }
    }

    MPI_Scatterv(input_buffer, (int *) context.block_size , (int *) context.displacements, MPI_DOUBLE,
                context.local_buffer , (int) context.block_size[context.rank], MPI_DOUBLE,
                0, MPI_COMM_WORLD);


    if (context.rank == 9){
        for (uint r = 0; r < context.n_processors; ++r){
            printf("RANK %u : [",r);
            for (uint i = 0; i < context.block_size[r]; ++i){
                printf(" %f, ", context.local_buffer[i]);
            }
            printf("]\n");
        }
    }



    MPI_Finalize();
    return 0;
}
