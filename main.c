#include <stdio.h>
#include <stdlib.h>
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
    double *input_buffer;
    int complete;
    int colour;
} typedef context_t;


void print_array(void* input_buffer, uint array_size) {
    double *in = (double *) input_buffer;
    for (uint y = 0; y < array_size; ++y) {
        for (int x = 0; x < array_size; ++x) {
            printf("%f ", in[array_size * y + x]);
        }
        printf("\n");
    }
}

double calculate_average(const double *input_buffer, const int y, const int x, const uint array_size){

    double above = input_buffer[array_size * (y-1) + x];
    double left = input_buffer[array_size * y + (x-1)];
    double below = input_buffer[array_size * (y+1) + x];
    double right = input_buffer[array_size * y + (x+1)];

    return (double) (above + left + below + right) / (double) 4;
}

void array_relaxation(context_t *context){

    int array_offset = context->displacements[context->rank];

    for (int i = 0; i < context->block_size[context->rank] ; ++i){

        // figure out what position we are at in the initial input buffer.
        int y = (array_offset + i) / (int) context->array_size;
        int x = (array_offset + i) % (int) context->array_size;

        // check for borders
        if (y == 0 || x == 0 || y == context->array_size - 1 || x == context->array_size - 1){
            continue;
        }

        double old_value = context->local_buffer[i];
        double new_value = calculate_average(context->input_buffer ,y , x, context->array_size);
        context->local_buffer[i] = new_value;

        if (fabs(new_value - old_value) > context->precision){
            context->complete = 0;
        }
    }
}

int main(int argc, char **argv) {

    int rc = MPI_Init(&argc, &argv);

    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI test program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    context_t *context = malloc(sizeof(context_t));
    context->array_size = (int) strtol(argv[1], 0,10);
    context->precision = strtof(argv[2], 0);

    context->colour = context->rank / 4;
    MPI_Comm new_comm;
    MPI_Comm_split(MPI_COMM_WORLD, context->colour, context->rank, &new_comm);

    MPI_Comm_rank(new_comm, &context->rank);
    MPI_Comm_size(new_comm, &context->n_processors);

    context->block_size = malloc((ssize_t) sizeof(double) * context->n_processors);
    context->displacements = malloc((ssize_t) sizeof(double) * context->n_processors);
    context->input_buffer = malloc((ssize_t) sizeof(double) * ((ssize_t) pow(context->array_size,2)));

    if (context->rank == 0) {
        uint remainder = (uint) (pow(context->array_size,2)) % context->n_processors;
        int sum = 0;
        for (uint i = 0; i < context->n_processors; ++i) {
            context->block_size[i] = (int) (pow(context->array_size,2)) / context->n_processors;
            if (remainder > 0) {
                context->block_size[i]++;
                remainder--;
            }
            context->displacements[i] = sum;
            sum += context->block_size[i];
        }
    }

    MPI_Bcast(context->block_size, context->n_processors, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(context->displacements, context->n_processors, MPI_INT, 0, MPI_COMM_WORLD);

    context->local_buffer = malloc((ssize_t) sizeof(double) * context->block_size[context->rank]);

    // make one processor allocate the array
    if (context->rank == 0) {
        for (int y = 0; y < context->array_size; ++y) {
            if (y == 0) {
                for (int x = 0; x < context->array_size; ++x) {
                    context->input_buffer[context->array_size * y + x] = 1;
                }
            } else {
                for (int x = 0; x < context->array_size; ++x) {
                    if (x == 0) {
                        context->input_buffer[context->array_size * y + x] = 1;
                    } else {
                        context->input_buffer[context->array_size * y + x] = 2;
                    }
                }
            }
        }
        context->complete = 1;
    }

    MPI_Bcast(&context->complete, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(context->input_buffer,  (int) pow(context->array_size,2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int result = 0;

    while(1) {
        // subdivide the initial input buffer
        MPI_Scatterv(context->input_buffer, (int *) context->block_size, (int *) context->displacements, MPI_DOUBLE,
                     context->local_buffer, (int) context->block_size[context->rank], MPI_DOUBLE,
                     0, MPI_COMM_WORLD);

        array_relaxation(context);

        // gather all the values that have been altered.
        MPI_Gatherv(context->local_buffer, context->block_size[context->rank],
                    MPI_DOUBLE, context->input_buffer,
                    (int*) context->block_size, (int*)context->displacements,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Bcast(context->input_buffer,  (int) pow(context->array_size,2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        MPI_Allreduce(&context->complete, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (result == context->n_processors) {
            break;
        } else {
            // not all processors are complete , reset the flag and do another cycle.
            context->complete = 1;
        }
    }

    if (context->rank == 0) {
        print_array(context->input_buffer, context->array_size);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
