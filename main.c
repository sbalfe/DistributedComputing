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

double set_average(const double *arr, int y, int x, uint arrSize){
    double above = arr[arrSize * (y-1) + x];
    double left = arr[arrSize * y + (x-1)];
    double below = arr[arrSize * (y+1) + x];
    double right = arr[arrSize * y + (x+1)];

    return (double) (above + left + below + right) / (double) 4;
}

void array_passthrough(context_t *context){

    int array_offset = context->displacements[context->rank];

   // printf("rank: %d\n", context->rank);
   // printf("array_offset: %d\n", array_offset);
    for (int i = 0; i < context->block_size[context->rank] ; ++i){

        // figure out what position we are at in the array in terms of rows and columns
        int y = (array_offset + i) / (int) context->array_size;
        int x = (array_offset + i) % (int) context->array_size;

        // check for borders
        if (y == 0 || x == 0 || y == context->array_size - 1 || x == context->array_size - 1){
            //printf("y: %d\n", y);
            //printf("x: %d\n", x);
            //printf("border check succeeded, continuing\n");
            continue;
        }

      //  printf("border check failed, changing value\n");

        double old_value = context->local_buffer[i];
        double new_value = set_average(context->input_buffer , y,x, context->array_size);
        //printf("changing local buffer value\n");
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

    int name_len;

    context_t *context = malloc(sizeof(context_t));
    context->array_size = (int) strtol(argv[1], 0,10);
    context->precision = strtof(argv[2], 0);

    MPI_Comm_rank(MPI_COMM_WORLD, &context->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &context->n_processors);

    char name[MPI_MAX_PROCESSOR_NAME];


    context->block_size = malloc((ssize_t) sizeof(double) * context->n_processors);
    context->displacements = malloc((ssize_t) sizeof(double) * context->n_processors);
    context->input_buffer = malloc((ssize_t) sizeof(double) * ((ssize_t) pow(context->array_size,2)));

    uint remainder = (uint) (pow(context->array_size,2)) % context->n_processors;

    if (context->rank == 0) {
        int sum = 0;
        for (uint i = 0; i < context->n_processors; ++i) {
            context->block_size[i] = (int) (pow(context->array_size,2)) / context->n_processors;
            if (remainder > 0) {
                context->block_size[i]++;
                remainder--;
            }
            context->displacements[i] = sum;
            printf("sum: %d\n", sum);
            sum += context->block_size[i];
        }
    }

    MPI_Bcast(context->block_size, context->n_processors, MPI_INT, 0, MPI_COMM_WORLD);

    if (context->rank == 6){
        printf("displacement value: %d\n", context->displacements[context->rank]);
        printf("processor count: %d\n", context->n_processors);
    }
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
                        context->input_buffer[context->array_size * y + x] = 0;
                    }
                }
            }
        }
    }

    context->complete = 1;


   while(1) {
        MPI_Scatterv(context->input_buffer, (int *) context->block_size, (int *) context->displacements, MPI_DOUBLE,
                     context->local_buffer, (int) context->block_size[context->rank], MPI_DOUBLE,
                     0, MPI_COMM_WORLD);

        array_passthrough(context);

        if (context->rank == 0) {
        MPI_Gatherv(context->local_buffer, context->block_size[context->rank],
                    MPI_DOUBLE, context->input_buffer,
                    (int*) context->block_size, (int*)context->displacements,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
        MPI_Gatherv(context->local_buffer, context->block_size[context->rank],
                    MPI_DOUBLE, NULL,
                    (int*)context->block_size, (int*)context->displacements,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        if (context->complete == 0) {
            context->complete = 1;
            MPI_Bcast(&context->complete, 1, MPI_INT, context->rank, MPI_COMM_WORLD);
        }
        else {
            break;
        }
    }

    if (context->rank == 0) {
        print_array(context->input_buffer, context->array_size);
    }

    MPI_Finalize();
    return 0;
}
