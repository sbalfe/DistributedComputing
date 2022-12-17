#include <stdio.h>
#include <mpi.h>
// test
int main(int argc, char **argv) {
    int rc, myrank, nproc, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];

    // set up the system. This must always be done. batch processing system (e.g. SLURM) starts the
    // processes on all the processors. This does the setup for connections between them.
    rc = MPI_Init(&argc, &argv);

    // always check rc to ensure there was no errors
    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // MPI_COMM_WORLD, system can be divided into subsets of processors called communicators
    // The world communicator is all processors
    // MPI_COMM_SELF refers to just the calling processor.

    // MPI_Comm_rank , each process in a communicator has a unique rank within that communicator
    // This is just an integer from 0 to size of the communicator - 1, therefore for world the rank ranges from 0 to total number of processors -1 .
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // obtain the size of the communicator
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int n[5];
    // all processors run the same code (SPMD).
    // This is how we get different things happening on different processors.
    if (myrank == 0) {
        printf("main reports %d procs\n", nproc);
        MPI_Send(n,5, MPI_INT, 1, 99, MPI_COMM_WORLD);
    }
    else if (myrank == 1){
        MPI_Status stat;
        MPI_Recv(n,5, MPI_INT, 0,99, MPI_COMM_WORLD, &stat);
    }

    namelen = MPI_MAX_PROCESSOR_NAME;
    MPI_Get_processor_name(name, &namelen);
    printf("hello world %d from ’%s’\n", myrank, name);

    // all processors must always call this to tidy up their MPI state.
    MPI_Finalize();

    return 0;
}
