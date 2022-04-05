// No need to modify this file.

#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

int main(int argc, char *argv[]) {
    int procID;
    int nproc;
    double startTime;
    double endTime;
    double prob = 0.1;
    int numIterations = 5;
    char *inputFilename = NULL;
    int opt = 0;

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Read command line arguments
    do {
        opt = getopt(argc, argv, "f:p:i:");
        switch (opt) {
        case 'f':
            inputFilename = optarg;
            break;

        case 'p':
            prob = atof(optarg);
            break;

        case 'i':
            numIterations = atoi(optarg);
            break;

        case -1:
            break;

        default:
            break;
        }
    } while (opt != -1);

    if (inputFilename == NULL) {
        printf("Usage: %s -f <filename> [-p <P>] [-i <N_iters>]\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    // Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Run computation
    startTime = MPI_Wtime();
    compute(procID, nproc, inputFilename, prob, numIterations);
    endTime = MPI_Wtime();

    // Cleanup
    MPI_Finalize();
    printf("Elapsed time for proc %d: %f\n", procID, endTime - startTime);
    return 0;
}