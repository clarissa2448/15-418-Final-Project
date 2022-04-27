#include <mpi.h>
#include <unistd.h>
#include <assert.h>
#include <algorithm> 
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <math.h> 
#include <cstring>
#include <stdio.h>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#define MAX 100
#define DEBUG false
using namespace std;

static int _argc;
static const char **_argv;
int    num_procs;
double x[MAX];			// Input array


void initialize() {
    for (int i = 0; i < MAX; i++) {
        x[i] = i;
        // printf("i %d\n", i);
    }
}

int main(int argc, char *argv[])
{
    int start, stop;
    int myid;
    int my_min;
    int others_min[100];		    // Save minimum separately
    int my_max;
    int others_max[100];		    // Save maximum separately
    
    MPI_Request rq_min[100], rq_max[100];  // Status variables

    MPI_Init(&argc,&argv);                      // Initialize

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // Get # processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);       // Get my rank (id)

    initialize();

    int edges = 20;
    // if (myid == 0) edges = 10 

    MPI_Barrier(MPI_COMM_WORLD);

    int n = MAX/num_procs;

    start = myid * n;
    printf("My id %d %d, %d\n", myid, n, start);

    if ( myid != (num_procs-1) )
    {
        stop = start + n;
    }
    else
    {
        stop = MAX;
    }
    

    /* --------------------------------------
    Now find the max. among my numbers
    -------------------------------------- */
    my_max = x[start];
    printf("PROCID: %d, 1Start: %d, Stop: %d\n", myid, start, stop);
    printf("PROCID: %d, 2Start: %d Stop: %d max: %d \n", myid, start, stop, my_max);
    if (start == 0) {
        printf("IT EQUALS 0\n");
    }
    for (int i = start+1; i < stop; i ++ )
    {
        // printf("x[i] %d", x[i]);
        if ( x[i] > my_max )
                my_max = x[i];
    }
    printf("PROCID: %d, 3Start: %d Stop: %d max: %d \n", myid, start, stop, my_max);

    for (int i = 0; i < edges; i++)
    {
        int v_process = myid ? 0: 1;
        MPI_Irecv(&others_max[i], 1, MPI_DOUBLE, v_process, 0, MPI_COMM_WORLD,       
        &rq_max[i]);
        MPI_Request send_request;
        MPI_Isend(&my_max, 1, MPI_DOUBLE, v_process, 0, MPI_COMM_WORLD,
                &send_request);

        
    }
    
    /* --------------------------------------
    Now synchronize to compute results
    -------------------------------------- */
    for (int i = 0; i < edges; i++)
    {
        MPI_Wait( &rq_max[i], NULL );
        printf("PROCID %d others max %d \n", myid, others_max[i]);
        if ( others_max[i] > my_max )
            my_max = others_max[i];
    }

    printf("PROCID %d max %d\n", myid, my_max);

    MPI_Finalize();
}
  
int main2(int argc, char *argv[])
{
    int start, stop;
    int myid;
    int my_min;
    int others_min[100];		    // Save minimum separately
    int my_max;
    int others_max[100];		    // Save maximum separately
    
    MPI_Request rq_min[100], rq_max[100];  // Status variables

    MPI_Init(&argc,&argv);                      // Initialize

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);  // Get # processors
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);       // Get my rank (id)

    initialize();
    MPI_Barrier(MPI_COMM_WORLD);

    int n = MAX/num_procs;

    start = myid * n;
    printf("My id %d %d, %d\n", myid, n, start);

    if ( myid != (num_procs-1) )
    {
        stop = start + n;
    }
    else
    {
        stop = MAX;
    }
    

    /* --------------------------------------
    Now find the max. among my numbers
    -------------------------------------- */
    my_max = x[start];
    printf("PROCID: %d, 1Start: %d, Stop: %d\n", myid, start, stop);
    printf("PROCID: %d, 2Start: %d Stop: %d max: %d \n", myid, start, stop, my_max);
    if (start == 0) {
        printf("IT EQUALS 0\n");
    }
    for (int i = start+1; i < stop; i ++ )
    {
        // printf("x[i] %d", x[i]);
        if ( x[i] > my_max )
                my_max = x[i];
    }
    printf("PROCID: %d, 3Start: %d Stop: %d max: %d \n", myid, start, stop, my_max);

    if ( myid == 0 )
    {
        /* -------------------------------------
            Get the max from others and compare
            ------------------------------------- */

            
        for (int i = 1; i < num_procs; i++)
        {
            MPI_Irecv(&others_max[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD,       
                &rq_max[i-1]);
        }
    }
    else
    {
        MPI_Request send_request;
        MPI_Isend(&my_max, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
            &send_request);
        // MPI_Isend(&my_max, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
        //     &rq_max[0]);
    }

    /* --------------------------------------
    Now synchronize to compute results
    -------------------------------------- */
    if ( myid == 0 )
    {
    // for (int i = 1; i < num_procs; i++)
    // {
    //     MPI_Wait( &rq_min[i], NULL );

    //     printf
    //     if ( others_min[i] < my_min )
    //         my_min = others_min[i];
    // }

    for (int i = 1; i < num_procs; i++)
    {
        MPI_Wait( &rq_max[i-1], NULL );
        printf("others max %d \n", others_max[i]);
        if ( others_max[i] > my_max )
            printf("others max bigger %d \n", others_max[i]);
            my_max = others_max[i];
    }

    printf("max %d\n", my_max);

    }
    else
    {  // The other processes must wait until their messages
        // has been received before exiting !!!
        // printf("PROCID %d Here?\n", myid);
        // printf("min %d\n", rq_min[0]);
        // MPI_Wait( &rq_min[0], NULL );
        // printf("seg %d\n", rq_max[0]);
        // MPI_Wait( &rq_max[0], NULL );
        printf("quack?\n");
    }

    MPI_Finalize();
}
