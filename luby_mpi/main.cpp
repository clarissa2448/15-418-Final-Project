// No need to modify this file.

#include "mpi.h"
#include "luby_mpi.h"
#include "luby_mpi_blocked_assignment.h"
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

#define BUFFER_LENGTH 300000
#define DEBUG false
using namespace std;

static int _argc;
static const char **_argv;

const char *get_option_string(const char *option_name, const char *default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return _argv[i + 1];
    return default_value;
}

int get_option_int(const char *option_name, int default_value, int _argc, char **_argv) {
  for (int i = _argc - 2; i >= 0; i -= 2){
        if (strcmp(_argv[i], option_name) == 0)
            return atoi(_argv[i + 1]);
    }
    return default_value;
}

float get_option_float(const char *option_name, float default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return (float)atof(_argv[i + 1]);
    return default_value;
}

/* Gen Dummy Adjacency List
0: <3>
1: <2, 3>
2: <1, 3> 
3: <0, 1, 2>
*/
set<int> * gen_dummy_adj_list() {
    set<int> * adj_list = (std::set<int> *) calloc(4, sizeof(set<int>));
    set<int> zero; 
    zero.insert(3);
    set<int> one; 
    one.insert(2);
    one.insert(3);
    set<int> two; 
    two.insert(1);
    two.insert(3);
    set<int> three;
    three.insert(0);
    three.insert(1);
    three.insert(2);
    adj_list[0] = zero;
    adj_list[1] = one;
    adj_list[2] = two;
    adj_list[3] = three;
    return adj_list;
}


// Write Input Adjacency List to File
void write_adj_list_to_file(set<int>* adj_list, int n, int E, int nproc) {
    ofstream adj_list_file;
    int N = 1 << n;
    char buffer[BUFFER_LENGTH];
    snprintf(buffer, BUFFER_LENGTH, "outputs/adj_list_%d_%d_%d.txt", n, E, nproc);
    adj_list_file.open(buffer);
    for (int i = 0; i < N; i++) {
        for (auto j = adj_list[i].begin(); j != adj_list[i].end(); j++) {
            adj_list_file << *j << " ";
        }
        adj_list_file << "\n";

        
    }
    adj_list_file.close();
}

// Write Output Independent Set to File
void write_mis_to_file(set<int> mis, int n, int E, int nproc) {
    ofstream mis_file;
    char buffer[BUFFER_LENGTH];
    snprintf(buffer, BUFFER_LENGTH, "outputs/maximal_indep_set_%d_%d_%d.txt", n, E, nproc);
    mis_file.open(buffer);
    for (auto j = mis.begin(); j != mis.end(); j++) {
        mis_file << *j << "\n";
    }
    mis_file.close();
}

int main(int argc, char *argv[]) {
    int procID;
    int nproc;
    double startTime;
    double endTime;

    int n = get_option_int("-n", 10, argc, argv);
    int E = get_option_int("-E", 100, argc, argv);
    double a = get_option_float("-a", 0.1f);
    double b = get_option_float("-b", 0.3f);
    double d = get_option_float("-d", 0.3f);
    int version = get_option_int("-v", 1, argc, argv);

    printf("Number of Nodes: %d Number of Edges: %d\n", n, E);
    printf("Probability Params: %lf %lf %lf.\n", a, b, d);
    printf("Version: %d\n", version);

    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    printf("Num processors: %d\n", nproc);

    // 1. Generate random graph (adjacency matrix format) using 
    // R-MAT (random graph model due to Chakrabarti, Zhan and 
    // Faloutsos, and used in the work of Blelloch, Fineman and 
    // Shun)
    // Description of R-MAT method available at:
    // https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf
    int root = 0;
    int N = (1 << n);
    set<int> * adj_list = (std::set<int> *) calloc(N, sizeof(set<int>)); 
    
    if (procID == root) {
        adj_list = generate_rmat_graph(n, E, a, b, d);
    }
    
    for (int i = 0; i < N; i++) {
        // Bcast length first
        int num_neighbors;
        if (procID == root) {
            num_neighbors = adj_list[i].size();
        }
        MPI_Bcast(&num_neighbors, 1, MPI_INT, root, MPI_COMM_WORLD);
        // Bcast set as array
        int * neighbors = (int *) calloc(num_neighbors, sizeof(int));
        if (procID == root) {
            neighbors = setToArray(adj_list[i]);
        }
        MPI_Bcast(neighbors, num_neighbors, MPI_INT, root, MPI_COMM_WORLD);
        set<int> neighbors_set = arrayToSet(neighbors, num_neighbors);
        adj_list[i] = neighbors_set;
    }

    if (DEBUG) {
        for (int i = 0; i < N; i ++) {
            printf("%d: ", i);
            for (auto itr = adj_list[i].begin(); itr != adj_list[i].end(); itr ++) {
                printf("%d ", *itr);
            }
            printf("\n");
        }
    }


    // Run computation
    startTime = MPI_Wtime();
    set<int> M;
    if (version == 1) M = luby_algorithm(procID, nproc, n, E, adj_list);
    else M = luby_algorithm_blocked_assignment(procID, nproc, n, E, adj_list);
    endTime = MPI_Wtime();

    if (DEBUG)
    {
        printf("Nodes in Maximal Independent Set\n");
        set<int, greater<int> >::iterator itr;
        for (itr = M.begin(); itr != M.end(); itr++) {
            printf("%d ", *itr);
        }
    }


    // Write to Files
    write_adj_list_to_file(adj_list, n, E, nproc);
    write_mis_to_file(M, n, E, nproc);

    // Cleanup
    MPI_Finalize();
    printf("Elapsed time for proc %d: %f\n", procID, endTime - startTime);
    return 0;
}