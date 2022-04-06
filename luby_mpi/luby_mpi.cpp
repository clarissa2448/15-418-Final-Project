#include <assert.h>
#include <algorithm> 
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <math.h> 
#include <cstring>
#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <random>
#include <iostream>
#include <vector>
#include <set>


#define BUFFER_LENGTH 300000
  
using namespace std;




// Parameters:
// adj_matrix: the adjacency matrix to which we are adding an edge
// i_start, i_end: the i range of the adjacency matrix to which we are
//                  adding an edge (not including i_end)
// j_start, j_end: the j range of the adjacency matrix to which we are
//                  adding an edge (not including j_end)
// a, b, c, d: As described in generate_rmat_graph
// We assuming i_end - i_start and j_end - j_start are powers of 2.
// Returns:
// True if a new edge has been added to the adjacency matrix, and false otherwise.
std::pair<int, int> generate_edge(int i_start, int i_end, int j_start, int j_end, 
                                    double a, double b, double c, double d, 
                                    std::uniform_real_distribution<double> rng, std::mt19937 mersenne_twister) {

    if (i_end == i_start + 1 && j_end == j_start + 1) return std::make_pair(i_start, j_start);

    // Choose which square to put the edge in randomly.
    double prob = rng(mersenne_twister);
    if (prob < a) {
        // Upper left
        i_end = i_start + (i_end - i_start)/2;
        j_end = j_start + (j_end - j_start)/2;
    }
    else if (prob < a + b) {
        // Upper right
        i_end = i_start + (i_end - i_start)/2;
        j_start = j_start + (j_end - j_start)/2;
    }
    else if (prob < a + b + c) {
        // Lower left
        i_start = i_start + (i_end - i_start)/2;
        j_end = j_start + (j_end - j_start)/2;
    }
    else {
        // Lower right
        i_start = i_start + (i_end - i_start)/2;
        j_start = j_start + (j_end - j_start)/2;
    }
    return generate_edge(i_start, i_end, j_start, j_end,
                            a, b, c, d, rng, mersenne_twister);
}



// n: Log (base 2) of the number of vertices
// E: number of edges in the generated graph
// a: Probability of edge going into upper left of adj. matrix
// b: Probability of edge going into upper right of adj. matrix
// c: Probability of edge going into lower left of adj. matrix (same as b)
// d: Probability of edge going into lower right of adj. matrix
// This creates the adjacency list of an undirected graph.
set<int> * generate_rmat_graph(int n, int E, double a, double b, double d) {

    std::random_device seed;
    std::mt19937 mersenne_twister(seed());
    std::uniform_real_distribution<double> rng(0.0, 1.0);
    double c = b; //probability of the lower left square of adjacency matrix.

    int N = (1 << n); // number of vertices
    set<int> * adjacency_list = (std::set<int> *) calloc(N, sizeof(set<int>));
    int num_edges_left = E;
    while (num_edges_left > 0) {

        std::pair<int, int> new_edge = generate_edge(0, N, 0, N, a, b, c, d, rng, mersenne_twister);
        int u = new_edge.first;
        int v = new_edge.second;
        if (u < v) continue; // Ignore all the edges that are above the main diagonal.

        // Check if this is a new edge
        bool edge_seen = false;
        set<int, greater<int> >::iterator itr;
        for (itr = adjacency_list[u].begin(); itr != adjacency_list[u].end(); itr++) {
            int neighbor = *itr;
            if (neighbor == v) {
                edge_seen = true;
                break;
            }
        }

        printf("%d %d\n", u, v);

        // Add if this is edge has not been seen
        if (!edge_seen) {
            adjacency_list[u].insert(v);
            adjacency_list[v].insert(u);
            num_edges_left--;
        }
    }

    return adjacency_list;
}

// Need set to array function
int* setToArray(set<int> s) {
    int* arr = (int*) calloc(s.size(), sizeof(int));
    return arr;
}

// Need array to set function
set<int> arrayToSet(int* arr, int size) {
    set<int> A;
    for (int i = 0; i < size; i ++){
        A.insert(arr[i]);
    }
    return A;
}

void set_union_yes(set<int> &A, set<int> &B) {
    set<int, greater<int> >::iterator itr;
    for (itr = B.begin(); itr != B.end(); itr++) {
        A.insert(*itr);
    }
}

set<int> set_intersect(set<int> &A, set<int> &B) {
    set<int, greater<int> >::iterator itr;
    set<int> C;
    for (itr = B.begin(); itr != B.end(); itr++) {
        if ((A.count(*itr))){
            C.insert(*itr);
        }
    }
    return C;
}

void set_minus(set<int> &A, set<int>&B) {
    set<int, greater<int> >::iterator itr;
    for (itr = B.begin(); itr != B.end(); itr++) {
        A.erase(*itr);
    }
}

set<int> luby_algorithm(int procID, int nproc, int n, int E, set<int> * adj_list) {
    set<int> A; 
    int* A_arr = (int*)calloc(BUFFER_LENGTH, sizeof(int));
    int* M = (int*)calloc(BUFFER_LENGTH, sizeof(int));
    // set<int> M;
    int lenA;
    int c = 3;
    const int root = 0;
    for (int i = 0; i < n; i ++) {
        A.insert(i);
        A_arr[i] = i;
    }
    lenA = A.size();
    int M_prime_size;
    set<int> M_prime;
    while (lenA!= 0) {
        float* random = (float*) calloc(lenA, sizeof(float));
        int nc = pow(n, c);
        if (procID == root) {
            for (int i = 0; i < lenA; i++) {
                random[i] = rand() % nc + 1;
            }
        }
        // Broadcast A and random
        MPI_Bcast(random, lenA, MPI_FLOAT, root, MPI_COMM_WORLD);
        // For each v in A do: 
        int num_groups; 
        if (lenA % nproc == 0) num_groups = lenA / nproc;
        else num_groups = lenA / nproc + 1;
        // Compare rv to all ru in N(v) and A and if rv is bigger, add v to M prime 
        for (int group_idx = 0; group_idx < num_groups/nproc; group_idx ++) {
            int v = group_idx * nproc + procID; 
            bool flag = true;

            set<int> intersect = set_intersect(adj_list[v], A);
            set<int, greater<int> >::iterator itr;
            for (itr = intersect.begin(); itr != intersect.end(); itr++) { // Check all u in N(v)
                if (random[v] <= random[*itr]) { // u is a neighbor
                    flag = false;
                    break; 
                }
            }
            if (flag) M_prime.insert(v); // Add v to M prime because it passes all checks
        }

        // Use this source for set union: https://www.cplusplus.com/reference/algorithm/set_union/
        // Use this source for send/recv: https://mpitutorial.com/tutorials/mpi-send-and-receive/
        // Use this source for set intersect: https://www.cplusplus.com/reference/algorithm/set_intersection/

        // Send M prime to root, who unions them all together and sets equal to M
        int* M_prime_arr;
        if (procID != root) {
            M_prime_size = M_prime.size();
            M_prime_arr = setToArray(M_prime);
             // Send size of set first 
            MPI_Send(&M_prime_size, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
            MPI_Send(M_prime_arr, M_prime_size, MPI_INT, root, 0, MPI_COMM_WORLD);
            
        } else {
            for (int i = root + 1; i < nproc; i ++) {
                // Receive length first
                MPI_Recv(&M_prime_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                M_prime_arr = (int*) calloc(M_prime_size, sizeof(int));
                // Receive Set 
                MPI_Recv(M_prime_arr, M_prime_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                set<int> M_prime_received = arrayToSet(M_prime_arr, M_prime_size);
                // Do a set union 
                set_union_yes(M_prime, M_prime_received);
            }
        }
        if (procID == root) {
            // Get new definition of A and broadcast
            // TODO: I dont think M is needed
            set<int> mPrimeNeighbors;
            set<int, greater<int> >::iterator itr;
            for (itr = M_prime.begin(); itr != M_prime.end(); itr++) {
                set<int, greater<int> >::iterator itrTwo;
                for (itrTwo = adj_list[*itr].begin(); itrTwo != adj_list[*itr].end(); itrTwo++){
                    mPrimeNeighbors.insert(*itrTwo);
                }
            }
            set_union_yes(mPrimeNeighbors, M_prime);

            // Get Set intersection
            set_minus(A, mPrimeNeighbors);
            // Broadcast lenA 
            lenA = A.size();
            A_arr = setToArray(A);

        }
         MPI_Bcast(&lenA, 1, MPI_INT, root, MPI_COMM_WORLD);
         MPI_Bcast(&A_arr, lenA, MPI_INT, root, MPI_COMM_WORLD);
         A = arrayToSet(A_arr, lenA);
        
    }
    // Broadcast Final M to output
    int final_m;
    if (procID == root) {
        final_m = M_prime.size();
        M = setToArray(M_prime);
    }
    
    MPI_Bcast(&final_m, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(M, final_m, MPI_INT, root, MPI_COMM_WORLD);
    set<int> M_set = arrayToSet(M, final_m);
    return M_set;
}


// int main(int argc, const char *argv[]) {
//     using namespace std::chrono;
//     typedef std::chrono::high_resolution_clock Clock;
//     typedef std::chrono::duration<double> dsec;

//     auto init_start = Clock::now();
//     double init_time = 0;

//     _argc = argc - 1;
//     _argv = argv + 1;

//     int n = get_option_int("-n", 10);
//     int E = get_option_int("-E", 100);
//     double a = get_option_float("-a", 0.1f);
//     double b = get_option_float("-b", 0.3f);
//     double d = get_option_float("-d", 0.3f);

//     printf("Number of Nodes: %d Number of Edges: %d\n", n, E);
//     printf("Probability Params: %lf %lf %lf.\n", a, b, d);

//     // 1. Generate random graph (adjacency matrix format) using 
//     // R-MAT (random graph model due to Chakrabarti, Zhan and 
//     // Faloutsos, and used in the work of Blelloch, Fineman and 
//     // Shun)
//     // Description of R-MAT method available at:
//     // https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf
//     bool ** adj_matrix = generate_undirected_graph(n, E, a, b, d);

//     for (int i = 0; i < n; i ++) {
//         for (int j = 0; j < n; j ++) {
//             printf("%d ", adj_matrix[i][j]);
//         }
//         printf("\n");
//     }

//     init_time += duration_cast<dsec>(Clock::now() - init_start).count();
//     printf("Initialization Time: %lf.\n", init_time);

//     auto compute_start = Clock::now();
//     double compute_time = 0;


//     compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
//     printf("Computation Time: %lf.\n", compute_time);


//     return 0;
// }