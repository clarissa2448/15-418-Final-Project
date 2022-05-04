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
#define DEBUG false
#define PROB_POWER 3
  
using namespace std;


// Parameters:
// adj_matrix: the adjacency matrix to which we are adding an edge
// i_start, i_end: the i range of the adjacency matrix to which we are
//                  adding an edge (not including i_end)
// j_start, j_end: the j range of the adjacency matrix to which we are
//                  adding an edge (not including j_end)
// a, b, c, d: As described in generate_rmat_graph
// We assume i_end - i_start and j_end - j_start are powers of 2.
// Returns:
// An edge which may be added to the graph.
std::pair<int, int> generate_edge(int i_start, int i_end, int j_start, int j_end,
                                    double a, double b, double c, double d,
                                    std::uniform_real_distribution<double>& rng, std::mt19937& mersenne_twister) {
    
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
    double c = b; // probability of going into lower left square

    int N = (1 << n); // number of vertices
    set<int> * adjacency_list = (std::set<int> *) calloc(N, sizeof(set<int>));
    for (int i = 0; i < N; i++) {
        set<int> s;
        adjacency_list[i] = s;
    }

    int num_edges_left = E;
    while (num_edges_left > 0) {
        std::pair<int, int> new_edge = generate_edge(0, N, 0, N, a, b, c, d, rng, mersenne_twister);
        int u = new_edge.first;
        int v = new_edge.second;
        if (u <= v) continue; // Ignore all the edges that are above the main diagonal.

        // Check if this is a new edge
        if (adjacency_list[u].count(v)) continue;

        adjacency_list[u].insert(v);
        adjacency_list[v].insert(u);
        num_edges_left--;
    }
    return adjacency_list;
}

// Convert set to array
int* setToArray(set<int> s) {
    int* arr = (int*) calloc(s.size(), sizeof(int));
    int idx = 0;
    set<int, greater<int>>::iterator itr;
    for (itr = s.begin(); itr != s.end(); itr++) {
        arr[idx] = *itr;
        idx++;
    }
    return arr;
}

// Convert array to set
set<int> arrayToSet(int* arr, int size) {
    set<int> A;
    for (int i = 0; i < size; i++) {
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

// Returns A \ B
set<int> set_minus(set<int> A, set<int> B) {
    set<int> result;
    set<int, greater<int>>::iterator itr;
    for (itr = A.begin(); itr != A.end(); itr++) {
        if (B.count(*itr)) continue;
        result.insert(*itr);
    }
    return result;
}

/*
void set_minus(set<int> &A, set<int>&B) {
    set<int, greater<int> >::iterator itr;
    for (itr = B.begin(); itr != B.end(); itr++) {
        A.erase(*itr);
    }
}
*/

void print_set(set<int> A) {
    set<int, greater<int>>::iterator itr;
    for (itr = A.begin(); itr != A.end(); itr++) {
        printf("%d ", *itr);
    }
    printf("\n");
}

// See Algorithm 2 at http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf
set<int> luby_algorithm(int procID, int nproc, int n, int E, set<int> * adj_list) {
    const int root = 0;

    int N = (1 << n);   // number of vertices
    set<int> active_set;
    int * active_set_array = (int *) calloc(N, sizeof(int));
    for (int i = 0; i < N; i++) {
        active_set.insert(i);
        active_set_array[i] = i;
    }

    // Timing
    double time_inc_start;
    double time_inc;
    // 1. Assignment of random priorities
    double total_time_step1 = 0;
    // 2. Communicate priorities to all processes
    double total_time_step2 = 0;
    // 3. Computing max priority out of all neighbors.
    double total_time_step3 = 0;
    // 4. Determine which vertices should be added to independent set and send to root
    double total_time_step4 = 0;
    // 5. Communicate newly added vertices to all processes.
    double total_time_step5 = 0;
    // 6. Communicate active set to all processes
    double total_time_step6 = 0;
    // 7. Broadcasts final independent set to all processes
    double total_time_step7 = 0;

    set<int> maximal_independent_set;
    int * maximal_independent_set_array;
    while (active_set.size() > 0) {

        // 1. Assign random priorities to all of the vertices in the active set.
        time_inc_start = MPI_Wtime();
        int nc = pow(N, PROB_POWER);
        float * random = (float *) calloc(N, sizeof(float));
        active_set_array = (int *) calloc(active_set.size(), sizeof(int));
        if (procID == root) {
            for (int i = 0; i < N; i++) {
                if (active_set.count(i)) random[i] = rand() % nc + 1;
            }

            std::set<int, greater<int>>::iterator itr;
            int idx = 0;
            for (itr = active_set.begin(); itr != active_set.end(); itr++) {
                active_set_array[idx] = *itr;
                idx++;
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step1 += time_inc;

        if (DEBUG) printf("PROCID: %d After 1\n", procID);
        

        // 2. Communicate the random priorities (and vertex-index pairs) to all
        // the vertices.
        time_inc_start = MPI_Wtime();
        MPI_Bcast(random, N, MPI_FLOAT, root, MPI_COMM_WORLD);
        MPI_Bcast(active_set_array, active_set.size(), MPI_INT, root, MPI_COMM_WORLD);
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step2 += time_inc;
        
        if (DEBUG) printf("PROCID: %d After 2\n", procID);
        

        // 3. For each vertex in the active set, compare its priority to that
        // of its neighbors and add it to the maximal independent set if its
        // priority is bigger than those of its neighbors.
        time_inc_start = MPI_Wtime();
        int vertices_per_group;
        if (active_set.size() % nproc == 0) vertices_per_group = active_set.size() / nproc;
        else vertices_per_group = active_set.size() / nproc + 1;
        int start = vertices_per_group * procID;
        int active_set_size = active_set.size();
        int end = std::min(vertices_per_group * (procID + 1), active_set_size);
        for (int idx = start; idx < end; idx++) {
            int v = active_set_array[idx];
            bool flag = true;
            set<int> neighbors = set_intersect(adj_list[v], active_set);
            set<int, greater<int>>::iterator itr;
            for (itr = neighbors.begin(); itr != neighbors.end(); itr++) {
                // Check all u in N(v)
                if (random[v] <= random[*itr]) {
                    flag = false;
                    break;
                }
            }
            if (flag) maximal_independent_set.insert(v);
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step3 += time_inc;

        if (DEBUG) printf("PROCID %d After 3\n", procID);

        MPI_Barrier(MPI_COMM_WORLD);
        
        // 4. Processors combine their values of maximal_independent_set.
        // Process 0 will hold the new maximal independent set.
        time_inc_start = MPI_Wtime();
        int set_size;
        if (procID != root) {
            set_size = maximal_independent_set.size();
            maximal_independent_set_array = setToArray(maximal_independent_set);
            MPI_Send(&set_size, 1, MPI_INT, root, 0, MPI_COMM_WORLD);
            MPI_Send(maximal_independent_set_array, set_size, MPI_INT, root, 0, MPI_COMM_WORLD);
        }
        else {
            for (int i = root + 1; i < nproc; i++) {
                // Receive length first
                MPI_Recv(&set_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                maximal_independent_set_array = (int *) calloc(set_size, sizeof(int));

                // Receive set
                MPI_Recv(maximal_independent_set_array, set_size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                set<int> received_independent_set = arrayToSet(maximal_independent_set_array, set_size);

                // Take union
                set_union_yes(maximal_independent_set, received_independent_set);
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step4 += time_inc;

        if (DEBUG) printf("PROCID %d After 4\n", procID);

        // 5. Process 0 broadcasts its maximal independent set to all of the other processes.
        time_inc_start = MPI_Wtime();
        if (procID == root) set_size = maximal_independent_set.size();
        MPI_Bcast(&set_size, 1, MPI_INT, root, MPI_COMM_WORLD);
        int * maximal_independent_set_array_recv;
        if (procID == root) 
        {
            maximal_independent_set_array_recv = setToArray(maximal_independent_set);
        }
        else {
            maximal_independent_set_array_recv = (int *) calloc(set_size, sizeof(int));
        }
        MPI_Bcast(maximal_independent_set_array_recv, set_size, MPI_INT, root, MPI_COMM_WORLD);
        maximal_independent_set = arrayToSet(maximal_independent_set_array_recv, set_size);
        
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step5 += time_inc;

        if (DEBUG) printf("procID %d After 5\n", procID);

        // 6. Process 0 computes the new active set and broadcasts this as well.
        time_inc_start = MPI_Wtime();
        int * active_set_array;
        if (procID == root) {
            set<int> killed_vertices;
            set<int, greater<int>>::iterator itr;
            for (itr = maximal_independent_set.begin(); itr != maximal_independent_set.end(); itr++) {
                set<int, greater<int>>::iterator itrTwo;
                for (itrTwo = adj_list[*itr].begin(); itrTwo != adj_list[*itr].end(); itrTwo++) {
                    killed_vertices.insert(*itrTwo);
                }
                killed_vertices.insert(*itr);
            }
            
            // if (DEBUG) printf("Current active set: ");
            // if (DEBUG) print_set(active_set);

            // Remove killed_vertices from active set
            active_set = set_minus(active_set, killed_vertices);
            active_set_size = active_set.size();
            active_set_array = setToArray(active_set);
        }
        MPI_Bcast(&active_set_size, 1, MPI_INT, root, MPI_COMM_WORLD);

        if (procID != root) {
            active_set_array = (int*) calloc(active_set_size, sizeof(int));
        }
        MPI_Bcast(active_set_array, active_set_size, MPI_INT, root, MPI_COMM_WORLD);
        if (procID != root) active_set = arrayToSet(active_set_array, active_set_size);
        
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step6 += time_inc;

        if (DEBUG) printf("PROCID %d After 6\n", procID);
        
    }

    // Broadcast final maximal independent set to all processes
    time_inc_start = MPI_Wtime();
    int output_size;
    if (procID == root) {
        output_size = maximal_independent_set.size();
        maximal_independent_set_array = setToArray(maximal_independent_set);
    }

    MPI_Bcast(&output_size, 1, MPI_INT, root, MPI_COMM_WORLD);
    if (procID != root) {
        maximal_independent_set_array = (int*) calloc(output_size, sizeof(int));
    }
    MPI_Bcast(maximal_independent_set_array, output_size, MPI_INT, root, MPI_COMM_WORLD);
    if (procID != root) maximal_independent_set = arrayToSet(maximal_independent_set_array, output_size);
    
    time_inc = MPI_Wtime() - time_inc_start;
    total_time_step7 += time_inc;

    // Running time
    if (procID == root) {
        printf("Process: %d, time for step %d: %f\n", procID, 1, total_time_step1);
        printf("Process: %d, time for step %d: %f\n", procID, 2, total_time_step2);
        printf("Process: %d, time for step %d: %f\n", procID, 3, total_time_step3);
        printf("Process: %d, time for step %d: %f\n", procID, 4, total_time_step4);
        printf("Process: %d, time for step %d: %f\n", procID, 5, total_time_step5);
        printf("Process: %d, time for step %d: %f\n", procID, 6, total_time_step6);
        printf("Process: %d, time for step %d: %f\n", procID, 7, total_time_step7);
    }
    
    return maximal_independent_set;
}