/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

// #include "wireroute.h"

#include <assert.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <stdio.h>
#include <vector>
#include <random>

static int _argc;
static const char **_argv;

const char *get_option_string(const char *option_name, const char *default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return _argv[i + 1];
    return default_value;
}

int get_option_int(const char *option_name, int default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return atoi(_argv[i + 1]);
    return default_value;
}

float get_option_float(const char *option_name, float default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return (float)atof(_argv[i + 1]);
    return default_value;
}


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
                                    std::uniform_real_distribution<double>& rng, std::mt19937& mersenne_twister) {

    if (i_end == i_start + 1 && j_end == j_start + 1) {
        return std::make_pair(i_start, j_start);
    }

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
std::vector<int> * generate_rmat_graph(int n, int E, double a, double b, double d) {

    std::random_device seed;
    std::mt19937 mersenne_twister(seed());
    std::uniform_real_distribution<double> rng(0.0, 1.0);
    double c = b; //probability of the lower left square of adjacency matrix.

    int N = (1 << n); // number of vertices
    std::vector<int> * adjacency_list = (std::vector<int> *) calloc(N, sizeof(std::vector<int>));
    int num_edges_left = E;
    while (num_edges_left > 0) {

        std::pair<int, int> new_edge = generate_edge(0, N, 0, N, a, b, c, d, rng, mersenne_twister);
        int u = new_edge.first;
        int v = new_edge.second;
        if (u < v) continue; // Ignore all the edges that are above the main diagonal.

        // Check if this is a new edge
        bool edge_seen = false;
        for (uint i = 0; i < adjacency_list[u].size(); i++) {
            int neighbor = adjacency_list[u][i];
            if (neighbor == v) {
                edge_seen = true;
                break;
            }
        }

        // Add if this is edge has not been seen
        if (!edge_seen) {
            adjacency_list[u].push_back(v);
            adjacency_list[v].push_back(u);
            num_edges_left--;
        }
    }

    return adjacency_list;
}

int main(int argc, const char *argv[]) {
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;

    auto init_start = Clock::now();
    double init_time = 0;

    _argc = argc - 1;
    _argv = argv + 1;

    int n = get_option_int("-n", 10);
    int E = get_option_int("-E", 100);
    double a = get_option_float("-a", 0.1f);
    double b = get_option_float("-b", 0.3f);
    double d = get_option_float("-d", 0.3f);

    printf("Number of Nodes: %d Number of Edges: %d\n", (1 << n), E);
    printf("Probability Params: %lf %lf %lf.\n", a, b, d);

    // 1. Generate random graph (adjacency list format) using 
    // R-MAT (random graph model due to Chakrabarti, Zhan and 
    // Faloutsos, and used in the work of Blelloch, Fineman and 
    // Shun)
    // Description of R-MAT method available at:
    // https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf
    std::vector<int> * adjacency_list = generate_rmat_graph(n, E, a, b, d);

    printf("done\n");
    
    int N = (1 << n); // number of vertices
    for (int u = 0; u < N; u++) {
        printf("Neighbors of %d: ", u);
        for (uint i = 0; i < adjacency_list[u].size(); i++) {
            printf("%d ", adjacency_list[u][i]);
        }
        printf("\n");
    }

    init_time += duration_cast<dsec>(Clock::now() - init_start).count();
    printf("Initialization Time: %lf.\n", init_time);

    auto compute_start = Clock::now();
    double compute_time = 0;


    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("Computation Time: %lf.\n", compute_time);


    return 0;
}