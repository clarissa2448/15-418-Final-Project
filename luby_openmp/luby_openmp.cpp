#include <omp.h>
#include <stdio.h>
#include <vector>
#include <random>

typedef std::pair<int, int> edge;

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
bool generate_edge(bool** adj_matrix, int i_start, int i_end, int j_start, int j_end, 
                    double a, double b, double c, double d, 
                    std::uniform_real_distribution<double> rng, std::mt19937 mersenne_twister) {
    if (i_end == i_start + 1 && j_end == j_start + 1) {
        if (adj_matrix[i_start][j_start]) return false;
        adj_matrix[i_start][j_start] = true;
        return true;
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
    return generate_edge(adj_matrix, i_start, i_end, j_start, j_end,
                            a, b, c, d, rng, mersenne_twister);
}



// n: Log (base 2) of the number of vertices
// E: number of edges in the generated graph
// a: Probability of edge going into upper left of adj. matrix
// b: Probability of edge going into upper right of adj. matrix
// c: Probability of edge going into lower left of adj. matrix
// d: Probability of edge going into lower right of adj. matrix
// Note that this generates adjacency matrix of a directed graph,
// but we can use this to generate the adjacency matrix of an
// undirected graph as described in Section 3.4 of the work of
// Chakrabarti, Zhan and Faloutsos.
bool** generate_rmat_graph(int n, int E, double a, double b, double c, double d) {

    std::random_device seed;
    std::mt19937 mersenne_twister(seed());
    std::uniform_real_distribution<double> rng(0.0, 1.0);

    int N = (1 << n); // number of vertices
    bool* adj_matrix_1d = (bool*) calloc(bool, N * N);
    bool** adj_matrix = (bool**) calloc(bool*, N);
    for (int i = 0; i < N; i++) adj_matrix[i] = (adj_matrix_1d + i * N);

    int num_edges_left = E;
    while (num_edges_left > 0) {
        bool new_edge = generate_edge(adj_matrix, 0, N, 0, N, a, b, c, d, rng, mersenne_twister);
        if (new_edge) num_edges_left--;
    }

    return adj_matrix;
}



int main() {

    int n = 10;
    int E = 200;
    int a = 0.1;
    int b = 0.3;
    int c = 0.3;
    int d = 0.3;

    // 1. Generate random graph (adjacency matrix format) using 
    // R-MAT (random graph model due to Chakrabarti, Zhan and 
    // Faloutsos, and used in the work of Blelloch, Fineman and 
    // Shun)
    // Description of R-MAT method available at:
    // https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf
    bool ** adj_matrix = generate_rmat_graph(n, E, a, b, c, d);

    return 0;
}