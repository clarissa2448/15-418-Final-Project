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
#include <vector>
#include <set>

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

// set<int> seq_algo(int n, int E, set<int> * adj_list) {
//     if 
// }

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
        if (u < v) continue; // Ignore all the edges that are above the main diagonal.

        // Check if this is a new edge
        if (adjacency_list[u].count(v)) continue;

        adjacency_list[u].insert(v);
        adjacency_list[v].insert(u);
        num_edges_left--;
    }
    return adjacency_list;
}

// https://arxiv.org/pdf/1202.3205.pdf Page 3
set<int> sequential_mis(int n, int E, set<int> * G, set<int> & avail) {
    if (!avail.size()) {
        return set<int>();
    }
    int v = *avail.begin();
    avail.erase(avail.begin());
    for(auto itr = G[v].begin(); itr != G[v].end(); itr++) {
        avail.erase(*itr);
    }
    set<int> res = sequential_mis(n, E, G, avail);
    res.insert(v);
    return res;
}

int main(int argc, char *argv[]) {
    // double startTime;
    // double endTime;

    int n = get_option_int("-n", 10, argc, argv);
    int E = get_option_int("-E", 100, argc, argv);
    double a = get_option_float("-a", 0.1f);
    double b = get_option_float("-b", 0.3f);
    double d = get_option_float("-d", 0.3f);

    printf("Number of Nodes: %d Number of Edges: %d\n", n, E);
    printf("Probability Params: %lf %lf %lf.\n", a, b, d);

    // 1. Generate random graph (adjacency matrix format) using 
    // R-MAT (random graph model due to Chakrabarti, Zhan and 
    // Faloutsos, and used in the work of Blelloch, Fineman and 
    // Shun)
    // Description of R-MAT method available at:
    // https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf
    set<int> * adj_list = generate_rmat_graph(n, E, a, b, d);
  
    printf("done gen graph\n");
    int N = (1 << n);
    for (int i = 0; i < N; i ++) {
        printf("%d: ", i);
        for (auto itr = adj_list[i].begin(); itr != adj_list[i].end(); itr ++) {
            printf("%d, ", *itr);
        }
        printf("\n");
    }

    set<int> avail;
    for (int i = 0; i < 1<<n; i++)
        avail.insert(avail.end(), i);
    set<int> res = sequential_mis(n, E, adj_list, avail);
    for (auto itr = res.begin(); itr != res.end(); itr++)
        printf("res vert %d\n", *itr);
    
    // printf("Elapsed time: %f\n", endTime - startTime);
    return 0;
}
