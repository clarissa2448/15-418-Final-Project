#include "mpi.h"
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
#include <chrono>
#include <sstream>
#include <string>

#define BUFFER_LENGTH 300000
#define DEBUG true
#define PROB_POWER 3

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

void print_set(set<int> A) {
    set<int, greater<int>>::iterator itr;
    for (itr = A.begin(); itr != A.end(); itr++) {
        printf("%d ", *itr);
    }
    printf("\n");
}


// Write Input Adjacency List to File
void write_adj_list_to_file(set<int>* adj_list, int n, int E, int nproc) {
    ofstream adj_list_file;
    int N = 1 << n;
    char buffer[BUFFER_LENGTH];
    snprintf(buffer, BUFFER_LENGTH, "outputs/adj_list_%d_%d.txt", n, E);
    adj_list_file.open(buffer);
    for (int i = 0; i < N; i++) {
        for (auto j = adj_list[i].begin(); j != adj_list[i].end(); j++) {
            adj_list_file << *j << " ";
        }
        adj_list_file << "\n";

        
    }
    adj_list_file.close();
}

// Read Adj from file
set<int>* read_adj_list_from_file(int n, int E, string filename) {
    int N = 1 << n;
    set<int> * adj_list = (std::set<int> *) calloc(N, sizeof(set<int>));
    std::ifstream file(filename);
    for (int i = 0; i < N; i++) {
        set<int> l;
        string line, tmp;
        std::getline(file, line);
        stringstream ss(line);
        while(getline(ss, tmp, ' ')){
            l.insert(stoi(tmp));
        }
        adj_list[i] = l;

    }
    return adj_list;

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
        if (u < v) continue; // Ignore all the edges that are above the main diagonal.

        // Check if this is a new edge
        if (adjacency_list[u].count(v)) continue;

        adjacency_list[u].insert(v);
        adjacency_list[v].insert(u);
        num_edges_left--;
    }
    return adjacency_list;


}

typedef std::pair<int, int> edge;

void print_boolean_array(bool * array, int length, int start) {
    for (int i = 0; i < length; i++) {
        if (array[i]) printf("%d " , i + start);
    }
    printf("\n");
}




int main(int argc, char *argv[]) {
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;
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
    auto start = Clock::now();
    set<int> * adj_list = generate_rmat_graph(n, E, a, b, d);
  
    printf("done gen graph\n");

    
    double end = duration_cast<dsec>(Clock::now() - start).count();

    write_adj_list_to_file(adj_list, n, E, 1);

    if (DEBUG) {
        set<int>* adj_list_again = read_adj_list_from_file(n, E, "outputs/adj_list_4_100.txt");
        int N = 1 << n;
        for (int i = 0; i < N; i ++) {
            set<int> l = adj_list_again[i];
            printf("%d: ", i);
            for (auto j = l.begin(); j != l.end(); j++) {
                printf("%d ", *j);
            }
            printf("\n");
        }
    }
    
    printf("Elapsed time: %f\n", end);
    return 0;
}
