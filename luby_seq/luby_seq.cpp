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

#define BUFFER_LENGTH 300000
#define DEBUG false
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

// Convert array to set
set<int> arrayToSet(int* arr, int size) {
    set<int> A;
    for (int i = 0; i < size; i++) {
        A.insert(arr[i]);
    }
    return A;
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

// See Algorithm 2 at http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf
set<int> luby_algorithm(int n, int E, set<int> * adj_list) {
    printf("luby got n: %d, E: %d \n", n, E);
    int N = (1 << n);   // number of vertices
    set<int> active_set;
    int * active_set_array = (int *) calloc(N, sizeof(int));
    for (int i = 0; i < N; i++) {
        active_set.insert(i);
        active_set_array[i] = i;
    }

    set<int> maximal_independent_set;
    // int * maximal_independent_set_array;
    while (active_set.size() > 0) {

        // 1. Assign random priorities to all of the vertices in the active set.
        int nc = pow(n, PROB_POWER);
        float * random = (float *) calloc(N, sizeof(float));
        active_set_array = (int *) calloc(active_set.size(), sizeof(int));
        for (int i = 0; i < N; i++) {
            if (active_set.count(i)) random[i] = rand() % nc + 1;
        }

        std::set<int, greater<int>>::iterator itr;
        int idx = 0;
        for (itr = active_set.begin(); itr != active_set.end(); itr++) {
            active_set_array[idx] = *itr;
            idx++;
        }

        // 2. For each vertex in the active set, compare its priority to that
        // of its neighbors and add it to the maximal independent set if its
        // priority is bigger than those of its neighbors.
        for (uint idx = 0; idx < active_set.size(); idx++) {
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

        // 3. Process 0 computes the new active set and broadcasts this as well.
        // int active_set_size;
        // int * active_set_array;
        set<int> killed_vertices;
        // set<int, greater<int>>::iterator itr;
        for (itr = maximal_independent_set.begin(); itr != maximal_independent_set.end(); itr++) {
            set<int, greater<int>>::iterator itrTwo;
            for (itrTwo = adj_list[*itr].begin(); itrTwo != adj_list[*itr].end(); itrTwo++) {
                killed_vertices.insert(*itrTwo);
            }
            killed_vertices.insert(*itr);
        }
        if (DEBUG) printf("Current active set: ");
        if (DEBUG) print_set(active_set);

        // Remove killed_vertices from active set
        active_set = set_minus(active_set, killed_vertices);
        // active_set_size = active_set.size();
        // active_set_array = setToArray(active_set);
    }

    // int output_size;
    // output_size = maximal_independent_set.size();
    // maximal_independent_set_array = setToArray(maximal_independent_set);


    return maximal_independent_set;
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
    auto start = Clock::now();
    set<int> M = luby_algorithm(n, E, adj_list);
    double end = duration_cast<dsec>(Clock::now() - start).count();

    printf("Nodes in Maximal Independent Set\n");
    set<int, greater<int> >::iterator itr;
    for (itr = M.begin(); itr != M.end(); itr++) {
        printf("%d ", *itr);
    }
    printf("\n");

    write_adj_list_to_file(adj_list, n, E, 1);
    write_mis_to_file(M, n, E, 1);
    
    printf("Elapsed time: %f\n", end);
    return 0;
}
