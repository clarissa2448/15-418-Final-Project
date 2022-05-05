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

#define BUFFER_LENGTH 300000
#define DEBUG false
#define PROB_POWER 3

using namespace std;

static int _argc;
static const char **_argv;

const char *get_option_string(const char *option_name, const char *default_value, int _argc, char **_argv) {
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
void write_adj_list_to_file(set<int>* adj_list, int n, int E) {
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
void write_mis_to_file(set<int> mis, int n, int E) {
    ofstream mis_file;
    char buffer[BUFFER_LENGTH];
    snprintf(buffer, BUFFER_LENGTH, "outputs/maximal_indep_set_%d_%d.txt", n, E);
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


set<int> luby_algorithm(int n, int E, set<int> * adj_list) {

    // 0. Initialize Variables
    int N = (1 << n);
    int NC = pow(N, PROB_POWER);
    int start = 0;
    int end = N;
    bool * active_set = (bool *) calloc(end - start, sizeof(bool));
    float * random = (float *) calloc(end - start, sizeof(float));
    bool * independent_set = (bool *) calloc(end - start, sizeof(bool));
    set<int> * active_set_neighbors = (set<int> *) calloc(end - start, sizeof(set<int>));
    for (int i = start; i < end; i++) {
        set<int> s(adj_list[i]);
        active_set_neighbors[i - start] = s;
        active_set[i - start] = true;
        random[i - start] = rand() % NC + 1;
    }
    float * max_priority_of_neighbors = (float *) calloc(end - start, sizeof(float));
    for (int i = start; i < end; i++) max_priority_of_neighbors[i - start] = -1;

    // removed_neighbors_from_active_set, in each round, contains the vertices in [start, end) which
    // are being removed from the active set, and are neighbors of some vertex in process i.
    // vertex_in_other_process contains the neighbors of the vertices in
    // removed_neighbors_from_active_set.
    // std::vector<int> removed_neighbors_from_active_set; 
    // std::vector<int> vertex_in_other_process;

    // Timing
    // 1. Count vertices in active set.
    double total_time_step1 = 0;
    // 2. Assignment of random priorities.
    double total_time_step2 = 0;
    // 3. Computing max priority out of all neighbors.
    double total_time_step3 = 0;
    // 4. Determine which vertices should be added to independent set.
    double total_time_step4 = 0;
    // 5. Update active_set
    double total_time_step5 = 0;
    // 6. Update active_set_neighbors.
    double total_time_step6 = 0;

    while (true) {

        // 1. First determine the number of vertices in the active
        // set across all of the processes.
        double time_inc_start = MPI_Wtime();
        int local_active_set_count = 0;
        for (int i = start; i < end; i++) {
            if (active_set[i - start]) local_active_set_count++;
        }
        if (local_active_set_count == 0) break;
        double time_inc = MPI_Wtime() - time_inc_start;
        total_time_step1 += time_inc;

        // 2. Assign random priorities to all vertices in the active set.
        time_inc_start = MPI_Wtime();
        for (int i = start; i < end; i++) {
            if (active_set[i - start]) random[i - start] = rand() % NC + 1;
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step2 += time_inc;

        // 3. Iterate through the edges in edge_list and update
        // max_priority_of_neighbors.
        time_inc_start = MPI_Wtime();
        for (int u = start; u < end; u ++) {
            if (!active_set[u - start]) continue;
            set<int> u_neighbors = active_set_neighbors[u - start];
            set<int, greater<int>>::iterator itr;
            for (itr = u_neighbors.begin(); itr != u_neighbors.end(); itr++) {
                int v = *itr;
                max_priority_of_neighbors[u - start] = std::max(max_priority_of_neighbors[u - start], random[v - start]);
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step3 += time_inc;

        // 4. Determine which vertices are going to be newly
        // added to the independent set.
        time_inc_start = MPI_Wtime();
        std::vector<int> newly_added_to_independent_set;
        for (int u = start; u < end; u++) {
            if (!active_set[u - start]) continue;
            if (max_priority_of_neighbors[u - start] >= random[u - start]) continue;

            // Adding to independent set
            independent_set[u - start] = true;
            active_set[u - start] = false;
            newly_added_to_independent_set.push_back(u);
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step4 += time_inc;

        // 5. Update the active set
        // Go through the larger independent set 
        // and go throughs all vertices newly added to independent set 
        // and sees if it's a neighbor of u
        time_inc_start = MPI_Wtime();
        std::vector<int> newly_removed_from_active_set;
        int newly_added_to_independent_set_size = newly_added_to_independent_set.size();
        for (int u = start; u < end; u ++) {
            for (int i = 0; i < newly_added_to_independent_set_size; i ++) {
                if (active_set_neighbors[u - start].count(newly_added_to_independent_set[i]) > 0) {
                    active_set[u - start] = false;
                    newly_removed_from_active_set.push_back(u);
                    break;
                }
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step5 += time_inc;

        // 6. Update active_set_neighbors. First determine
        // vertices in this process which are removed from
        // the independent set, then for each other process i,
        // send the removed vertices which are neighbors of i.
        time_inc_start = MPI_Wtime();
        for (uint i = 0; i < newly_removed_from_active_set.size(); i++) {
            int u = newly_removed_from_active_set[i];
            std::set<int, greater<int>>::iterator itr;
            for (itr = active_set_neighbors[u - start].begin(); itr != active_set_neighbors[u - start].end(); itr++) {
                int v = *itr;
                active_set_neighbors[v - start].erase(u);
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step6 += time_inc;

    }
    // Running time
    printf("Time for step %d: %f\n", 1, total_time_step1);
    printf("Time for step %d: %f\n", 2, total_time_step2);
    printf("Time for step %d: %f\n", 3, total_time_step3);
    printf("Time for step %d: %f\n", 4, total_time_step4);
    printf("Time for step %d: %f\n", 5, total_time_step5);
    printf("Time for step %d: %f\n", 6, total_time_step6);

    // Gather the independent set
    set<int> final_local_independent_set;
    for (int u = start; u < end; u++) {
        if (independent_set[u - start]) {
            final_local_independent_set.insert(u);
        }
    }
    return final_local_independent_set;



}




int main(int argc, char *argv[]) {
    double startTime;
    double endTime;

    // If the following option is 1, we create a graph and return
    // If the following option is 0, we read a graph from an input file
    int create_graph = get_option_int("-g", 1, argc, argv);
    int n = get_option_int("-n", 10, argc, argv);
    int N = (1 << n);
    set<int> * adj_list;
    // set<int> * adj_list = (std::set<int> *) calloc(N, sizeof(set<int>)); 
    int E = get_option_int("-E", 100, argc, argv);
    double a = get_option_float("-a", 0.1f);
    double b = get_option_float("-b", 0.3f);
    double d = get_option_float("-d", 0.3f);
    string input_filename = get_option_string("-f", "", argc, argv);


    printf("Number of Nodes: %d Number of Edges: %d\n", n, E);
    printf("Probability Params: %lf %lf %lf.\n", a, b, d);


    if (create_graph == 1) {
        // 1. Generate random graph (adjacency matrix format) using 
        // R-MAT (random graph model due to Chakrabarti, Zhan and 
        // Faloutsos, and used in the work of Blelloch, Fineman and 
        // Shun)
        // Description of R-MAT method available at:
        // https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf
        adj_list = generate_rmat_graph(n, E, a, b, d);

        if (DEBUG) {
            int num_edges = 0;
            for (int i = 0; i < N; i ++) {
                printf("%d: ", i);
                for (auto itr = adj_list[i].begin(); itr != adj_list[i].end(); itr ++) {
                    printf("%d ", *itr);
                }
                num_edges += adj_list[i].size();
                printf("\n");
            }
            printf("NUMBER EDGES: %d \n", num_edges);
        }

        write_adj_list_to_file(adj_list, n, E);
        printf("Graph created\n");
        return 0;
    } else {
        printf("reading\n");
         adj_list = read_adj_list_from_file(n, E, input_filename);
    }

    startTime = MPI_Wtime();
    set<int> M = luby_algorithm(n, E, adj_list);
    endTime = MPI_Wtime();
    // printf("Nodes in Maximal Independent Set\n");
    // set<int, greater<int> >::iterator itr;
    // for (itr = M.begin(); itr != M.end(); itr++) {
    //     printf("%d ", *itr);
    // }
    // printf("\n");
    write_mis_to_file(M, n, E);
    
    printf("Elapsed time: %f\n", endTime-startTime);
    return 0;
}
