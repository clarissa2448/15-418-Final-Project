#include "luby_mpi_blocked_assignment.h"
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

typedef std::pair<int, int> edge;

void print_boolean_array(bool * array, int length, int start) {
    for (int i = 0; i < length; i++) {
        if (array[i]) printf("%d " , i + start);
    }
    printf("\n");
}

// See Algorithm 2 at http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf
set<int> luby_algorithm_blocked_assignment(int procID, int nproc, int n, int E, set<int> * adj_list) {

    const int root = 0;

    // 0. This process will calculate the block of vertices for 
    // which it is responsible. It will keep track of which
    // vertices in that block are in the active set and the
    // random priorities for each one. In addition, for each
    // vertex in this block, we keep track of which of the neighbors
    // are still in the active set.
    int N = (1 << n);
    int NC = pow(N, PROB_POWER);
    int vertices_per_process;
    if (N % nproc == 0) vertices_per_process = N / nproc;
    else vertices_per_process = N / nproc + 1;
    int start = vertices_per_process * procID;
    int end = std::min(vertices_per_process * (procID + 1), N);
    bool * active_set = (bool *) calloc(end - start, sizeof(bool));
    float * random = (float *) calloc(end - start, sizeof(float));
    bool * independent_set = (bool *) calloc(end - start, sizeof(float));
    set<int> * active_set_neighbors = (set<int> *) calloc(end - start, sizeof(set<int>));
    for (int i = start; i < end; i++) {
        set<int> s(adj_list[i]);
        active_set_neighbors[i - start] = s;
        active_set[i - start] = true;
        random[i - start] = rand() % NC + 1;
    }
    float * max_priority_of_neighbors = (float *) calloc(end - start, sizeof(float));
    for (int i = start; i < end; i++) max_priority_of_neighbors[i - start] = -1;

    int* larger_independent_set = (int *) calloc(N, sizeof(int));
    int* independent_set_sizes = (int *) calloc(nproc, sizeof(int));
    int* independent_set_displ = (int *) calloc(nproc, sizeof(int));

    // removed_neighbors_from_active_set[i], in each round, contains the vertices in [start, end) which
    // are being removed from the active set, and are neighbors of some vertex in process i.
    // vertex_in_other_process[i] contains the neighbors, in process i, of the vertices in
    // removed_neighbors_from_active_set[i].
    std::vector<int> * removed_neighbors_from_active_set = (std::vector<int> *) calloc(nproc, sizeof(std::vector<int>)); 
    std::vector<int> * vertex_in_other_process = (std::vector<int> *) calloc(nproc, sizeof(std::vector<int>));

    int* array_to_receive_neighbors = (int *) calloc(2 * E, sizeof(int));

    // Timing
    // 1. Count vertices in active set.
    double total_time_step1 = 0;
    // 2. Assignment of random priorities.
    double total_time_step2 = 0;
    // 3. Creating sorted list of edges.
    double total_time_step3 = 0;
    // 4. Computing max priority out of all neighbors.
    double total_time_step4 = 0;
    // 5. Determine which vertices should be added to independent set.
    double total_time_step5 = 0;
    // 6. Communicate newly added vertices to all processes.
    double total_time_step6 = 0;
    // 7. Update active_set.
    double total_time_step7 = 0;
    // 8. Update active_set_neighbors.
    double total_time_step8 = 0;

    while (true) {

        // 1. First determine the number of vertices in the active
        // set across all of the processes.
        double time_inc_start = MPI_Wtime();
        int local_active_set_count = 0;
        for (int i = start; i < end; i++) {
            if (active_set[i - start]) local_active_set_count++;
        }
        int overall_active_set_count;
        MPI_Reduce(&local_active_set_count, &overall_active_set_count, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Bcast(&overall_active_set_count, 1, MPI_INT, root, MPI_COMM_WORLD);
        if (overall_active_set_count == 0) break;
        double time_inc = MPI_Wtime() - time_inc_start;
        total_time_step1 += time_inc;

        if (DEBUG) {
            printf("active set: ");
            print_boolean_array(active_set, end - start, start);

            printf("independent set: ");
            print_boolean_array(independent_set, end - start, start);
        }

        // 2. Assign random priorities to all vertices in the active set.
        time_inc_start = MPI_Wtime();
        for (int i = start; i < end; i++) {
            if (active_set[i - start]) random[i - start] = rand() % NC + 1;
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step2 += time_inc;

        // 3. Create a list of edges (u, v) where u < v, sorted
        // first by u and then by v.
        time_inc_start = MPI_Wtime();
        std::vector<edge> edge_list;
        for (int u = start; u < end; u++) {
            if (!active_set[u - start]) continue;
            set<int> u_neighbors = active_set_neighbors[u - start];
            set<int, greater<int>>::iterator itr;
            for (itr = u_neighbors.begin(); itr != u_neighbors.end(); itr++) {
                int v = *itr;
                edge_list.push_back(std::make_pair(u, v));
            }
        }
        std::sort(edge_list.begin(), edge_list.end(), [](edge a, edge b){
            float minA = min(a.first, a.second); 
            float minB = min(b.first, b.second);
            float maxA = max(a.first, a.second); 
            float maxB = max(b.first, b.second);
            if (minA < minB) return true;
            else if (minA == minB) return maxA < maxB;
            else return false;
        });
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step3 += time_inc;

        // 4. Iterate through the edges in edge_list and update
        // max_priority_of_neighbors.
        time_inc_start = MPI_Wtime();
        for (int i = 0; i < edge_list.size(); i++) {
            edge e = edge_list[i];
            int u = e.first;
            int v = e.second;
            int u_process = u / vertices_per_process;
            int v_process = v / vertices_per_process;


            if (v_process == procID && u_process == procID) {
                if (DEBUG) printf("PROCID: %d, V process = procID | v_process: %d | u: %d | v: %d \n", procID, v_process, u, v);
                // Swap u and v - at this point, we can assume
                // that u belongs to procID.
                e = std::make_pair(v, u);
                u = e.first;
                v = e.second;
                u_process = u / vertices_per_process;
                v_process = v / vertices_per_process;
                max_priority_of_neighbors[u - start] = std::max(max_priority_of_neighbors[u - start], random[v - start]);
            }
            else if (u < v) {
                
                // First send u's priority to v_process, then
                // receive v's priority from v_process.
                float communicated_priority = random[u - start];
                if (DEBUG) printf("PROCID: %d, u < v | v_process: %d | u: %d v: %d, priority: %f \n", procID, v_process, u, v, communicated_priority);
                MPI_Send(&communicated_priority, 1, MPI_FLOAT, v_process, 0, MPI_COMM_WORLD);
                MPI_Recv(&communicated_priority, 1, MPI_FLOAT, v_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                max_priority_of_neighbors[u - start] = std::max(max_priority_of_neighbors[u - start], communicated_priority);
                if (DEBUG) printf("ProcID: %d, Received priority: %f\n", procID, communicated_priority);
            }
            else if (u > v) {
                if (DEBUG) printf("PROCID: %d u >= v | v_process: %d | u: %d | v: %d \n", procID, v_process, u, v);
                // First receive v's priority from v_process, then
                // send u's priority to v_process.
                float communicated_priority;
                MPI_Recv(&communicated_priority, 1, MPI_FLOAT, v_process, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (DEBUG) printf("procID: %d, Received priority: %f\n", procID, communicated_priority);
                max_priority_of_neighbors[u - start] = std::max(max_priority_of_neighbors[u - start], communicated_priority);
                communicated_priority = random[u - start];
                MPI_Send(&communicated_priority, 1, MPI_FLOAT, v_process, 0, MPI_COMM_WORLD);
                if (DEBUG) printf("procID: %d, Sent priority: %f\n", procID, communicated_priority);
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step4 += time_inc;

        // 5. Determine which vertices are going to be newly
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
        total_time_step5 += time_inc;

        // 6. Gather all the newly added vertices at root
        // and communicate those to all the processes.
        time_inc_start = MPI_Wtime();
        int newly_added_to_independent_set_size = newly_added_to_independent_set.size();
        int * newly_added_to_independent_set_array = (int *) calloc(newly_added_to_independent_set_size, sizeof(int));
        for (int i = 0; i < newly_added_to_independent_set_size; i++) {
            newly_added_to_independent_set_array[i] = newly_added_to_independent_set[i];
        }

        // Send the size of independent set to root 
        MPI_Gather(
            &newly_added_to_independent_set_size, 
            1, 
            MPI_INT, 
            independent_set_sizes, 
            1, 
            MPI_INT,
            root, 
            MPI_COMM_WORLD
        );
        int prefix_sum = 0;
        if (procID == root) {
            
            for (int i = 0; i < nproc; i ++) {
                independent_set_displ[i] = prefix_sum;
                prefix_sum += independent_set_sizes[i];
            }
        }

        MPI_Gatherv(
            newly_added_to_independent_set_array,
            newly_added_to_independent_set_size, 
            MPI_INT,
            larger_independent_set, 
            independent_set_sizes,
            independent_set_displ,
            MPI_INT,
            root, 
            MPI_COMM_WORLD
        );
        free(newly_added_to_independent_set_array);

        // Root sends independent set to everyone
        MPI_Bcast(&prefix_sum, 1, MPI_INT, root, MPI_COMM_WORLD);
        MPI_Bcast(larger_independent_set, prefix_sum, MPI_INT, root, MPI_COMM_WORLD);
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step6 += time_inc;

        // 7. Update the active set
        // Each process goes through the larger independent set 
        // and go throughs all vertices newly added to independent set 
        // and sees if it's a neighbor of u
        time_inc_start = MPI_Wtime();
        std::vector<int> newly_removed_from_active_set;
        for (int u = start; u < end; u ++) {
            for (int i = 0; i < prefix_sum; i ++) {
                if (active_set_neighbors[u - start].count(larger_independent_set[i]) > 0) {
                    active_set[u - start] = false;
                    newly_removed_from_active_set.push_back(u);
                    break;
                }
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step7 += time_inc;
        
        // 8. Update active_set_neighbors. First determine
        // vertices in this process which are removed from
        // the independent set, then for each other process i,
        // send the removed vertices which are neighbors of i.
        time_inc_start = MPI_Wtime();
        for (int i = 0; i < nproc; i++) {
            removed_neighbors_from_active_set[i].clear();
            vertex_in_other_process[i].clear();
        }
        for (int i = 0; i < newly_removed_from_active_set.size(); i++) {
            int u = newly_removed_from_active_set[i];
            std::set<int, greater<int>>::iterator itr;
            for (itr = active_set_neighbors[u - start].begin(); itr != active_set_neighbors[u - start].end(); itr++) {
                int v = *itr;
                int v_process = v / vertices_per_process;
                removed_neighbors_from_active_set[v_process].push_back(u);
                vertex_in_other_process[v_process].push_back(v);
            }
        }

        for (int i = 0; i < nproc; i++) {
            // Communicate with process i, the neighbors of i which are
            // being removed from the active set.
            if (i == procID) {
                for (int j = 0; j < removed_neighbors_from_active_set[procID].size(); j++) {
                    int removed_vertex = removed_neighbors_from_active_set[procID][j];
                    int neighbor = vertex_in_other_process[procID][j];
                    active_set_neighbors[neighbor - start].erase(removed_vertex);
                }
            }
            else if (i < procID) {
                // First receive the vertices from process i which are
                // neighbors of this process. array_to_receive_neighbors
                // contains a vertex u in process i, followed by its
                // neighbors in process procID.
                int num_vertices_communicated;
                MPI_Recv(&num_vertices_communicated, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(array_to_receive_neighbors, 2 * num_vertices_communicated, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Remove all of those neighbors from the adjacency sets.
                for (int j = 0; j < num_vertices_communicated; j++) {
                    int vertex_from_i = array_to_receive_neighbors[2 * j];
                    int neighbor_in_this_process = array_to_receive_neighbors[2 * j + 1];
                    active_set_neighbors[neighbor_in_this_process - start].erase(vertex_from_i);
                }

                // Now send the vertices from procID which are neighbors
                // of process i, but are being removed from active set.
                num_vertices_communicated = removed_neighbors_from_active_set[i].size();
                for (int j = 0; j < num_vertices_communicated; j++) {
                    int u = removed_neighbors_from_active_set[i][j];
                    int v = vertex_in_other_process[i][j];
                    array_to_receive_neighbors[2 * j] = u;
                    array_to_receive_neighbors[2 * j + 1] = v;
                }
                MPI_Send(&num_vertices_communicated, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(array_to_receive_neighbors, 2 * num_vertices_communicated, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            else {
                // Do the above but in reverse order.

                // Send the vertices from procID which are neighbors
                // of process i, but are being removed from the active set.
                int num_vertices_communicated = removed_neighbors_from_active_set[i].size();
                for (int j = 0; j < num_vertices_communicated; j++) {
                    int u = removed_neighbors_from_active_set[i][j];
                    int v = vertex_in_other_process[i][j];
                    array_to_receive_neighbors[2 * j] = u;
                    array_to_receive_neighbors[2 * j + 1] = v;
                }
                MPI_Send(&num_vertices_communicated, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(array_to_receive_neighbors, 2 * num_vertices_communicated, MPI_INT, i, 0, MPI_COMM_WORLD);

                // Then receive the vertices from process i which
                // are neighbors of procID and are removed from active set.
                MPI_Recv(&num_vertices_communicated, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(array_to_receive_neighbors, 2 * num_vertices_communicated, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // Remove those neighbors from the active set.
                for (int j = 0; j < num_vertices_communicated; j++) {
                    int vertex_from_i = array_to_receive_neighbors[2 * j];
                    int neighbor_in_this_process = array_to_receive_neighbors[2 * j + 1];
                    active_set_neighbors[neighbor_in_this_process - start].erase(vertex_from_i);
                }
            }
        }
        time_inc = MPI_Wtime() - time_inc_start;
        total_time_step8 += time_inc;
    }

    // Running time
    printf("Process: %d, time for step %d: %f\n", procID, 1, total_time_step1);
    printf("Process: %d, time for step %d: %f\n", procID, 2, total_time_step2);
    printf("Process: %d, time for step %d: %f\n", procID, 3, total_time_step3);
    printf("Process: %d, time for step %d: %f\n", procID, 4, total_time_step4);
    printf("Process: %d, time for step %d: %f\n", procID, 5, total_time_step5);
    printf("Process: %d, time for step %d: %f\n", procID, 6, total_time_step6);
    printf("Process: %d, time for step %d: %f\n", procID, 7, total_time_step7);
    printf("Process: %d, time for step %d: %f\n", procID, 8, total_time_step8);

    // Gather the independent set.
    int local_ind_set_size = 0;
    for (int i = start; i < end; i++) {
        if (independent_set[i - start]) local_ind_set_size++;
    }
    int overall_ind_set_size;
    MPI_Reduce(&local_ind_set_size, &overall_ind_set_size, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    MPI_Bcast(&overall_ind_set_size, 1, MPI_INT, root, MPI_COMM_WORLD);

    int * final_local_independent_set = (int *) calloc(local_ind_set_size, sizeof(int));
    int * final_overall_independent_set = (int *) calloc(overall_ind_set_size, sizeof(int));

    int local_idx = 0;
    for (int u = start; u < end; u++) {
        if (independent_set[u - start]) {
            final_local_independent_set[local_idx] = u;
            local_idx++;
        }
    }

    // Communicate sizes for allgatherv.
    int * local_ind_set_sizes = (int *) calloc(nproc, sizeof(int));
    int * local_ind_set_displs = (int *) calloc(nproc, sizeof(int));
    MPI_Allgather(&local_ind_set_size, 1, MPI_INT, local_ind_set_sizes, 1, MPI_INT, MPI_COMM_WORLD);
    int prefix_sum = 0;
    for (int i = 0; i < nproc; i++) {
        local_ind_set_displs[i] = prefix_sum;
        prefix_sum += local_ind_set_sizes[i];
    }

    MPI_Allgatherv(
        final_local_independent_set, 
        local_ind_set_size, 
        MPI_INT, 
        final_overall_independent_set, 
        local_ind_set_sizes, 
        local_ind_set_displs, 
        MPI_INT,
        MPI_COMM_WORLD
    );

    // Convert final overall independent set to a set.
    set<int> result;
    for (int i = 0; i < overall_ind_set_size; i++) {
        result.insert(final_overall_independent_set[i]);
    }
    return result;

}