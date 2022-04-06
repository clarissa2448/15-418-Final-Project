#include "mpi.h"
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
#include <vector>
#include <set>

using namespace std;

#ifndef _LUBY_MPI_H
#define _LUBY_MPI_H

set<int> luby_algorithm(int procID, int nproc, int n, int E, set<int> * adj_list);
set<int> * generate_rmat_graph(int n, int E, double a, double b, double d);

#endif // _LUBY_MPI_H