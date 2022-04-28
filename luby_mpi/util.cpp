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
#include "util.h"

void print_boolean_array(bool * array, int length, int start) {
    for (int i = 0; i < length; i++) {
        if (array[i]) printf("%d " , i + start);
    }
    printf("\n");
}