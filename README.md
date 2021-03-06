# Summary

We will implement parallel algorithms for finding maximal independent sets in the shared memory, message passing, and data-parallel models.

# Background

We will be considering the problem of finding a maximal independent set of a graph G = (V, E), that is, a subset S is a subset V of vertices such that S is an independent set, and moreover, for any v in V setminus S, S union {v} is not an independent set. As described in Algorithm 1 of [this survey by Shi, Wang, Shang](http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf), there is a simple sequential algorithm for this problem: we can iterate through the vertices of G and maintain an independent set S, adding new vertices v to S if v is not adjacent to any of the vertices in S.

In this project, we will implement existing parallel algorithms for maximal independent set. One such algorithm is Luby's algorithm (shown in Algorithm 2 of Shi, Wang and Shang's survey) which maintains an independent set S over the course of several rounds, as well as a set of vertices U is a subset of V setminus S of G which are guaranteed to not be adjacent to any vertices in S. In each round, the vertices in U are assigned (in parallel) a random priority, and the vertices in U which have a higher priority than all of their neighbors are added to S --- these vertices and their neighbors are then removed from U. It can be shown that O(log n) rounds suffice with high probability to process all the vertices (here n = |V|).

We also plan to implement the algorithms of [Blelloch, Fineman and Shun](https://arxiv.org/pdf/1202.3205.pdf). Blelloch, Fineman and Shun observe that if the vertices are processed according to a random permutation (fixed at the beginning of the algorithm) then the span of a parallel version of the greedy sequential algorithm is O(log^2 n) with high probability over the initial random permutation. This algorithm also effectively maintains an independent set S over the course of several rounds, and a set U is a subset of V setminus S consisting of vertices which are not adjacent to S. In each round, for all the vertices in U whose neighbors in G do not come before them in the initial random permutation, the algorithm will take those vertices and add them to S, and then remove those vertices and their neighbors from U (and then recurse on U). This algorithm can be shown to require at most O(log^2 n) rounds with high probability over the initial permutation. The work required by this algorithm is O(m log^2 n), but Blelloch, Fineman and Shun also propose algorithms with O(m) work.

# Challenge

Due to the iterative nature of these algorithms, one key challenge in achieving good speedup is reducing the overhead of synchronization. Minimizing communication in the message-passing and data-parallel models is also likely to be nontrivial.

## Resources

Conceptual resources we are using include the survey of Shi, Wang, Shang mentioned above, as well as the work of Blelloch, Fineman and Shun. We are likely to implement these algorithms in the message-passing and data-parallel models from scratch. For implementations of these algorithms in the shared memory model, we will potentially use [this implementation of Luby's Algorithm](https://github.com/ldhulipala/Maximal-Independent-Sets) or [the Ligra framework](https://github.com/jshun/ligra) as starter code. We will be making use of the PSC machines and NVIDIA GPUs.

# Goals and Deliverables

Overall, we hope to determine how to minimize communication and synchronization overhead --- we will be showing speedup graphs at the poster session. Our goals are as follows:

## 75% Goals

In case the work goes slowly, we would aim to implement Luby's algorithm and Algorithms 2 and 3 in the paper by Blelloch, Fineman and Shun in the shared memory model (openMP).

## 100% Goals

Implement both algorithms in the message passing model (MPI).

## 125% Goals

If we make rapid progress on the above goals, we will empirically evaluate the following heuristic in the message-passing model: in each round, we will simply sort the vertices by degree and then remove vertices with low degree. This heuristic could potentially have the benefit of low communication, and may result in larger independent sets compared to the algorithms mentioned above which make progress by including vertices with large degrees.

# Schedule

- During the week of 3/27, we will aim to get a working shared memory implementation of Luby's algorithm.
- During the week of 4/3, we will aim to get a working shared memory implementation of the algorithms of Blelloch, Fineman and Shun.
- During the week of 4/10, we will aim to get a working message passing implementation of Luby's algorithm.
- During the week of 4/17, we will aim to implement the algorithms of Blelloch, Fineman and Shun in the message passing model.
- During the week of 4/24, we will attempt to implement our proposed heuristic in the message passing model.

# References
- Shi, Wang and Shang's survey [http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf](http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf)
- Guy E. Blelloch, Jeremy T. Fineman, Julian Shun. Greedy Sequential Maximal Independent Set and Matching are Parallel on Average. Published in SPAA 2012. [https://arxiv.org/pdf/1202.3205.pdf](https://arxiv.org/pdf/1202.3205.pdf)

# Milestone

## Work completed so far

We have decided to implement Luby's algorithm, and apply optimizations which may lead to good speedup. We chose Luby's algorithm since it is a classical parallel algorithm for this problem. There exist later algorithms such as that of Blelloch, Fineman, and Shun, but we choose Luby's algorithm due to its conceptual simplicity, as it has less parameters.

We used the [R-MAT](https://www.cs.cmu.edu/~christos/PUBLICATIONS/siam04.pdf) random graph model of Chakrabarti, Faloutsos and Zhan as our input dataset. We also have a preliminary implementation of Luby's algorithm in MPI. Our implementation is as follows:
- First, process 0 assigns random priorities to all vertices sequentially.
- Then, the entire array of priorities is communicated to all processes.
- Each process deals with a subset of the vertices in the active set (i.e. those which can still be added to the independent set). If these have higher priorities than their neighbors, then they are added to the process's independent set.
- The processes then combine their independent sets.
- Process 0 then computes the new active set and sends it to the rest.

## Progress on goals and deliverables, and updated schedule

We have a working MPI implementation of Luby's algorithm (which was originally our goal for the current week), though our implementation can be made more efficient. In particular, the assignment of random priorities to active vertices should be parallelized. In addition, instead of broadcasting the random priorities of all active vertices to all the processes, we can have communication only when it is necessary. In other words, we will partition the vertices between the processes, such that if a process is responsible for a vertex v that is in the active set, then this process will only communicate the priority of v to the processes which are responsible for the neighbors of v. The process which is responsible for v will also be the one to assign a priority to v.

We will proceed according to the following schedule:

### Second half of the week of April 11th:
- Implement a sequential greedy algorithm for maximal independent set (Arvind)
- Improve MPI implementation of Luby's algorithm so that in each iteration, only the new vertices added to the independent sets in that round are communicated (Clarissa)

### First half of the week of April 18th:
- Modify MPI implementation of Luby's algorithm so that each process is responsible for a contiguous subset of the vertices from 1 to N. This strategy may lead to better cache behavior. (Arvind)
- Read papers on graph partitioning to select an additional strategy (Clarissa and Arvind)

### Second half of the week of April 18th:
- Modify MPI implementation of Luby's algorithm so that after each iteration, the vertices in the active set are divided evenly between processes. The strategy may lead to more balanced workload. (Clarissa)
- Read papers on graph partitioning to select additional strategy (Clarissa and Arvind)

### Week of April 25th:
- Select one other known strategy for graph partitioning and modify the MPI implementation to use it (Clarissa and Arvind)

### MPI Results Midway
| n | E  | time  |
| :---:   | :-: | :-: |
| 10 | 2048 | 0.01 |
| 10 | 4096 | 0.01 |
| 10 | 131072 | 0.075 |
| 10 | 251072 | 0.201 |
| 10 | 520000 | 1.27 |
| 12 | 2048 | 0.01 |
| 12 | 130000 | 0.32 |
| 12 | 160000 | 0.37 |
| 12 | 240000 | 0.44 |
| 12 | 300000 | 0.51 |
| 12 | 500000 | 0.66 |
| 12 | 750000 | 0.84 |
| 12 | 1000000 | 1.13 | 
| 12 | 2000000 | 1.41 | 
| 12 | 4000000 | 2.65 | 
| 14 | 100000 | 2.251654| 
| 14 | 500000 | 3.66| 
| 14 | 1000000 | 5| 
| 14 | 2000000 | 6.47| 
| 14 | 4000000 | 8.43| 
| 14 | 10000000 | 12.4| 
