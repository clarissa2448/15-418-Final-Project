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
- During the week of 4/24, we will attempt to implement our proposed heuristic in the shared-memory model.

# References
- Shi, Wang and Shang's survey [http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf](http://web.mit.edu/jeshi/www/public/papers/parallel_MIS_survey.pdf)
- Guy E. Blelloch, Jeremy T. Fineman, Julian Shun. Greedy Sequential Maximal Independent Set and Matching are Parallel on Average. Published in SPAA 2012. [https://arxiv.org/pdf/1202.3205.pdf](https://arxiv.org/pdf/1202.3205.pdf)
