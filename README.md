# Traveling Salesman Problem (TSP) Algorithms Collection

## Overview

A comprehensive collection of exact, approximation, and metaheuristic algorithms for solving the Traveling Salesman Problem (TSP). This repository provides implementations of:  
- Exact methods (brute‑force, branch‑and‑bound)  
- A simple greedy approximation (Nearest Neighbour)  
- A classic metaheuristic (Christofides)  
- A benchmarking script to compare Nearest Neighbour vs. Christofides 

Algorithm Comparison:

While exact methods guarantee an optimal tour, they suffer from factorial time complexity (O(n!)) and only handle very small instances. Approximation and metaheuristic methods trade optimality for speed:

  Nearest Neighbour:
        Greedy approximation
        Time complexity: O(n^2)
        No worst‑case approximation bound; typically yields tours within 20–50% of optimum on random instances.
        Pros: trivial to implement, very fast.
        Cons: early greedy choices may lead to long detours.

  Christofides:
        Metaheuristic for metric TSP
        Steps: MST -> perfect matching on odd‑degree vertices -> Eulerian tour shortcutting
        Time complexity: dominated by matching, O(n^3)
        Approximation guarantee: at most 1.5x the optimal tour length.
        Pros: provable quality bound, good practical performance.
        Cons: more complex to implement, higher constant factors.

Use the comparison script to evaluate solution quality and runtime across your own TSP instances (or on some TSPLIB95 TSP instances).
