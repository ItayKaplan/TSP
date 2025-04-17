"""
This is a comparison between the Nearest Neighbour algorithm and the Christofides
algorithm. Neither algorithm produces the exact solution; both deliver
approximate tours with far better time complexity than exact methods
(e.g., brute‑force has O(n!)) complexity).

- Nearest Neighbour is a greedy approximation algorithm with
  O(n^2) time complexity.
- Christofides is a metaheuristic combining MST construction,
  minimum‑weight perfect matching, and tour shortcutting, running in
  O(n^3) and guaranteeing at most a 1.5× approximation
  ratio under the triangle inequality. 

******************ESSENTIAL!******************

Before running this code, download the TSPLIB95 folder containing the
benchmark TSP instances and their known solutions from:
http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/

**********************************************

In this file, we benchmark and compare the solution quality and runtime
performance of both algorithms on various TSP instances. 
                                    
                        ( ͡° ͜ʖ ͡°)_/¯
"""

import os
import numpy as np
import tsplib95
from tqdm import tqdm
import networkx as nx
import re
from collections import defaultdict

# --- Nearest Neighbor Algorithm ---
def nearest_neighbour(distance_matrix):
    """Simple greedy implementation of the Nearest Neighbor heuristic for TSP."""
    n = len(distance_matrix)
    unvisited = list(range(1, n))
    route = [0]
    
    while unvisited:
        current = route[-1]
        nearest = min(unvisited, key=lambda x: distance_matrix[current][x])
        route.append(nearest)
        unvisited.remove(nearest)
    
    route.append(0)  # Return to start
    total_distance = calculate_tour_cost(route, distance_matrix)
    return route, total_distance

# --- Christofides Algorithm ---
def christofides_tsp(distance_matrix):
    """Solve TSP using Christofides algorithm (guaranteed 3/2-approximation for metric TSP)."""
    # Step 1: Minimum Spanning Tree (MST)
    mst = minimum_spanning_tree(distance_matrix)
    
    # Step 2: Find odd-degree vertices in MST
    odd_vertices = find_odd_degree_vertices(mst)
    
    # Step 3: Minimum-weight perfect matching on odd vertices
    min_matching = minimum_weight_matching(distance_matrix, odd_vertices)
    
    # Step 4-6: Combine MST and matching, find Eulerian circuit, create Hamiltonian tour
    multigraph = mst + min_matching
    eulerian_circuit = find_eulerian_circuit(multigraph)
    tour = create_hamiltonian_tour(eulerian_circuit)
    
    # Add return to starting node and calculate cost
    tour.append(tour[0])
    total_cost = calculate_tour_cost(tour, distance_matrix)
    
    return tour, total_cost

def minimum_spanning_tree(distance_matrix):
    """Create a Minimum Spanning Tree using Prim's algorithm."""
    n = len(distance_matrix)
    mst_edges = []
    visited = [False] * n
    
    # Start from the first city
    visited[0] = True
    
    while len(mst_edges) < n - 1:
        min_distance = float('inf')
        best_edge = None
        
        for current in range(n):
            if visited[current]:
                for next_city in range(n):
                    if not visited[next_city] and distance_matrix[current][next_city] < min_distance:
                        min_distance = distance_matrix[current][next_city]
                        best_edge = (current, next_city)
        
        if best_edge:
            mst_edges.append(best_edge)
            visited[best_edge[1]] = True
    
    return mst_edges

def find_odd_degree_vertices(mst_edges):
    """Find vertices with odd degree in the Minimum Spanning Tree."""
    degree = defaultdict(int)
    for u, v in mst_edges:
        degree[u] += 1
        degree[v] += 1
    
    return [vertex for vertex, deg in degree.items() if deg % 2 != 0]

def minimum_weight_matching(distance_matrix, odd_vertices):
    """Find minimum-weight perfect matching for odd-degree vertices."""
    if len(odd_vertices) <= 1:
        return []
    
    n_odd = len(odd_vertices)
    
    # Create cost matrix for the subgraph of odd vertices
    cost_matrix = np.zeros((n_odd, n_odd))
    for i in range(n_odd):
        for j in range(n_odd):
            if i != j:
                u, v = odd_vertices[i], odd_vertices[j]
                cost_matrix[i][j] = distance_matrix[u][v]
            else:
                cost_matrix[i][j] = float('inf')
    
    # Run Edmonds' algorithm
    match = edmonds_blossom(cost_matrix)
    
    # Convert matching to edges in the original graph
    matching_edges = []
    for i in range(n_odd):
        j = match[i]
        if i < j:  # Only add each edge once
            matching_edges.append((odd_vertices[i], odd_vertices[j]))
    
    return matching_edges

def edmonds_blossom(cost_matrix):
    """Optimized implementation of Edmonds' Blossom algorithm for minimum-weight perfect matching."""
    n = len(cost_matrix)
    
    # Initialize dual variables and matching
    dual_u = np.min(cost_matrix, axis=1)
    dual_v = np.zeros(n)
    match = [-1] * n
    
    # Continue until all vertices are matched
    while -1 in match:
        root = match.index(-1)
        
        # BFS tree initialization
        tree = [root]
        parent = [-1] * n
        
        # Initialize slacks
        slack = np.full(n, float('inf'))
        slack_u = np.zeros(n, dtype=int)
        
        for i in range(n):
            if i != root:
                slack[i] = cost_matrix[root][i] - dual_u[root] - dual_v[i]
                slack_u[i] = root
        
        # Find augmenting path
        while True:
            # Find vertex with minimum slack
            unassigned_vertices = [i for i in range(n) if parent[i] == -1 and i != root]
            if not unassigned_vertices:
                break
                
            j = min(unassigned_vertices, key=lambda i: slack[i])
            delta = slack[j]
            
            # Update dual variables
            for i in tree:
                dual_u[i] += delta
            for i in range(n):
                if parent[i] != -1:
                    dual_v[i] -= delta
                elif i != root:
                    slack[i] -= delta
            
            # Add j to the tree
            parent[j] = slack_u[j]
            if match[j] == -1:
                # Found augmenting path - perform path augmentation
                while j != -1:
                    k = parent[j]
                    next_matched = match[k]
                    match[k] = j
                    match[j] = k
                    j = next_matched
                break
            
            # Extend the tree
            k = match[j]
            tree.append(j)
            tree.append(k)
            parent[k] = j
            
            # Update slacks
            for i in range(n):
                if parent[i] == -1:
                    new_slack = cost_matrix[k][i] - dual_u[k] - dual_v[i]
                    if new_slack < slack[i]:
                        slack[i] = new_slack
                        slack_u[i] = k
    
    return match

def find_eulerian_circuit(edges):
    """Find an Eulerian circuit in the multigraph using Hierholzer's algorithm."""
    if not edges:
        return []
        
    # Build adjacency list
    graph = defaultdict(list)
    for u, v in edges:
        graph[u].append(v)
        graph[v].append(u)
    
    if not graph:
        return []
    
    # Hierholzer's algorithm
    circuit = []
    
    def dfs(vertex):
        while graph[vertex]:
            next_vertex = graph[vertex].pop()
            graph[next_vertex].remove(vertex)
            dfs(next_vertex)
        circuit.append(vertex)
    
    # Start from any vertex
    start_vertex = next(iter(graph))
    dfs(start_vertex)
    
    return list(reversed(circuit))

def create_hamiltonian_tour(eulerian_circuit):
    """Convert Eulerian circuit to Hamiltonian tour by removing duplicates."""
    tour = []
    visited = set()
    
    for vertex in eulerian_circuit:
        if vertex not in visited:
            tour.append(vertex)
            visited.add(vertex)
    
    return tour

def calculate_tour_cost(tour, distance_matrix):
    """Calculate the total cost of a TSP tour."""
    total_cost = 0
    for i in range(len(tour) - 1):
        total_cost += distance_matrix[tour[i]][tour[i+1]]
    return total_cost

def main():
    """Example usage of Christofides algorithm."""
    distances = np.array([
        [0, 10, 15, 20],
        [10, 0, 35, 25],
        [15, 35, 0, 30],
        [20, 25, 30, 0]
    ])
    
    tour, total_cost = christofides_tsp(distances)
    
    print("Approximate optimal tour:", tour)
    print("Total tour cost:", total_cost)

# --- Benchmarking Functions ---
def benchmark_tsp_algorithms(problem_file):
    """Benchmark NN and Christofides algorithms on a single TSP problem."""
    try:
        problem = tsplib95.load(problem_file)
        n = problem.dimension
        distance_matrix = np.zeros((n, n), dtype=float)
        
        # Create distance matrix from problem
        for i in range(n):
            for j in range(n):
                if i != j:
                    distance_matrix[i][j] = problem.get_weight(i+1, j+1)
        
        # Run algorithms
        nn_solution, nn_cost = nearest_neighbour(distance_matrix)
        christofides_solution, christofides_cost = christofides_tsp(distance_matrix)
        
        # Find optimal tour cost if available
        optimal_length = find_optimal_tour_length(problem_file, distance_matrix, problem)
        
        # Prepare results
        results = {
            'problem_name': os.path.basename(problem_file),
            'num_cities': n,
            'nn_cost': nn_cost,
            'christofides_cost': christofides_cost,
        }
        
        if optimal_length is not None:
            nn_deviation = (nn_cost - optimal_length) / optimal_length * 100
            christofides_deviation = (christofides_cost - optimal_length) / optimal_length * 100
            results.update({
                'optimal_length': optimal_length,
                'nn_deviation': nn_deviation,
                'christofides_deviation': christofides_deviation
            })
        
        return results
    
    except Exception as e:
        print(f"Error processing {os.path.basename(problem_file)}: {e}")
        return None

def find_optimal_tour_length(problem_file, distance_matrix, problem):
    """Extract optimal tour length from tour files or problem metadata."""
    # Try multiple possible tour file extensions
    basename = problem_file[:-4]  # Remove .tsp extension
    possible_extensions = ['.opt.tour', '.tour', '.opt']
    
    for ext in possible_extensions:
        tour_file = basename + ext
        if os.path.exists(tour_file):
            try:
                tour_problem = tsplib95.load(tour_file)
                if hasattr(tour_problem, 'tours') and tour_problem.tours:
                    optimal_tour = [i-1 for i in tour_problem.tours[0]]
                    if optimal_tour[0] != optimal_tour[-1]:
                        optimal_tour.append(optimal_tour[0])
                    return calculate_tour_cost(optimal_tour, distance_matrix)
            except Exception:
                pass
    
    # Check if optimal value is in comments
    if hasattr(problem, 'comment'):
        match = re.search(r'optimal[^\d]*(\d+)', problem.comment.lower())
        if match:
            return float(match.group(1))
    
    return None

def batch_benchmark_tsp_problems(path):
    """Benchmark multiple TSP problems and summarize results."""
    directory_path = path
        
    tsp_files = [
        os.path.join(root, file) 
        for root, _, files in os.walk(directory_path) 
        for file in files if file.lower().endswith('.tsp')
    ]
        
    all_results = []
    
    # Process each TSP file
    for file in tqdm(tsp_files, desc="Processing TSP Problems"):
        result = benchmark_tsp_algorithms(file)
        if result:
            all_results.append(result)
            print(f"\nProblem: {result['problem_name']}")
            print(f"NN Cost: {result['nn_cost']}")
            print(f"Christofides Cost: {result['christofides_cost']}")
            if 'optimal_length' in result:
                print(f"Optimal Length: {result['optimal_length']}")
                print(f"NN Deviation: {result['nn_deviation']:.2f}%")
                print(f"Christofides Deviation: {result['christofides_deviation']:.2f}%")
    
    # Display summary statistics
    if all_results:
        display_summary_statistics(all_results)

def display_summary_statistics(results):
    """Display summary statistics for all benchmark results."""
    nn_costs = [r['nn_cost'] for r in results]
    christofides_costs = [r['christofides_cost'] for r in results]
    
    print("\n--- General Summary ---")
    print(f"Average NN Cost: {np.mean(nn_costs):.2f}")
    print(f"Average Christofides Cost: {np.mean(christofides_costs):.2f}")
        
    # Print comprehensive summary table
    print("\n" + "="*80)
    print("{:<15} {:<10} {:<15} {:<15} {:<15}".format(
        "Graph Name", "Cities", "Optimal Cost", "NN Cost", "Christofides Cost"))
    print("="*80)
    
    for result in sorted(results, key=lambda x: x['num_cities']):
        optimal_cost = result.get('optimal_length', 'N/A')
        optimal_str = f"{optimal_cost:.2f}" if optimal_cost != 'N/A' else 'N/A'
        
        problem_name = result['problem_name']
        if problem_name.lower().endswith('.tsp'):
            problem_name = problem_name[:-4]
        
        print("{:<15} {:<10} {:<15} {:<15.2f} {:<15.2f}".format(
            problem_name[:15], 
            result['num_cities'],
            optimal_str,
            result['nn_cost'],
            result['christofides_cost']
        ))
    print("="*80)
    
    # Print deviations for instances with known optimal solutions
    optimal_available = [r for r in results if 'optimal_length' in r]
    if optimal_available:
        print("\n--- Performance Analysis (instances with known optimal solutions) ---")
        print("{:<15} {:<10} {:<15} {:<15}".format(
            "Graph Name", "Cities", "NN Dev %", "Christofides Dev %"))
        print("-"*60)
        
        for result in optimal_available:
            problem_name = result['problem_name']
            if problem_name.lower().endswith('.tsp'):
                problem_name = problem_name[:-4]
            
            print("{:<15} {:<10} {:<15.2f} {:<15.2f}".format(
                problem_name[:15],
                result['num_cities'],
                result['nn_deviation'],
                result['christofides_deviation']
            ))
            
        avg_nn_dev = np.mean([r['nn_deviation'] for r in optimal_available])
        avg_christofides_dev = np.mean([r['christofides_deviation'] for r in optimal_available])
        print("\nAverage NN Deviation: {:.2f}%".format(avg_nn_dev))
        print("Average Christofides Deviation: {:.2f}%".format(avg_christofides_dev))

if __name__ == "__main__":
    # Enter The path to the TSPLIB95 folder
    tsplib_path = r"ENTER THE PATH TO THE TSPLIB95 FOLDER HERE"
    batch_benchmark_tsp_problems(tsplib_path)
