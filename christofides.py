import numpy as np
import networkx as nx

def christofides_tsp(distance_matrix):

    # Minimum Spanning Tree (MST)
    mst = prim_mst(distance_matrix)
    
    # Find odd-degree vertices in MST
    odd_vertices = find_odd_degree_vertices(mst)
    
    # Minimum-weight perfect matching on odd vertices
    min_matching = minimum_weight_matching(distance_matrix, odd_vertices)
    
    # Combine MST and matching
    multigraph = combine_mst_and_matching(mst, min_matching)
    
    # Find Eulerian circuit
    eulerian_circuit = find_eulerian_path(multigraph)
    
    # Convert Eulerian circuit to Hamiltonian cycle
    tour = create_hamiltonian_tour(eulerian_circuit)
    
    # Calculate total tour cost
    total_cost = calculate_tour_cost(distance_matrix, tour)
    
    return tour, total_cost

def prim_mst(distance_matrix):
    
    n = len(distance_matrix)
    mst_edges = []
    visited = [False] * n
    
    # Start from the first city
    start_city = 0
    visited[start_city] = True
    
    while len(mst_edges) < n - 1:
        min_distance = float('inf')
        best_edge = None
        
        for current_city in range(n):
            if visited[current_city]:
                for next_city in range(n):
                    if not visited[next_city] and distance_matrix[current_city][next_city] < min_distance:
                        min_distance = distance_matrix[current_city][next_city]
                        best_edge = (current_city, next_city)
        
        if best_edge:
            mst_edges.append(best_edge)
            visited[best_edge[1]] = True
    
    return mst_edges

# Find vertices with odd degree in the given graph.
def find_odd_degree_vertices(mst_edges):
    # Determine the maximum vertex index to set list size
    max_vertex = 0
    for u, v in mst_edges:
        max_vertex = max(max_vertex, u, v)
    
    # Initialize degree list with zeros
    degree = [0] * (max_vertex + 1)
    
    # Count degrees
    for u, v in mst_edges:
        degree[u] += 1
        degree[v] += 1
    
    # Return vertices with odd degree
    return [i for i, deg in enumerate(degree) if deg % 2 != 0]

# Find minimum-weight perfect matching for vertices with odd degree.
# Uses NetworkX's implementation of Edmonds' Blossom algorithm.
def minimum_weight_matching(distance_matrix, odd_vertices):
    
    # Return empty matching if fewer than 2 odd vertices
    if len(odd_vertices) <= 1:
        return []
    
    # Create complete graph on odd vertices
    G = nx.Graph()
    
    # Add edges with weights from the distance matrix
    for i, u in enumerate(odd_vertices):
        for j, v in enumerate(odd_vertices):
            if i < j:  # Add each edge only once
                G.add_edge(u, v, weight=distance_matrix[u][v])
    
    # Find minimum weight matching using NetworkX implementation of Blossom algorithm
    matching = nx.algorithms.matching.min_weight_matching(G)
    
    # Convert to our edge format
    matching_edges = list(matching)
    
    return matching_edges

# Combine edges from the minimum spanning tree and the matching.
def combine_mst_and_matching(mst_edges, matching_edges):
    
    return mst_edges + matching_edges

# Find an Eulerian path in the given graph using a depth-first approach.
def find_eulerian_path(multigraph_edges):

    circuit = []
    graph = {}
    
    # Build adjacency list
    for u, v in multigraph_edges:
        if u not in graph:
            graph[u] = []
        if v not in graph:
            graph[v] = []
        graph[u].append(v)
        graph[v].append(u)
    
    def dfs(vertex):
        while graph[vertex]:
            next_vertex = graph[vertex].pop()
            graph[next_vertex].remove(vertex)
            dfs(next_vertex)
        circuit.append(vertex)
    
    # Start from any vertex
    start_vertex = list(graph.keys())[0]
    dfs(start_vertex)
    
    return list(reversed(circuit))

# Convert an Eulerian circuit to a Hamiltonian cycle by removing repeated vertices.
def create_hamiltonian_tour(eulerian_circuit):
    
    tour = []
    visited = set()
    
    for vertex in eulerian_circuit:
        if vertex not in visited:
            tour.append(vertex)
            visited.add(vertex)
    
    # Add return to starting node
    tour.append(tour[0])  
    return tour

# Calculate the total cost of a tour.
def calculate_tour_cost(distance_matrix, tour):

    total_cost = 0
    for i in range(len(tour) - 1):
        total_cost += distance_matrix[tour[i]][tour[i+1]]
    
    return total_cost

def main():
    
    # Example distance matrix (symmetric)
    distances = np.array([
        [0, 19, 21, 15, 12],
        [19, 0, 5, 11, 7],
        [21, 5, 0, 11, 20],
        [15, 11, 10, 0, 31],
        [12, 7, 20, 31, 0] 
    ])
    
    # Get tour and total cost
    tour, total_cost = christofides_tsp(distances)
    
    return tour, total_cost

if __name__ == "__main__":
    print(main())

