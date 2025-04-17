import numpy as np

def nearest_neighbour(distance_matrix):
    
    n = distance_matrix.shape[0]  # Number of vertices
    INF = float('inf')
    
    # Start from vertex 0
    start_vertex = 0
    current_vertex = start_vertex
    unvisited = set(range(n)) - {start_vertex}
    tour = [start_vertex]
    total_distance = 0
    
    # Build the tour
    while unvisited:
        min_distance = INF
        nearest_vertex = None
        
        for vertex in unvisited:
            if distance_matrix[current_vertex, vertex] < min_distance:
                min_distance = distance_matrix[current_vertex, vertex]
                nearest_vertex = vertex
        
        tour.append(nearest_vertex)
        total_distance += min_distance
        current_vertex = nearest_vertex
        unvisited.remove(nearest_vertex)
    
    # Return to start vertex to complete the tour
    tour.append(start_vertex)
    total_distance += distance_matrix[current_vertex, start_vertex]
    
    return tour, total_distance

# Usage example with the matrix from your example
if __name__ == "__main__":
    INF = float('inf')  # Define infinity
    
    # The matrix specified in the example
    distance_matrix = np.array([
        [INF, 10, 15, 20],
        [10, INF, 35, 25],
        [15, 35, INF, 30],
        [20, 25, 30, INF]
    ])
    
    # Run the Nearest Neighbor algorithm
    tour, cost = nearest_neighbour(distance_matrix)
    print(f"Nearest Neighbor tour = {tour}")
    print(f"Nearest Neighbor cost = {cost}")
    
