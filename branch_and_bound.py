import numpy as np
from queue import PriorityQueue

# Define infinity
INF = float('inf')

def branch_and_bound_tsp(distance_matrix):
    
    num_cities = len(distance_matrix)
    
    # Priority queue to manage nodes (sorted by lower bound)
    priority_queue = PriorityQueue()
    
    # Create a deep copy of the distance matrix to avoid modifying the original
    initial_matrix = np.copy(distance_matrix)
    
    # Initial reduction and calculation of initial lower bound
    reduced_matrix, initial_lower_bound = reduce_matrix(initial_matrix)
    
    # Insert root node into the queue
    # (lower bound, current cost, path, remaining cities, reduced matrix)
    start_city = 0  # Start from city A (index 0)
    priority_queue.put((initial_lower_bound, 0, [start_city], list(range(1, num_cities)), reduced_matrix))
    
    best_tour = None
    best_cost = float('inf')
    
    # Continue as long as there are nodes in the queue
    while not priority_queue.empty():
        # Extract the node with the lowest lower bound
        lower_bound, current_cost, path, remaining_cities, current_matrix = priority_queue.get()
        
        # If the lower bound is higher than the best cost so far, no need to continue with this branch
        if lower_bound >= best_cost:
            continue
        
        # If there are no more cities to add to the path
        if not remaining_cities:
            # Check if we can return to the starting city
            last_city = path[-1]
            if current_matrix[last_city][start_city] != INF:
                total_cost = current_cost + distance_matrix[last_city][start_city]
                if total_cost < best_cost:
                    # Add the starting city to complete the circuit
                    best_cost = total_cost
                    best_tour = path + [start_city]
            continue
        
        # Examine each possibility of adding a city to the path
        current_city = path[-1]
        for next_city in remaining_cities:
            # The distance from the current city to the next city
            edge_cost = current_matrix[current_city][next_city]
            
            # If there is no way to reach the next city, skip it
            if edge_cost == INF:
                continue
            
            # Copy the matrix to avoid modifying the original
            new_matrix = np.copy(current_matrix)
            
            # Mark the row of the current city as infinity
            new_matrix[current_city, :] = INF
            
            # Mark the column of the next city as infinity
            new_matrix[:, next_city] = INF
            
            # Prevent early return to the starting city (unless it's the last move)
            if len(remaining_cities) > 1:
                new_matrix[next_city][start_city] = INF
            
            # Perform reduction and calculate additional cost
            reduced_matrix, reduction_cost = reduce_matrix(new_matrix)
            
            # Calculate current cost and new lower bound
            new_cost = current_cost + distance_matrix[current_city][next_city]
            new_lower_bound = new_cost + reduction_cost
            
            # If the lower bound is higher than the best cost so far, no need to add to the queue
            if new_lower_bound >= best_cost:
                continue
            
            # Prepare the new path and remaining cities
            new_path = path + [next_city]
            new_remaining = [city for city in remaining_cities if city != next_city]
            
            # Add the new node to the queue
            priority_queue.put((new_lower_bound, new_cost, new_path, new_remaining, reduced_matrix))
    
    return best_tour, best_cost

def reduce_matrix(matrix):
    
    num_cities = len(matrix)
    reduction_cost = 0
    reduced_matrix = np.copy(matrix)
    
    # Row reduction - subtract the minimum value from each row
    for i in range(num_cities):
        row_min = min(reduced_matrix[i])
        if row_min != INF and row_min > 0:
            reduced_matrix[i] -= row_min
            reduction_cost += row_min
    
    # Column reduction - subtract the minimum value from each column
    for j in range(num_cities):
        col_min = min(reduced_matrix[:, j])
        if col_min != INF and col_min > 0:
            reduced_matrix[:, j] -= col_min
            reduction_cost += col_min
    
    return reduced_matrix, reduction_cost

# Usage example
if __name__ == "__main__":
    
    distance_matrix = np.array([
        [INF, 19, 21, 15, 12],
        [19, INF, 5, 11, 7],
        [21, 5, INF, 11, 7],
        [15, 11, 10, INF, 31],
        [12, 7, 20, 31, INF]
    ])
    
    # Run the algorithm
    optimal_tour, optimal_cost = branch_and_bound_tsp(distance_matrix)
    print(f"optimal_tour = {optimal_tour}\n" + f"optimal_cost = {optimal_cost}")
