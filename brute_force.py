from itertools import permutations

def brute_force_tsp(start_node, adjacency_matrix):
    
    num_nodes = len(adjacency_matrix)

    # List of all nodes except the start node
    other_nodes = [node for node in range(num_nodes) if node != start_node]

    # Initialize shortest path and weight
    shortest_path = None
    shortest_weight = float('inf')

    # Generate all possible orderings of the remaining nodes
    for perm in permutations(other_nodes):
        current_weight = 0
        current_node = start_node
        path = [start_node]  # Start path from the given start node

        # Calculate the weight of the current path
        for next_node in perm:
            current_weight += adjacency_matrix[current_node][next_node]
            current_node = next_node
            path.append(next_node)

        # Complete the cycle by returning to the start node
        current_weight += adjacency_matrix[current_node][start_node]
        path.append(start_node)

        # Check if this is the shortest path found so far
        if current_weight < shortest_weight:
            shortest_weight = current_weight
            shortest_path = path

    return shortest_path, shortest_weight


if __name__ == "__main__":
    adj_matrix = [
        [0, 19, 21, 15, 12],
        [19, 0, 5, 11, 7],
        [21, 5, 0, 11, 7],
        [15, 11, 10, 0, 31],
        [12, 7, 20, 31, 0]
    ]

    path, weight = brute_force_tsp(0, adj_matrix)
    print(f"Shortest path: {path}")
    print(f"Total weight: {weight}")
