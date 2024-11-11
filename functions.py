import networkx as nx

def swapkeyval(indict):
    """
    Takes in a dictionary and returns a swapped key-value dictionary
    """
    outdict = dict()
    for k, v in indict.items():
        outdict[str(v)] = k
    return outdict


def density(cluster):
    """
    Calculate density of the entire graph
    For a cluster C with n nodes and m edges, density is defined as:
    Density(C) = 2m/n(n-1)
    """
    num_edges = 0
    num_nodes = len(cluster)
    
    # Count all edges
    for node, neighbors in cluster.items():
        num_edges += len(neighbors)
    
    num_edges //= 2
    
    if num_nodes < 2:
        return 0
    
    density = (2 * num_edges) / (num_nodes * (num_nodes - 1))
    return density

def check_connection(cluster):
    """
    Indicates if the given cluster (networkx class / 
    interaction_network class) is fully connected, 
    or if there are any other isolated vertices that aren't connected.
    """
    clusters = list()
    if nx.is_connected(cluster):
        return True, clusters
    else:
        clusters = list(nx.connected_components(cluster))
        return False, clusters
    

def shortest_path(vertex, all_vertices):
    # Setting up a dictionary of shortest paths and the start vertex is put in a list
    shortest_paths = dict()
    to_process = [str(vertex)]

    # While the "to_process" list is not empty, we do branch and bound
    while to_process:
        current_branch = to_process.pop()

        # We look through all neighbors. If it goes to a vertex already in the path, it's unoptimal and is discarded
        for neighbor in all_vertices[int(current_branch.split("_")[-1])]:
            if str(neighbor) in current_branch.split("_"):
                continue  # Skip if neighbor is already in the current branch

            # If the path doesn't loop, we add the neighbor to the current branch and check if it's a new shortest path
            new_path = current_branch + "_" + str(neighbor)
            to_process.append(new_path)
            start_end = current_branch.split("_")[0] + "->" + str(neighbor)
            
            if start_end in shortest_paths:
                # If there are more paths of equal length, we're working with a list
                if isinstance(shortest_paths[start_end], list):
                    if len(new_path.split("_")) < len(shortest_paths[start_end][0].split("_")):
                        shortest_paths[start_end] = new_path
                    elif len(new_path.split("_")) == len(shortest_paths[start_end][0].split("_")):
                        shortest_paths[start_end].append(new_path)
                else:
                    # If there's only one shortest path found until now
                    if len(new_path.split("_")) < len(shortest_paths[start_end].split("_")):
                        shortest_paths[start_end] = new_path
                    elif len(new_path.split("_")) == len(shortest_paths[start_end].split("_")):
                        shortest_paths[start_end] = [shortest_paths[start_end], new_path]
            else:
                # If no shortest path has been identified between the two points
                shortest_paths[start_end] = new_path
        
    # Returning the output
    return shortest_paths