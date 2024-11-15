import networkx as nx
import time

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

def check_connection(vertices):
    """
    Checks if the given vertices are all connected with every one of them.
    Give two outputs:
    - First output: returns True or False for the question "is it fully connected?"
    - Second output: returns the groups with vertices that are connected between each other.
    """
    # Get a list of all vertices in the graph
    all_vertices = list(vertices.keys())
    
    # To keep track of visited vertices
    visited = set()
    
    # List to store the connected components (groups)
    components = []
    
    def bfs(start_vertex):
        # Perform BFS starting from 'start_vertex'
        visited.add(start_vertex)
        queue = [start_vertex]
        component = [start_vertex]
        
        while queue:
            vertex = queue.pop(0)
            for neighbor in vertices[vertex]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
                    component.append(neighbor)
        return component
    
    # Iterate through all vertices and find all components
    for vertex in all_vertices:
        if vertex not in visited:
            component = bfs(vertex)
            components.append(component)
    
    return len(components) == 1, components    
    

def shortest_path(vertex, all_vertices, debug_mode=False, timed = False):
    # Setting up a dictionary of shortest paths and the start vertex is put in a list
    shortest_paths = dict()
    to_process = [str(vertex)]
    count = 0

    # While the "to_process" list is not empty, we do branch and bound
    while to_process:
        if debug_mode:
            print(count, len(to_process), len(shortest_paths))
        
        #Taking one of the short branches to process
        current_branch = to_process.pop(0)

        # We look through all neighbors. If it goes to a vertex already in the path, it's unoptimal and is discarded
        for neighbor in all_vertices[int(current_branch.split("_")[-1])]:
            
            # If the path doesn't loop, we add the neighbor to the current branch and check if it's a new shortest path
            new_path = current_branch + "_" + str(neighbor)
            start_end = new_path.split("_")[0] + "->" + str(neighbor)
            
            if start_end in shortest_paths:
                # If there are more paths of equal length, we're working with a list
                if isinstance(shortest_paths[start_end], list):
                    if len(new_path.split("_")) < len(shortest_paths[start_end][0].split("_")):
                        shortest_paths[start_end] = new_path
                        to_process.append(new_path)
                    elif len(new_path.split("_")) == len(shortest_paths[start_end][0].split("_")):
                        shortest_paths[start_end].append(new_path)
                        to_process.append(new_path)
                else:
                    # If there's only one shortest path found until now
                    if len(new_path.split("_")) < len(shortest_paths[start_end].split("_")):
                        shortest_paths[start_end] = new_path
                        to_process.append(new_path)
                    elif len(new_path.split("_")) == len(shortest_paths[start_end].split("_")):
                        shortest_paths[start_end] = [shortest_paths[start_end], new_path]
                        to_process.append(new_path)
            else:
                # If no shortest path has been identified between the two points
                shortest_paths[start_end] = new_path
                to_process.append(new_path)
        
    return shortest_paths