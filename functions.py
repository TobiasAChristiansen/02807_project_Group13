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



def count_occurances(counts_dict, edge_dict):
    for key in edge_dict:
        counts_dict[key] = counts_dict.get(key, 0) + edge_dict[key]
    return counts_dict


def lowest_first_from_to(u, v):
    return str(min([u,v])) + "->" + str(max([u,v]))

def lowest_first_from_to_edge(u, v):
    return str(min([u,v])) + "-" + str(max([u,v]))



def girvan_newman_modularity(graph, clusters):
    """
    Calculate modularity (Q) for a given interaction network and clusters, using Louvain Algorithm
    
    Parameters:
    - graph: a dictionary. Keys are protein IDs, and values are dictionaries of neighbors and their weights.
    - clusters: List of dictionaries, where each dictionary represents a cluster.
                Keys are protein IDs, and values are dictionaries of neighbors and their weights.
    
    Returns:
    - modularity (Q) value with range [-0.5, 1]
    """
    total_weight = sum(
        sum(neighbors.values()) for neighbors in graph.values()
    ) / 2  

    node_strength = {node: sum(neighbors.values()) for node, neighbors in graph.items()}

    modularity = 0.0

    for cluster in clusters:
        nodes_in_cluster = set(cluster.keys())

        in_cluster_weight = 0.0 
        total_strength = sum(node_strength[node] for node in nodes_in_cluster)  

        for node_i in nodes_in_cluster:
            neighbors = graph.get(node_i, {})
            for node_j, weight in neighbors.items():
                if node_j in nodes_in_cluster:
                    in_cluster_weight += weight

        in_cluster_weight /= 2.0

        modularity += (in_cluster_weight / total_weight) - (total_strength / (2 * total_weight)) ** 2

    return modularity


