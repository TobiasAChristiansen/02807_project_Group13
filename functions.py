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