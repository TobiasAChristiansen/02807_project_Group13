import networkx as nx
import time
from gprofiler.gprofiler import GProfiler
import numpy as np
import multiprocessing
from tqdm.notebook import tqdm


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

        
        try:
            modularity += (in_cluster_weight / total_weight) - (total_strength / (2 * total_weight)) ** 2
        except:
            print(total_weight, graph)
            raise ValueError()


    return modularity

def enrichment_analysis(protein_list, organism = "hsapiens", sign_level = 0.05):
    
    # Initialize GProfiler object
    print(protein_list)
    gp = GProfiler(return_dataframe=True)

    #Generate GProfiler with the information from protein sequences and organism
    results = gp.profile(
        organism=organism,
        query=protein_list,
        sources=['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC'],  # Include GO (BP, MF, CC), KEGG, and Reactome
        significance_threshold_method='fdr'  # FDR for multiple testing correction
    )
    
    if len(results) == 0:
        return "Missing more than one found protein"

    #Extract p-value from most significant term
    pval = results["p_value"].iloc[0]

    #Return most significant term + p-value
    if pval < sign_level:
        return [results["name"].iloc[0], results["p_value"].iloc[0]]

    else:
        return "No significant functions found"

def decode_edgesused(id : int):
    """
    Decodes edges used integers to string interactions.
    eg. 0 --> 0-1940 (not real interaction).

    For decoding protein name, use "decode_proteinname()"
    """

def decode_proteinname():
    """
    Decodes protein used integers to string names as given in stringDB.
    eg. 0 --> ENV09918398 (not real protein).

    For decoding edges in cluster, use "decode_edgesused()"
    """

def construct_graph(input, graph_name : str = "PPI_GraphNetwork", reformat : bool = False, debug_mode : bool = False):
        """
        Initializes graph network and constructs edges, based on input data.
        Input data should be in the format: dict(dict()), where the inner and outer key is a protein id, and the inner value is the normalized combined score.
        Example: interaction_data[0][9827] = 0.3; the interaction between protein 0 and 9827, has (normalized) probability = 0.3
        *reformat [WIP]*  = True: converts input data to appropriate netwulf object for nx.visualize() - Running reformat_clustervar(). Else, iterates through input adding edges and nodes one at a time
        """
        graph_network = nx.Graph(name=graph_name) #initialize empty graph
        if debug_mode: print("##### Splitting label and cluster dict #####")
        enrichment_labels = []
        data = []
        for entry in input:
            if debug_mode: print(entry)
            enrichment_labels.append(entry[0])
            data.append(entry[1])

        #Creating nodes in network
        if debug_mode: print("##### Adding nodes to graph #####")
        allproteins = set()
        for cluster in data:
            for root in cluster.keys():
                allproteins.add(root)
                for neighbor in cluster.keys():
                    allproteins.add(neighbor)
                    graph_network.add_node(neighbor, size = np.random.random())
                    if debug_mode: print(f"--- Node collected: {neighbor} ---")
        #graph_network.add_nodes_from(allproteins)

        #Simple adding of edges to graph network - Asbjørn
        if debug_mode: print("##### Adding edges to graph #####")
        for cluster in data:
            for root in cluster.keys():
                for neighbor in cluster[root].keys():
                    if root != neighbor:
                        graph_network.add_edge(root, neighbor)
                        if debug_mode: print(f"--- Edge added between: {root, neighbor} ---")

        return graph_network

def reformat_clustervar(data):
    """
    NB: not developed yet!
    Converts input data to appropriate netwulf object for nx.visualize()
    """

def custom_pool(job, target, ProgressBar : bool = True):
    """
    Keeps track of child processes in map reduce when multiprocessing.pool is used
    *job*: custom task to be performed
    *target*: target to perform task on
    """
    #Logfile for tracking - in dir <CustomPool_logs>
    #var init
    collected_results = []
    total_iters = len(target)
    
    #Run with tqdm()
    if ProgressBar: 
        with multiprocessing.Pool() as pool:
            for result in tqdm(pool.imap_unordered(job, target), total=len(target), desc="Child processes completed:"):
                collected_results.append(result)
    
    #Run without tqdm - Progress saved to CustomPool_logs
    else: 
        with multiprocessing.Pool() as pool:
            iter = pool.imap_unordered(job, target)
            
            while True:
                with open(f"./CustomPool_logs/CompletedContent.txt", "w") as logfile: #each ChildProcess writes to own logfile
                    try:
                        result = next(iter)
                    
                    except StopIteration:                            
                    # All jobs have been processed
                        print("All child jobs completed.")
                        break
                        
                    except ChildProcessError as err:
                        with open("./CustomPool_logs/FailedPools.txt", "w") as logfile_fail:
                            print(f"Processing of {err.args[0]} job failed.", file=logfile_fail)
                        
                    else:
                        collected_results.append(result)
                        with open(f"./CustomPool_logs/Pool_Iteration_Completed.txt", "w") as iterfile:
                            print(f"Job {len(collected_results)} completed!", file=iterfile)
                        print(f"Job {len(collected_results)} completed! \n{result}", file=logfile)

    if collected_results:
        with open("./CustomPool_logs/AllPool_logs.txt", "w") as logfile:
            print('All completed jobs:', collected_results, file=logfile)
    
    return collected_results


