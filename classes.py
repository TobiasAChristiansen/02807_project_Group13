import functions as f
import constants as c
import pandas as pd
import numpy as np
import networkx as nx 
import os, sys, joblib
from tqdm.notebook import tqdm
from functools import reduce
import multiprocessing
from IPython.display import clear_output



class interaction_network:
    """**Interaction network:**\n
    This class object is used to keep track of both vertices and edges in a connected network.
    Furthermore, it contains functions for creating clusters from raw data
    The class can be initialized with arguments for concurrent changes to several functions\n\n
    **Class functions:**\n
    load_data(data_path, filename): Loads network data from a tsv file\n\n
    **Class attributes:**\n
    self.vertices: a dictonary where the keys are the protein names, and the value is a list of protein names of interacting proteins\n
    self.testdataset: applies a smaller test data from string (not human proteins)\n
    self.threshold: given threshold in percent TO REMOVE - meaning threshold=0.05 removes lowest 5% values of combined_score. Must be between 0 and 1
    self.NewPython3912: NewPython3912 is given as true, if user is using python version 3.9.12 or older (in that case, tqdm() parsing of data wont run and progress bar will not appear!)
    """

    def __init__(self, 
                 testdataset : bool = False,
                 threshold : int = 0,
                 NewPython3912 : bool = False):
        self.vertices = dict()
        self.data = None
        self.encoding_dict = None
        self.encoding_edgesused = {}
        self.occurances_small = None
        self.graph_name = "PPI_GraphNetwork"
        self.graph_network = None
        self.shortest_paths = dict()
        self.testdataset = testdataset
        self.threshold = threshold
        self.NewPython3912 = NewPython3912
        self.finished_clusters = list()
        self.CollectMostTravelled = list() #Used for node size in graphing
        self.CurrentSeveranceScore = None #Used for node size in graphing


    def __str__(self):
        return "\n".join(f"{key} {value}" for key, value in self.vertices.items())




    def __repr__(self):
        return "\n".join(f"{key} {value}" for key, value in self.vertices.items())
    

    def create_encoding_dict(self, 
                             file_url: str = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz", 
                             compression : str = "gzip", 
                             sep : str = "\t"):
        if self.testdataset:
            file_url = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"
        
        #Message
        print("Creating encoding table")

        #Loading the relevant data
        self.encoding_dict = pd.read_csv(file_url, compression=compression, sep=sep)

        #Isolating the string id
        self.encoding_dict = self.encoding_dict[["#string_protein_id"]].values

        #Converting to a dictionary
        self.encoding_dict = {str(c):e[0] for c, e in enumerate(self.encoding_dict)}

        #Swapping keys and values
        self.encoding_dict = f.swapkeyval(self.encoding_dict)


    def load_data(self,
                  file_url: str = "https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz", 
                  compression : str = "gzip", 
                  sep : str = " ", 
                  req_experimental : bool = True):
        #input control threshold
        if self.threshold > 100 or self.threshold < 0:
            raise ValueError("Give threshold in percent, between 0 and 100%")

        #Test dataset
        if self.testdataset:
            file_url = "https://stringdb-downloads.org/download/stream/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.onlyAB.txt.gz"

        #Checking if an encoding dict have been created:
        if not self.encoding_dict:
            raise SyntaxError("self.create_encoding_dict() should be run prior to load_data")

        #Loading the data
        print("Fetching data")
        self.data = pd.read_csv(file_url, compression=compression, sep=sep)

        #If we only want interactions with experimental evidence:
        if req_experimental:
            self.data = self.data[self.data["experimental"] > 0]
        

        #Encoding the data frame:
        print("Cropping data")
        self.data = self.data.apply(lambda series: series.map( lambda x: self.encoding_dict[x] if x in self.encoding_dict else x)).reset_index()
        self.data = self.data[["protein1", "protein2", "combined_score"]]
        threshold_value = self.data["combined_score"].quantile(self.threshold) #calculate combined_score value associated with given threshold percentage
        self.data = self.data[self.data["combined_score"] > threshold_value] #filter threshold
        
        #Return info on threshold value
        if threshold_value > 0: 
            print(f"Given threshold: {self.threshold} filtered out combined_score values: {threshold_value} and below")
        else:
            print(f"Given threshold: {self.threshold} filtered out no values")
        
        #Min-max normalizing of "combined_score"
        self.data["combined_score"] = (self.data["combined_score"]-self.data["combined_score"].min())/(self.data["combined_score"].max()-self.data["combined_score"].min()) if self.data["combined_score"].min() != self.data["combined_score"].max() else self.data["combined_score"] / self.data["combined_score"].max()
        self.data["combined_score"] = self.data["combined_score"].round(2)


        #Flag for NewPython3912 (tqdm() wont run in newer versions)
        if self.NewPython3912:
            print("Parsing data without progress bar...")
            for i in range(len(self.data)):
                row = self.data.iloc[i]
                try:
                    linedata = [row["protein1"], row["protein2"], row["combined_score"]]
                    # Check if the main key exists; if not, initialize it as an empty dictionary
                    if int(linedata[0]) not in self.vertices:
                        self.vertices[int(linedata[0])] = {}
                    if int(linedata[1]) not in self.vertices:
                        self.vertices[int(linedata[1])] = {}
                    
                    # Update the sub-dictionary with the new key-value pair
                    self.vertices[int(linedata[0])][int(linedata[1])] = linedata[2]
                    self.vertices[int(linedata[1])][int(linedata[0])] = linedata[2]

                except:
                    raise ValueError(f"Row {row} was discarded")
        
        else: #in case python version is "old" and tqdm() can run   
            for i in tqdm(range(len(self.data)), desc="Parsing data"):
                row = self.data.iloc[i]
                try:
                    linedata = [row["protein1"], row["protein2"], row["combined_score"]]
                    # Check if the main key exists; if not, initialize it as an empty dictionary
                    if int(linedata[0]) not in self.vertices:
                        self.vertices[int(linedata[0])] = {}
                    if int(linedata[1]) not in self.vertices:
                        self.vertices[int(linedata[1])] = {}
                    
                    # Update the sub-dictionary with the new key-value pair
                    self.vertices[int(linedata[0])][int(linedata[1])] = linedata[2]
                    self.vertices[int(linedata[1])][int(linedata[0])] = linedata[2]

                except:
                    raise ValueError(f"Row {row} was discarded")
        
        #Deleting the data. We don't need it anymore, since we have loaded
        self.data = None

        #Putting the self.vertices item into a list. These will be the preliminary clusters
        self.vertices = [self.vertices]

        #Check connectivity of the cluster
        is_connected, connected_keys = f.check_connection(self.vertices[0])

        #If not all vertices are interconnected, we split the cluster
        if not is_connected:
            self.split_cluster(0, connected_keys)

        #Printing the size of the data structure
        print(f"Data loaded in var .vertices, with size: {sys.getsizeof(self.vertices)/1000000}mb")


    def construct_graph(self,
                        interaction_data: dict):
        """
        Initializes graph network and constructs edges, based on input data.
        Input data should be in the format: dict(dict()), where the inner and outer key is a protein id, and the inner value is the normalized combined score.
        Example: interaction_data[0][9827] = 0.3; the interaction between protein 0 and 9827, has (normalized) probability = 0.3
        """
        self.graph_network = nx.Graph(name=self.graph_name) #initialize empty graph
        
        #Creating nodes in network
        allproteins = set()
        for root in interaction_data.keys():
            allproteins.add(root)
            for neighbor in interaction_data.keys():
                allproteins.add(neighbor)
        self.graph_network.add_nodes_from(allproteins)

        #Simple adding of edges to graph network - Asbjørn
        for root in interaction_data.keys():
            for neighbor in interaction_data[root].keys():
                if root != neighbor:
                    self.graph_network.add_edge(root, neighbor)
    

    def shortest_path(self, vertex, debug_mode=False, method="dijkstra"):
        # Setting up a dictionary of shortest paths and the start vertex is put in a list
        shortest_paths = dict()
        to_process = [str(vertex)]

        if method.lower() == "bfs":
            raise NotImplementedError

            # While the "to_process" list is not empty, we do branch and bound
            while to_process:
                if debug_mode:
                    print(len(to_process), len(shortest_paths))

                #Taking one of the short branches to process
                current_branch = to_process.pop(0)

                # We look through all neighbors. If it goes to a vertex already in the path, it's unoptimal and is discarded
                for neighbor in self.vertices[int(current_branch.split("_")[-1])]:

                    # If the path doesn't loop, we add the neighbor to the current branch and check if it's a new shortest path
                    new_path = current_branch + "_" + str(neighbor)
                    start_end = new_path.split("_")[0] + "->" + str(neighbor)

                    if start_end in shortest_paths:
                        # If there are more paths of equal length, we're working with a list
                        if isinstance(shortest_paths[start_end], list):
                            if len(new_path.split("_")) < len(shortest_paths[start_end][0].split("_")):
                                shortest_paths[start_end] = new_path
                                to_process.append(new_path)
                            #elif len(new_path.split("_")) == len(shortest_paths[start_end][0].split("_")):
                            #    shortest_paths[start_end].append(new_path)
                            #    to_process.append(new_path)
                        else:
                            # If there's only one shortest path found until now
                            if len(new_path.split("_")) < len(shortest_paths[start_end].split("_")):
                                shortest_paths[start_end] = new_path
                                to_process.append(new_path)
                            #elif len(new_path.split("_")) == len(shortest_paths[start_end].split("_")):
                            #    shortest_paths[start_end] = [shortest_paths[start_end], new_path]
                            #    to_process.append(new_path)
                    else:
                        # If no shortest path has been identified between the two points
                        shortest_paths[start_end] = new_path
                        to_process.append(new_path)

            return shortest_paths

        elif method.lower() == "dijkstra":
            #Adding an initial distance of 0
            to_process = [[str(vertex), 0]]

            # While the "to_process" list is not empty, we do branch and bound
            while to_process:
                if debug_mode:
                    print(len(to_process), len(self.shortest_paths))

                #Taking one of the short branches to process
                current_branch = to_process.pop(0)

                if debug_mode:
                    print(current_branch[0])

                # We look through all neighbors. If it goes to a vertex already in the path, it's unoptimal and is discarded
                for neighbor in self.vertices[0][int(current_branch[0].split("_")[-1])]:
                    if str(neighbor) in current_branch[0].split("_"):
                        continue

                    # If the path doesn't loop, we add the neighbor to the current branch and check if it's a new shortest path
                    #Adding the neigbor
                    new_path = [current_branch[0] + "_" + str(neighbor), current_branch[1]]
                    start_end = f.lowest_first_from_to(new_path[0].split("_")[0], str(neighbor))
                    new_path_split = new_path[0].split("_")
                    new_path[1] = new_path[1] + (1 - self.vertices[0][int(new_path_split[-2])][int(new_path_split[-1])])



                    if start_end in self.shortest_paths:
                        if new_path[1] < self.shortest_paths[start_end][1]:
                            self.shortest_paths[start_end] = new_path
                            to_process.append(new_path)
                    else:
                        # If no shortest path has been identified between the two points
                        self.shortest_paths[start_end] = new_path
                        to_process.append(new_path)

            #Returning a list of all edges A-B
            edges_dict = dict()
            for item in self.shortest_paths:
                for i in range(1, len(split_path := self.shortest_paths[item][0].split("_"))):
                    edge = f.lowest_first_from_to_edge(split_path[i-1], split_path[i])
                    edges_dict[edge] = edges_dict.get(edge, 0) + 1
            #num_included = 1000
            #keys_to_include = [item[0] for item in sorted(edges_dict.items(), key = lambda x: x[1], reverse=True)[:num_included]]
            return edges_dict #{key:edges_dict[key] for key in keys_to_include}


    def evaluate_most_used_path(self, cluster, debug_mode=False):
        """This function uses mapreduce to count the number of times an edge is in a shortest path. 
        The function returns a dictionary where they keys are edges and the value is counts"""
        if debug_mode:
            print("Initializing MapReduce of shortest_path()...") #debug_mode
            counter = 0
            len_cluster = len(cluster)
            for vertex in cluster:
                counter += 1
                print(f"eval_most_used_path_progress = {counter}/{len_cluster}")
                results = self.shortest_path(vertex)
        else:
            #Mapping the shortest_path function onto the cluster
            print("-----Mapping-----")
            results = f.custom_pool(self.shortest_path, cluster, ProgressBar=True)

            #Reducing the result by using count_occurances
            print("-----Reducing-----")
            occurances = reduce(f.count_occurances, results, {})

        #tweaking data structure with encoding dict. This was not needed
        #if debug_mode: print("Initializing MapReduce of shortest_path()...") #debug_mode
        #occurances_small = {}
        #count=0
        #if self.NewPython3912: #New python - no progress bar
        #    print("Creating encoding dict without progress bar...") #debug_mode
        #    for key, value in occurances.items():
        #        self.encoding_edgesused[count] = key
        #        occurances_small[count] = value
        #        count += 1
        #else: #Old python - with progress bar
        #    for key, value in tqdm(occurances.items(), desc="Parsing encoding_edgesused data"):
        #        self.encoding_edgesused[count] = key
        #        occurances_small[count] = value
        #        count += 1

        #self.occurances_small = occurances_small
        return occurances# occurances_small
    

    def severance_score(self, shortest_paths):
        """Calculating the severance score by dividing by the weight. The weak edges will be preferred"""
        for key in shortest_paths:
            start_end = key.split("-")
            shortest_paths[key] = shortest_paths[key] / self.vertices[0][int(start_end[0])][int(start_end[1])]
        return shortest_paths
    
    def find_edge_to_cut(self, severance_scores):
        edge_to_cut = sorted(severance_scores.items(), key = lambda x:x[1], reverse=True)[0][0]
        return edge_to_cut

    def cut_edge(self, clusternumber, edge_to_cut):
        start_end = edge_to_cut.split("-")
        path_format = "_".join(start_end)

        #Deleting the edge from the data
        del self.vertices[clusternumber][int(start_end[0])][int(start_end[1])]
        del self.vertices[clusternumber][int(start_end[1])][int(start_end[0])]

        #Deleting shortest paths that rely on this edge.
        for key in self.shortest_paths:
            if path_format in self.shortest_paths[key]:
                del self.shortest_paths[key]
    
    def split_cluster(self, clusterindex, list_of_keylists):
        """Splitting the cluster using a list of lists of keys that should be grouped together"""
        new_clusters = list()
        for keylist in list_of_keylists:
            new_clusters.append(dict())
            for key in keylist:
                new_clusters[-1][key] = self.vertices[clusterindex][key]
        del self.vertices[clusterindex]
        for new_cluster in new_clusters:
            self.vertices.append(new_cluster)

    def edge_to_remove(self):
        """Function for determining which edge should be removed from the cluster"""
        print("-----Finding edge to remove-----")
        res = self.evaluate_most_used_path(self.vertices[0])
        res = self.severance_score(res)
        res = self.find_edge_to_cut(res)
        return res

    def cluster(self):
        #While there are still clusters to be processed, we run a loop of cutting and evaluating
        edge_to_remove = ""
        while len(self.vertices) > 0:
            print(f"NEW ITERATION --- Clusters left: {len(self.vertices)} --- Cluster length: {len(self.vertices[0])} --- Finished clusters: {len(self.finished_clusters)} --- Latest cluster length: {len(self.finished_clusters[-1][1].values()) if self.finished_clusters else 0} --- last edge severed: {edge_to_remove}")
            
            #Determine density of the cluster
            print("-----Density calculation-----")
            density = f.density(self.vertices[0])
            #print(f"-----Density: {density}-----")

            #If the conditions are met, we classify by GSEA
            if density >= 0.7:
                print("-----GSEA-----")
                decoded_cluster = self.vertices[0] #self.decode_edgesused(cluster=self.vertices[0])
                self.finished_clusters.append([f.enrichment_analysis(list(decoded_cluster.keys())), decoded_cluster])
                del self.vertices[0]

            #If it does not satisfy the density condition, we find an edge to cut
            else:
                print("-----Evaluate edge removal-----")
                copy_current_cluster = self.vertices[0].copy()   #Making a copy of the current cluster so we can revert the cut
                edge_to_remove = self.edge_to_remove()   #Identifying which edge to remove
                self.cut_edge(0, edge_to_remove) #Removing the edge

                #Check connectivity of the cluster
                is_connected, connected_keys = f.check_connection(self.vertices[0])

                #If not all vertices are interconnected, we split the cluster
                if not is_connected:
                    self.split_cluster(0, connected_keys)

                    #modularity calculation:
                    print("-----Modularity-----")
                    modularity = f.girvan_newman_modularity(copy_current_cluster, self.vertices[-len(connected_keys):])

                    #If the modularity will not be increased with the split, we restore the edge and finalize the cluster
                    if modularity < 0:
                        print("-----Finalized previous cluster-----")
                        self.CollectMostTravelled.append(self.CurrentSeveranceScore) #Collect severance scores if cluster is completed - used in graphing
                        decoded_cluster = copy_current_cluster #self.decode_edgesused(cluster=copy_current_cluster)
                        self.finished_clusters.append([f.enrichment_analysis(list(decoded_cluster.keys())), decoded_cluster])
                        self.vertices = self.vertices[:-3]
                    
                    #If the modularity is increased, we finalize the split
                    else:
                        print("-----Finalized split-----")


            clear_output(wait = False) #Clear output after each iteration
                        

        
        self.write_clusters()

    def write_clusters(self):
        joblib.dump(self.finished_clusters, "./joblib_vars/finished_clusters.joblib.gz")
        #outfile = open("results.txt", "w")
        #for finished_cluster in self.finished_clusters:
        #    outfile.write(finished_cluster[0] + "\n")
        #    for key in finished_cluster[1]:
        #        outfile.write(key + ": " + ", ".join(finished_cluster[1][key]))
        #outfile.close()

    def decode_edgesused(self, id : int = 0, cluster : dict = None, full_dict : bool = False):
        """
        Decodes edges used integers to string interactions.
        eg. 0 --> 0-1940 (not real interaction).

        cluster = insert cluster ; returns a decoded cluster in case cluster is given        
        full_dict = True ; Return full dictionary, not just interaction for id given

        For decoding protein name, use "decode_proteinname()"
        """
        return_dict = {}
        if full_dict: #return full dictionary, not just interaction for id given
            for key, value in self.occurances_small.items():
                return_dict[self.encoding_edgesused[key]] = value
            return return_dict
        
        elif cluster != None:
            for key, value in cluster.items():
                return_dict[self.encoding_edgesused[key]] = value
            return return_dict

        else:
            return self.encoding_edgesused[key]

    def decode_proteinname(self, id : int = 0, cluster : dict = None, full_dict : bool = False):
        """
        Decodes protein used integers to string names as given in stringDB.
        eg. 0 --> ENV09918398 (not real protein).

        cluster = insert cluster ; returns a decoded cluster in case cluster is given 
        full_dict = True ;Return full dictionary, not just protein name for id given

        For decoding edges in cluster, use "decode_edgesused()"
        """
        return_dict = {}
        if full_dict: #return full dictionary, not just protein name for id given
            for key, value in self.occurances_small.items():
                return_dict[self.encoding_dict[key]] = value
            return return_dict

        elif cluster != None:
            for key, value in cluster.items():
                return_dict[self.encoding_edgesused[key]] = value
            return return_dict

        else:
            return self.encoding_dict[key]

    def construct_graph(self, input, graph_name : str = "PPI_GraphNetwork", reformat : bool = False, debug_mode : bool = False):
        """
        Initializes graph network and constructs edges, based on input data.
        Input data should be in the format: dict(dict()), where the inner and outer key is a protein id, and the inner value is the normalized combined score.
        Example: interaction_data[0][9827] = 0.3; the interaction between protein 0 and 9827, has (normalized) probability = 0.3
        *reformat [WIP]*  = True: converts input data to appropriate netwulf object for nx.visualize() - Running reformat_clustervar(). Else, iterates through input adding edges and nodes one at a time
        """
        self.graph_network = nx.Graph(name=self.graph_name) #initialize empty graph
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
                    self.graph_network.add_node(neighbor, size = self.NodeSize(node=neighbor)) #Node size = sum of normalized combined score
                    if debug_mode: print(f"--- Node collected: {neighbor} ---")
        #graph_network.add_nodes_from(allproteins)

        #Simple adding of edges to graph network - Asbjørn
        if debug_mode: print("##### Adding edges to graph #####")
        for cluster in data:
            for root in cluster.keys():
                for neighbor in cluster[root].keys():
                    if root != neighbor:
                        self.graph_network.add_edge(root, neighbor)
                        if debug_mode: print(f"--- Edge added between: {root, neighbor} ---")


    def NodeSize(self, node):
        """
        Finding node size from sum of normalized combined score, using self.data
        """
        sum_NormCombinedScore = 0
        for prot_col in ["protein1", "protein2"]:
            curr_data = self.data[self.data[prot_col] == node]
            sum_NormCombinedScore += curr_data["combined_score"].sum()
        return sum_NormCombinedScore


