import functions as f
import constants as c
import pandas as pd
import numpy as np
import networkx as nx 
import os, sys
from tqdm.notebook import tqdm
from functools import reduce
import multiprocessing



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
        self.graph_name = "PPI_GraphNetwork"
        self.graph_network = None
        self.shortest_paths = dict()
        self.testdataset = testdataset
        self.threshold = threshold
        self.NewPython3912 = NewPython3912


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
        self.data["combined_score"] = (self.data["combined_score"]-self.data["combined_score"].min())/(self.data["combined_score"].max()-self.data["combined_score"].min())
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
                    
                    # Update the sub-dictionary with the new key-value pair
                    self.vertices[int(linedata[0])][int(linedata[1])] = linedata[2]

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
                    
                    # Update the sub-dictionary with the new key-value pair
                    self.vertices[int(linedata[0])][int(linedata[1])] = linedata[2]

                except:
                    raise ValueError(f"Row {row} was discarded")
        
        #Deleting the data. We don't need it anymore, since we have loaded
        self.data = None
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
    

    def shortest_path(self, vertex, debug_mode=False, method="Dijkstra"):
        # Setting up a dictionary of shortest paths and the start vertex is put in a list
        shortest_paths = dict()
        to_process = [str(vertex)]
        count = 0

        if method.lower() == "bfs":

            # While the "to_process" list is not empty, we do branch and bound
            while to_process:
                if debug_mode:
                    print(count, len(to_process), len(shortest_paths))

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
                    print(count, len(to_process), len(self.shortest_paths))

                #Taking one of the short branches to process
                current_branch = to_process.pop(0)

                # We look through all neighbors. If it goes to a vertex already in the path, it's unoptimal and is discarded
                for neighbor in self.vertices[int(current_branch[0].split("_")[-1])]:

                    # If the path doesn't loop, we add the neighbor to the current branch and check if it's a new shortest path
                    #Adding the neigbor
                    new_path = [current_branch[0] + "_" + str(neighbor), current_branch[1]]
                    start_end = f.lowest_first_from_to(new_path[0].split("_")[0], str(neighbor))
                    new_path_split = new_path[0].split("_")

                    if start_end in self.shortest_paths:
                        continue
                    
                    #Updating the length of the path
                    new_path[1] = new_path[1] + (1 - self.vertices[int(new_path_split[-2])][int(new_path_split[-1])])

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
        if debug_mode: print("Initializing MapReduce of shortest_path()...") #debug_mode
        with multiprocessing.Pool() as pool:
            results = pool.map(self.shortest_path, cluster)
        occurances = reduce(f.count_occurances, results, {})

        #tweaking data structure with encoding dict
        if debug_mode: print("Initializing MapReduce of shortest_path()...") #debug_mode
        occurances_small = {}
        count=0
        if self.NewPython3912: #New python - no progress bar
            print("Creating encoding dict without progress bar...") #debug_mode
            for key, value in occurances.items():
                self.encoding_edgesused[count] = key
                occurances_small[count] = value
                count += 1
        else: #Old python - with progress bar
            for key, value in tqdm(occurances.items(), desc="Parsing encoding_edgesused data"):
                self.encoding_edgesused[count] = key
                occurances_small[count] = value
                count += 1

        return occurances_small





