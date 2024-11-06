import functions as f
import constants as c
import pandas as pd
import numpy as np
from tqdm.notebook import tqdm




class interaction_network:
    """**Interaction network:**\n
    This class object is used to keep track of both vertices and edges in a connected network.
    Furthermore, it contains functions for creating clusters from raw data\n\n
    **Class functions:**\n
    load_data(data_path, filename): Loads network data from a tsv file\n\n
    **Class attributes:**\n
    self.vertices: a dictonary where the keys are the protein names, and the value is a list of protein names of interacting proteins"""

    def __init__(self):
        self.vertices = dict()
        self.data = None
        self.encoding_dict = None



    def __str__(self):
        return "\n".join(f"{key} {value}" for key, value in self.vertices.items())




    def __repr__(self):
        return "\n".join(f"{key} {value}" for key, value in self.vertices.items())
    



    def create_encoding_dict(self, 
                             file_url: str = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz", 
                             compression : str = "gzip", 
                             sep : str = "\t"):
        
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
        self.data["combined_score"] = (self.data["combined_score"]-self.data["combined_score"].min())/(self.data["combined_score"].max()-self.data["combined_score"].min()) #Min-max normalizing of "combined_score"

        for i in tqdm(range(len(self.data)), desc="Parsing data"):
            row = self.data.iloc[i]
            try:
                linedata = [row["protein1"], row["protein2"], row["combined_score"]]
                if linedata[0] in self.vertices:
                    self.vertices[linedata[0]].add(linedata[1])
                else:
                    self.vertices[linedata[0]] = {linedata[1]}
                if linedata[1] in self.vertices:
                    self.vertices[linedata[1]].add(linedata[0])
                else:
                    self.vertices[linedata[1]] = {linedata[0]}
            except:
                raise ValueError(f"Row {row} was discarded")
        
        #Deleting the data. We don't need it anymore, since we have loaded
        self.data = None