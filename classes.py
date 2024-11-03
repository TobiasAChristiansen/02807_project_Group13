import functions as f
import constants as c
import pandas as pd



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


    def __str__(self):
        return "\n".join(f"{key} {value}" for key, value in self.vertices.items())


    def __repr__(self):
        return "\n".join(f"{key} {value}" for key, value in self.vertices.items())
    

    def create_encoding_dict(self, file_url: str = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz", compression : str = "gzip", sep : str = "\t"):
        
        #Loading the relevant data
        self.encoding_dict = pd.read_csv(file_url, compression=compression, sep=sep)

        #Isolating the string id
        self.encoding_dict = self.encoding_dict[["#string_protein_id"]]

        #Converting to a dictionary
        self.encoding_dict = self.encoding_dict.to_dict()

        #Swapping keys and values
        self.encoding_dict = f.swapkeyval(self.encoding_dict)


    def load_data(self, df : pd.DataFrame, file_url: str = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz", compression : str = "gzip", sep : str = "\t", req_experimental : bool = True):
        
        #Loading the data
        self.data = pd.read_csv(file_url, compression=compression, sep=sep)

        #If we only want interactions with experimental evidence:
        if req_experimental:
            self.data = self.data[self.data["experimental"] > 0]
        
        #Encoding the data frame:
        self.data = self.data.map(lambda x: self.encoding_dict[x] if x in self.encoding_dict else x)

        for line in infile:
            try:
                linedata = line.split()[:2]
                if linedata[0] in self.vertices:
                    self.vertices[linedata[0]].add(linedata[1])
                else:
                    self.vertices[linedata[0]] = {linedata[1]}
                if linedata[1] in self.vertices:
                    self.vertices[linedata[1]].add(linedata[0])
                else:
                    self.vertices[linedata[1]] = {linedata[0]}
            except:
                print(f"Line '{line}' was discarded")
        