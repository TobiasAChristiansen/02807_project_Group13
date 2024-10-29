import functions as f
import constants as c



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


    def load_data(self, data_path : str = c.datadir, filename: str = "short_alzheimers.tsv") -> None:
        infile = open(data_path + filename)
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
        