# NetworKit Python Shell


from _NetworKit import *    # import the extension module

# standard library modules
import math
import textwrap
import collections

# non standard library modules
import networkx as nx
import tabulate

# local modules
import termgraph
import stopwatch

#--------- NetworKit Python Shell functions ----------------#

def readGraph(path, format=None):
    """    Read graph and return a NetworKit::Graph"""
    # TODO: detect file format by looking at the file content
    if format is None:    # look at file extension
        if path.endswith(".graph"):
            reader = METISGraphReader()
        else:
            raise Exception("unknown graph file format")
    else:
        if (format is "metis"):
            reader = METISGraphReader()
        elif (format is "edgelist"):
            reader = EdgeListIO('\t', 1)

    with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
        G = reader.read(path)
        return G

    return None
    



########  CONVERSION ########

def nx2nk(nxG, weightAttr=None):
    """ 
    Convert a networkx.Graph to a NetworKit.Graph
        :param weightAttr: the edge attribute which should be treated as the edge weight
     """
    # TODO: consider weights
    n = nxG.number_of_nodes()
    nkG = Graph(n)
    
    if weightAttr is not None:
        nkG.markAsWeighted()
        for (u, v) in nxG.edges():
            w = nxG.edge[u][v][weightAttr]
            nkG.addEdge(u, v, w)
    else:
        for (u, v) in nxG.edges():
            nkG.addEdge(u, v)
    
    return nkG


def nk2nx(nkG):
    """ Convert a NetworKit.Graph to a networkx.Graph """
    nxG = nx.Graph()
    if nkG.isMarkedAsWeighted():
        for (u, v) in nkG.edges():
            nxG.add_edge(u, v, weight=nkG.weight(u, v))
    else:
        for (u, v) in nkG.edges():
            nxG.add_edge(u, v)
    return nxG
    

########  PROPERTIES ########


def nm(nkG):
    n = nkG.numberOfNodes()
    m = nkG.numberOfEdges()
    return (n, m)

def degrees(nkG):
    minMaxDeg = GraphProperties.minMaxDegree(nkG)
    avgDeg = GraphProperties.averageDegree(nkG)
    return (minMaxDeg[0], minMaxDeg[1], avgDeg)

def components(nxG):
    """Analyze connected components"""
    components = nx.connected_components(nxG)
    nComponents = len(components)
    componentSizes = [len(component) for component in components]
    componentSizes.sort(reverse=True) # sort in descending order
    sizeLargestComponent = componentSizes[0]
    return (nComponents, sizeLargestComponent)


def properties(nkG, settings):
    nxG = None
    if settings["networkx"]:
        print("[...] converting to NetworX.Graph for some properties....")
        nxG = nk2nx(nkG)

    print("[...] calculating basic properties")
    
    # size
    n, m = nm(nkG)    # n and m

    # degree
    minDeg, maxDeg, avgDeg = degrees(nkG)

        # number of isolated nodes
    isolates = len(nx.isolates(nxG)) if nxG else None

    # number of cliques
    cliques = len(list(nx.find_cliques(nxG))) if nxG else None


    # number of self-loops
    loops = len(nxG.selfloop_edges()) if nxG else None
    
    # density
    dens = nx.density(nxG) if nxG else None

    # diameter
    dia = None
    if settings["diameter"] and (n < 1000):
        print("calculating diameter...")
        dia = nx.diameter(nxG)


    # calculate eccentricity
    ecc = None
    if settings["eccentricity"] and (n < 1000):
        print("calculating eccentricity...")
        eccentricities = nx.eccentricity(nxG)
        ecc = sum(val for val in eccentricities.values()) / n


    # community detection

    ncomPLP, modPLP = None, None
    ncomPLM, modPLM = None, None
    if settings["communities"]:
        print("[...] detecting communities")
        # perform PLP and PLM community detection
        print("performing community detection: PLP")
        # TODO: avoid printout of status bar
        plp = LabelPropagation()
        zetaPLP = plp.run(nkG)
        ncomPLP = zetaPLP.numberOfClusters()
        modPLP = Modularity().getQuality(zetaPLP, nkG)
        print("performing community detection: PLM")
        PLM = Louvain("balanced")
        zetaPLM = PLM.run(nkG)
        ncomPLM = zetaPLM.numberOfClusters()
        modPLM = Modularity().getQuality(zetaPLM, nkG)

    # degree histogram
    
    labels, histo = None, None
    if settings["degreeDistribution"]:
        print("[...] calculating degree histogram")    
        histo = nx.degree_histogram(nxG)
        (labels, histo) = compressHistogram(histo, nbins=25)

    # connected components
    nComponents, sizeLargestComponent = None, None
    if settings["components"]:
        print("[...] finding connected components")    
        nComponents, sizeLargestComponent = components(nxG)

    # clustering
    avglcc = None
    if settings["clustering"]:
        avglcc = GraphProperties.averageLocalClusteringCoefficient(nkG)

    # degree assortativity
    assort = None
    if settings["assortativity"]:
        assort = nx.degree_assortativity_coefficient(nxG)




    # betweenness centrality
    # TODO: average betweenness centrality?

    props = {
         "name": nkG.getName(),
         "n": n,
         "m": m,
         "minDeg": minDeg,
         "maxDeg": maxDeg,
         "avgDeg": avgDeg,
         "avglcc": avglcc,
         "nComponents": nComponents,
         "sizeLargestComponent": sizeLargestComponent,
         "dia": dia,
         "ecc": ecc,
         "isolates": isolates,
         "loops": loops,
         "ncomPLP": ncomPLP,
         "modPLP": modPLP,
         "ncomPLM": ncomPLM,
         "modPLM": modPLM,
         "dens": dens,
         "assort": assort,
         "cliques": cliques,
         "histo": (labels, histo),
         }

    return props


def showProperties(nkG, settings=collections.defaultdict(lambda: True)):
    props = properties(nkG, settings)
    basicProperties = [
        ["nodes (n)", props["n"]],
        ["edges (m)", props["m"]],
        ["min. degree", props["minDeg"]],
        ["max. degree", props["maxDeg"]],
        ["avg. degree", props["avgDeg"]],
        ["isolated nodes", props["isolates"]],
        ["self-loops", props["loops"]],
        ["density", "{0:.6f}".format(props["dens"]) if props["dens"] else None]
    ]
    pathStructure = [
        ["connected components", props["nComponents"]],
        ["size of largest component", props["sizeLargestComponent"]],
        ["diameter", props["dia"]],
        ["avg. eccentricity", props["ecc"]],
    ]
    
    miscProperties = [
        ["degree assortativity", "{0:.6f}".format(props["assort"]) if props["assort"] else None],
        ["cliques", props["cliques"]]
    ]

    communityStructure = [
        ["avg. local clustering coefficient", "", "{0:.6f}".format(props["avglcc"]) if props["avglcc"] else None],
        ["PLP community detection", "", ""],
        ["", "communities", props["ncomPLP"]],
        ["", "modularity", "{0:.6f}".format(props["modPLP"]) if props["modPLP"] else None],
        ["PLM community detection", "", ""],
        ["", "communities", props["ncomPLM"]],
        ["", "modularity", "{0:.6f}".format(props["modPLM"]) if props["modPLM"] else None],
    ]

    print()
    print("Network Properties")
    print("==================")
    print("Basic Properties")
    print(tabulate.tabulate(basicProperties))
    print("Path Structure")
    print(tabulate.tabulate(pathStructure))
    print("Miscellaneous")
    print(tabulate.tabulate(miscProperties))
    print("Community Structure")
    print(tabulate.tabulate(communityStructure))
    print("Degree Distribution")
    print("-------------------")
    (labels, histo) = props["histo"]
    if labels and histo:
        termgraph.graph(labels, histo)



def showPropertiesOld(nkG):

    propertiesTextBlock = """
    Graph Properties: {name}
    ========================

    Basic Properties
    ----------------
    - nodes (n){n:16}
    - edges (m)            {m}
    - min. degree             {minDeg}
    - max. degree             {maxDeg}
    - isolated nodes         {isolates}
    - self-loops            {loops}
    - density            {dens:10.6f}


    Path Structure
    --------------
    - connected components         {components}
    - diameter            {dia}

    Miscellaneous
    -------------
    - degree assortativity            {assort:10.6f}
    - cliques                {cliques}

    Community Structure
    -------------------
    - avg. local clustering coefficient         {avglcc:10.6f}
    - PLP community detection
        - communities                {ncomPLP}
        - modularity             {modPLP:10.6f}
    - PLM community detection
        - communities            {ncomPLM}
        - modularity             {modPLM:10.6f}


    Distributions
    -------------

    - degree distribution

    """

    props = properties(nkG)
    print(textwrap.dedent(propertiesTextBlock.format(**props)))
    (labels, histo) = props["histo"]
    termgraph.graph(labels, histo)

    

def compressHistogram(hist, nbins=20):
    """ Compress a histogram to a number of bins"""
    compressed = [None for i in range(nbins)]
    labels = [None for i in range(nbins)]
    nbinsprev = len(hist)
    binwidth = math.ceil(nbinsprev / nbins)
    for i in range(nbins):
        l = i * binwidth
        u = (i+1) * binwidth
        compressed[i] = sum(hist[l : u])
        labels[i] = "{0}-".format(l)
    return (labels, compressed)
        

    
def printDegreeHistogram(nxG, nbins=25):
    """ Prints a degree histogram as a bar chart to the terminal"""
    hist = nx.degree_histogram(nxG)
    
    (labels, hist) = compressHistogram(hist, nbins)
    termgraph.graph(labels, hist)
    


def hpProperties(nkG):
    """ For large graphs: get an overview of some properties"""
    print("min/max degree:")
    
    
    
# NetworKit algorithm engineering workflows

# class DynamicCommunityDetectionWorkflow:
#     
#     def __init__(self):
#         clusterings = []    # list of clusterings
#     
#     
#     def start(self, nMax, deltaT):
#         
#         self.G = Graph(0)
#         self.Gproxy = GraphEventProxy(self.G)
#         #self.generator = DynamicBarabasiAlbertGenerator(self.Gproxy)
#         self.dcd = DynamicLabelPropagation(self.Gproxy)
#         
#         while (self.G.numberOfNodes() < nMax):
#             self.generator.generate()
#             if (self.G.time() % deltaT) == 0:
#                 zeta = self.dcd.run()
#                 self.clusterings.append(zeta)
    
class GraphConverter:
    
    def __init__(self, reader, writer):
        self.reader = reader
        self.writer = writer
        
    def convert(self, inPath, outPath):
        G = self.reader.read(inPath)
        self.writer.write(G, outPath)
        
    def __str__(self):
        return "GraphConverter: {0} => {0}".format(self.reader, self.writer)

def getConverter(fromFormat, toFormat):
    
    readers = {"metis": METISGraphReader, "edgelist" : EdgeListIO}    
    writers = {"edgelist": EdgeListIO}
    
    reader = readers[fromFormat]()
    writer = writers[toFormat]()
    
    return GraphConverter(reader, writer)



def inspectCommunities(zeta, G):
    """ Show information about communities"""
    communitySizes = zeta.clusterSizes()
    mod = Modularity().getQuality(zeta, G)
    commProps = [
        ["# communities", zeta.numberOfClusters()],
        ["min community size", min(communitySizes)],
        ["max community size", max(communitySizes)],
        ["avg. community size", sum(communitySizes) / len(communitySizes)],
        ["imbalance", zeta.getImbalance()],
        ["modularity", mod],
    ]
    print(tabulate.tabulate(commProps))
    

def evalCommunityDetection(algo, G):
    """ Evaluate a community detection algorithm """
    
    t = stopwatch.Timer()
    zeta = algo.run(G)
    t.stop()
    results = [
        ["time [s]", t.elapsed],
        ["# communities", zeta.numberOfClusters()],
        ["modularity", Modularity().getQuality(zeta, G)]
    ]
    print(tabulate.tabulate(results))
    
    
    
