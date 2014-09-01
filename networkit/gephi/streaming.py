from . import pyclient as _gephipyclient   # we want to hide these from the user
import urllib as _urllib

"""
This module provides methods to export data from NetworKit directly to Gephi using the
Gephi Streaming Plugin.
"""

class GephiStreamingClient:
    """
     A python singleton providing access to a running Master instance of the Gephi
     Streaming plugin
    """

    def __init__(self, url='http://localhost:8080/workspace0'):
        #Disabling Flushing means quite a good performance boost.
        self._pygephi = _gephipyclient.GephiClient(url, autoflush=10000)
        self.graphExported = False

    def _urlError(self, e):
        print("Could not connect to the gephi streaming plugin. Did you start the streaming master server in gephi?")

    def exportGraph(self, graph):
        """ Exports the given graph to gephi. No attributes or weights are exported.
        Please note that only graphs with indexed edges can be exported, so
        edges will be indexed if neccessary. """

        try:
            self._pygephi.clean()

            #Index edges if neccessary
            graph.indexEdges()

            nAttrs = {}

            for node in graph.nodes():
                self._pygephi.add_node(str(node), **nAttrs)

            for edge in graph.edges():
                edgeId = str(min(edge[0],edge[1])) + '-' + str(max(edge[0],edge[1]))# graph.edgeId(edge[0], edge[1])
                self._pygephi.add_edge(edgeId, edge[0], edge[1], False)

            self._pygephi.flush()
            self.graphExported = True
        except _urllib.error.URLError as e:
            self._urlError(e)

    def exportAdditionalEdge(self, u, v):
        """ Adds an edge in an already exported graph."""
        if self.graphExported != True:
            print("Error: Cannot add edges. Export Graph first!")
            return      
        try:
            edgeId = str(min(u,v)) + '-' + str(max(u,v))# graph.edgeId(u,v)
            self._pygephi.add_edge(edgeId, u, v, False)
            self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e)        

    def removeExportedEdge(self, u, v):
        """ Removes an edge from an already exported graph."""
        if self.graphExported != True:
            print("Error: Cannot remove edges. Export Graph first!")
            return
        try:
            edgeId = str(min(u,v)) + '-' + str(max(u,v))#  graph.edgeId(u,v)
            self._pygephi.delete_edge(edgeId)
            self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e)

    def exportEventStream(self, stream):
        if self.graphExported != True:
            print("Error: Cannot export event stream. Export Graph first!")
            return

        nAttrs = {}
        try:
            for ev in stream:
                if ev.type == ev.NODE_ADDITION:
                    self._pygephi.add_node(str(ev.u), **nAttrs)
                elif ev.type == ev.NODE_REMOVAL:
                    self._pygephi.delete_node(str(ev.u))
                elif ev.type == ev.EDGE_ADDITION:
                    edgeId = str(min(ev.u,ev.v)) + '-' + str(max(ev.u,ev.v))
                    self._pygephi.add_edge(edgeId, ev.u, ev.v, False)
                elif ev.type == ev.EDGE_REMOVAL:
                    edgeId = str(min(ev.u,ev.v)) + '-' + str(max(ev.u,ev.v))
                    self._pygephi.delete_edge(edgeId)
                elif ev.type == ev.EDGE_WEIGHT_UPDATE:
                    print("Edge weights not yet supported in gephi streaming!")
                elif ev.type == ev.TIME_STEP:
                    self._pygephi.flush()
                self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e)

    def exportNodeValues(self, graph, values, attribute_name):
        """
        This method exports node values (e.g. community information, betwenness centrality values)
        to gephi using the Gephi Streaming Plugin. Use exportGraph(Graph) first to export
        the graph itself.

        Parameters:
        - values: python list or Partition object that contains the values to be exported.
        - attribute_name: name of the node attribute in gephi
        """

        try:
            if len(values) != graph.numberOfNodes():
                print("Warning: #Nodes (", graph.numberOfNodes(), ") does not match #Values (", len(values), ").")

            for i in range(0, min(len(values), graph.numberOfNodes())):
                nAttrs = {attribute_name:values[i]}
                self._pygephi.change_node(str(graph.nodes()[i]), **nAttrs)

            self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e)

    def exportCoordinates(self, graph, scale=1):
        try:
            xcoords = [scale*graph.getCoordinate(v)[0] for v in graph.nodes()]
            ycoords = [scale*graph.getCoordinate(v)[1] for v in graph.nodes()]
            self.exportNodeValues(graph, xcoords, 'x')
            self.exportNodeValues(graph, ycoords, 'y')
            self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e) 

    def exportEdgeValues(self, graph, values, attribute_name):
        """
        This method exports an edge attribute to gephi using the Gephi Streaming Plugin.
        Use exportGraph(Graph) first to export graph itself.

        Parameters:
        - values: python list contains the values to be exported.
        - attribute_name: name of the edge attribute in gephi
        """
        try:
            if len(values) != graph.numberOfEdges():
                print("Warning: #Nodes (", graph.numberOfEdges(), ") does not match #Values (", len(values), ").")

            idx = 0
            for edge in graph.edges():
                edgeId = str(min(edge[0],edge[1])) + '-' + str(max(edge[0],edge[1])) #graph.edgeId(edge[0], edge[1]) #
                eAttrs = {attribute_name:values[graph.edgeId(edge[0], edge[1])], "Type":"Undirected"}#still need to use the old edge to access the graph array
                self._pygephi.change_edge(edgeId, edge[0], edge[1], False, **eAttrs)
                idx += 1

            self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e)

    def clearGraph(self):
        """ Cleans the first workspace in Gephi. """

        try:
            self._pygephi.clean()
            self._pygephi.flush()
            self.graphExported = False

        except _urllib.error.URLError as e:
            self._urlError(e)
