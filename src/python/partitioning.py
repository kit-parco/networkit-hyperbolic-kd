import math
import numpy as np

from algebraic import laplacianEigenvectors
from _NetworKit import Partition, Modularity

try:
	import tabulate
except ImportError:
	print(""" WARNING: module 'tabulate' not found, please install it to use the full functionality of NetworKit """)

def computeEdgeCut(partition, graph):
	cut = 0

	for (n1, n2) in graph.edges():
		if partition[n1] != partition[n2]:
			if (graph.isWeighted()):
				cut += graph.weight(n1,n2)
			else:
				cut += 1

	return cut

def computeImbalance(partition, graph):
	desired = graph.numberOfNodes() / float(partition.numberOfSubsets())

	maximum = max(partition.subsetSizes())

	return maximum / desired

def inspectPartitions(partition, graph):
	partitionSizes = partition.subsetSizes()
	mod = Modularity().getQuality(partition, graph)
	props = [
		["# partitions", partition.numberOfSubsets()],
		["min partition size", min(partitionSizes)],
		["max partition size", max(partitionSizes)],
		["avg. partition size", sum(partitionSizes) / len(partitionSizes)],
		["imbalance", computeImbalance(partition, graph)],
		["edge cut", computeEdgeCut(partition, graph)],
		["modularity", mod],
	]
	print(tabulate.tabulate(props))

class SpectralPartitioner(object):
	"""

	Please note that the code in this class assumes the nodes of a graph to be numbered
	from 0 to n.

	"""
	def __init__(self, graph, count, balanced=True):
		self.graph = graph
		self.count = count

		self.balanced = balanced

	def prepareSpectrum(self):
		spectrum = laplacianEigenvectors(self.graph, cutoff = (math.ceil(math.log(self.count, 2)) + 1))
		self.eigenvectors = spectrum[1]
		self.eigenvalues = spectrum[0]

	def getQuantiles(self, eigv, vertices, count = 1):
		values = [eigv[i] for i in vertices]
		values.sort()

		sections = count + 1
		quantiles = []

		for i in range(1, sections):
			quantile = values[math.floor(len(values) * i / sections)]
			quantiles.append(quantile)

		return quantiles

	def getMean(self, eigv, vertices):
		values = [eigv[i] for i in vertices]
		mean = np.mean(values)

		return mean

	def trisect(self, partition=None, iteration=1):
		if partition is None:
			vertices = self.graph.nodes()
		else:
			vertices = self.partitions[partition]


		eigv = self.eigenvectors[iteration]

		quantiles = self.getQuantiles(eigv, vertices, count = 2)


		partA = self.nextPartition
		partB = self.nextPartition + 1
		partC = self.nextPartition + 2
		self.nextPartition += 3		# TODO this is not thread-safe


		self.partitions[partA] = []
		self.partitions[partB] = []
		self.partitions[partC] = []

		for vertex in vertices:
			if (eigv[vertex] < quantiles[0]):
				self.partitions[partA].append(vertex)
			elif (eigv[vertex] < quantiles[1]):
				self.partitions[partB].append(vertex)
			else: 
				self.partitions[partC].append(vertex)

		if (not (partition is None)):
			del self.partitions[partition]

	def bisect(self, count, partition=None, iteration=1):
		if count == 1:
			return

		if count == 3:
			self.trisect()
			return

		if partition is None:
			vertices = self.graph.nodes()
		else:
			vertices = self.partitions[partition]


		eigv = self.eigenvectors[iteration]

		if (self.balanced):
			split = self.getQuantiles(eigv, vertices)[0]
		else:
			split = self.getMean(eigv, vertices)

		partA = self.nextPartition
		partB = self.nextPartition + 1
		self.nextPartition += 2		# TODO this is not thread-safe


		self.partitions[partA] = []
		self.partitions[partB] = []

		for vertex in vertices:
			if (eigv[vertex] < split):
				self.partitions[partA].append(vertex)
			else:
				self.partitions[partB].append(vertex)

		if (not (partition is None)):
			del self.partitions[partition]

		if count > 2:
			if (count % 2 == 0):
				self.bisect(count / 2, partition = partA, iteration = iteration + 1)
				self.bisect(count / 2, partition = partB, iteration = iteration + 1)
			else:
				nextCount = (count - 1) / 2
				if nextCount > 2:
					self.bisect(nextCount, partition = partA, iteration = iteration + 1)
					self.bisect(nextCount + 1, partition = partB, iteration = iteration + 1)
				else:
					self.bisect(nextCount, partition = partA, iteration = iteration + 1)
					self.trisect(partition = partB, iteration = iteration + 1)


	def generatePartition(self):
		partition = Partition(size=self.graph.numberOfNodes())

		for partIndex in self.partitions:
			vertices = self.partitions[partIndex]

			if (len(vertices) < 1):
				continue

			firstItem = vertices[0] # TODO come on.. there has to be an easier way of doing this...
			partition.toSingleton(firstItem)
			subsetID = partition[firstItem]

			for vertex in vertices[1:]:
				partition.addToSubset(subsetID, vertex)

		self.partition = partition
		return partition

	def run(self):
		self.nextPartition = 0
		self.partitions = {}
		self.prepareSpectrum()

		self.bisect(self.count)

		self.generatePartition()

	def getPartition(self):
		return self.partition