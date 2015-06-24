# standard library modules
from base64 import b64encode
import logging


# networkit submodules
from . import community
from . import centrality
from . import termgraph
from . import auxiliary
from . import nxadapter
from . import stopwatch

# TODO: refactor imports
from . import properties

# external modules
try:
	import tabulate
except ImportError:
	logging.warning("""WARNING: module 'tabulate' not installed, which is required by some functions.""")

try:
	import pandas
except ImportError:
	logging.warning("""WARNING: module 'pandas' not installed, which is required by some functions.""")

try:
	import powerlaw
except ImportError:
	logging.warning("""WARNING: module 'powerlaw' not installed, which is required by some
						functions.""")

try:
	import matplotlib.pyplot as plt
	from matplotlib._pylab_helpers import Gcf
except ImportError:
	logging.warning("""WARNING: module 'matplotlib' not installed, which is required for plotting.""")

try:
	import seaborn
except ImportError:
	logging.warning("""WARNING: module 'seaborn' not installed, which is required by some plotting functions. Also your plots will look better with it.""")

try:
	from IPython.core.pylabtools import print_figure
	from IPython.core.display import HTML
except ImportError:
	logging.warning("""WARNING: module 'IPython' not installed, which is required by some functions.""")




def asImage(plotFunction, plotArgs=[], plotKwargs={}, size=(8,6)):
	"""
	Call any plot function with the given argument and return the image in an HTML <img> tag.
	"""
	plt.figure(figsize=size)
	plotFunction(*plotArgs, **plotKwargs)
	# Get a handle for the plot that was just generated
	fig = Gcf.get_all_fig_managers()[-1].canvas.figure
	# Generate a data URL for the image
	imageData = "data:image/png;base64,{0}".format(b64encode(print_figure(fig)).decode("utf-8"))
	# Remove the plot from the list of plots for the current cell
	Gcf.destroy_fig(fig)
	# generate img tag
	image = "<img src='{0}'\>".format(imageData)
	return image

def computeNetworkProperties(G):
	"""
	"""
	networkProperties = [
			["nodes, edges", "{0}, {1}".format(G.numberOfNodes(), G.numberOfEdges())],
			["directed?", "{0}".format(G.isDirected())],
			["weighted?", "{0}".format(G.isWeighted())],
			["density", "{0:.6f}".format(properties.density(G))],
			["diameter range", "{0}".format(properties.Diameter.estimatedDiameterRange(G, error=0.1))]
		]
	return networkProperties


def powerLawStats(centralities):
	powerLawStats = {}
	for (centralityName, centralityScores) in centralities.items():
		fit = powerlaw.Fit(centralityScores)
		R, p = fit.distribution_compare("power_law", "exponential", normalized_ratio=True)
		gamma = fit.alpha
		powerLawStats[centralityName] = ((R > 0), R, gamma)
	return powerLawStats

def computeRankCorrelations(centralities : pandas.DataFrame, method="spearman"):
	return centralities.corr(method=method)






def profile(G):
	"""
	Output profile page of network as HTML
	"""
	# settings
	defaultSize = (5,2)
	histArgs = {"bins" : 100, "figsize" : (12,8)}

	# compute global network attributes
	networkProperties = computeNetworkProperties(G)
	networkPropertiesTable = tabulate.tabulate(networkProperties, tablefmt="html")

	hopPlot = asImage(plot.hopPlot, plotArgs=[G], size=defaultSize)

	# compute node properties
	nodeProperties = perties(G)


	# compute figures
	(plaw, _, gamma) = properties.degreePowerLaw(G)

	# compute images
	ddPlot  = asImage(plot.degreeDistribution, plotArgs=[G], size=defaultSize)
	ddHist = asImage(nodeProperties["degree"].hist, plotKwargs=histArgs, size=defaultSize)
	ccPlot = asImage(plot.nodeProperty, plotKwargs={"data" : nodeProperties["clustering"], "label": "local clustering coefficient"}, size=defaultSize)
	ccHist = asImage(nodeProperties["clustering"].hist, plotKwargs=histArgs, size=defaultSize)
	compPlot = asImage(plot.connectedComponentsSizes, [G], size=(1,1))

	page = HTML(profileTemplate.format(**locals()))
	return page


class Profile:
	""" This class computes and presents a structural profile of a networks"""

	pageTemplate = """
		<style media="screen" type="text/css">
		#wrapper {
		    width: 500px;
		    border: 1px solid black;
		}
		#first {
		    width: 300px;
		    border: 1px solid red;
		}
		#second {
		    border: 1px solid green;
		}
		</style>

		<div id="page">
		<h1>Network Profile</h1>

		<h2>Network Properties</h2>
			<div id="wrapper">
				<div id="first">
					{networkPropertiesTable}
				</div>
				<div id="second">
					{hopPlot}
				</div>
			</div>

		<h2>Network Partitions</h2>

		<h2>Node Centrality Measures</h2>

		<h3>Degree</h3>
		{ddPlot}
		{ddHist}

		<h3>Local Clustering Coefficient</h3>
		{ccPlot}
		{ccHist}




		power law distribution: {plaw} {gamma}

		{compPlot}
		</div>
	"""


	def __init__(self, G, settings={}):
		if (G.isDirected()):
			raise Exception("Profiling currently only supported for undirected graphs")
		self.G = G
		self.settings = settings

	def computePartitions(self):
		G = self.G

		# TODO: refactor module membership of component algorithms
		partitionAlgos = 	{ 	"components":	(properties.ConnectedComponents, 	(G,)),
								"communities":	(community.PLM, 					(G,))
							}

		partitions = {}
		for (algoName, (algoClass, params)) in partitionAlgos.items():
			algo = algoClass(*params)
			t = stopwatch.Timer()
			algo.run()
			logging.info("{algoName} computed in {time} s".format(algoName=algoName, time="{:.2E}".format(t.elapsed)))
			partitions[algoName] = algo.getPartition()

		# store partition vectors in data frame so pandas can operate on them
		self.partitionVectors = pandas.DataFrame(dict((name, partition.getVector()) for (name, partition) in partitions.items()))


	def computeNodeCentralities(self):
		G = self.G
		(n, m) = G.size()

		# TODO: normalization?

		nodeCentralityAlgos = {
								"degree":		(centrality.DegreeCentrality, 			(G, )),
								"coreness":		(centrality.CoreDecomposition, 			(G, )),
								"clustering":	(centrality.LocalClusteringCoefficient, (G, )),
								"pagerank":		(centrality.PageRank, 					(G, )),
								"kpath":		(centrality.KPathCentrality,			(G, )),
								"katz":			(centrality.KatzCentrality,				(G, )),
								"betweenness":	(centrality.ApproxBetweenness2,			(G, max(42, n / 1000), False)),
								"closeness":	(centrality.ApproxCloseness,			(G, max(42, n / 1000), False))
								}

		centralityScores = {}
		for (algoName, (algoClass, params)) in nodeCentralityAlgos.items():
			algo = algoClass(*params)
			t = stopwatch.Timer()
			algo.run()
			logging.info("{algoName} computed in {time} s".format(algoName=algoName, time="{:.2E}".format(t.elapsed)))
			centralityScores[algoName] = algo.scores()
		self.nodeCentralities = pandas.DataFrame(centralityScores)



	def compute(self):
<<<<<<< local
		# compute node centralities
		self.computeNodeCentralities()
		# compute rank correlations between node centrality values
		self.nodeCentralityCorrelations = self.nodeCentralities.corr(method="spearman")
		self.describeNodeCentralities()

	def describeNodeCentralities(self):
		stdStats = self.nodeCentralities.describe()

		powerLawStats = pandas.DataFrame(columns=self.nodeCentralities.columns)
		for (centralityName, centralityScores) in self.nodeCentralities.items():
			fit = powerlaw.Fit(centralityScores)
			R1, p1 = fit.distribution_compare("power_law", "exponential", normalized_ratio=True)#
			R2, p2 = fit.distribution_compare("power_law", "lognormal")
			gamma = fit.alpha
			powerLawStats[centralityName] = pandas.Series([(R1 > 0) and (R2 > 0), gamma], index=["fitspowerlaw", "powerlawexponent"])

		self.nodeCentralityStats = stdStats.append(powerLawStats)

=======
		# compute node centralities
		self.computeNodeCentralities()
		# compute rank correlations between node centrality values
		self.nodeCentralityCorrelations = self.nodeCentralities.corr(method="spearman")
>>>>>>> other

	def getPage(self):
		raise NotImplementedError("TODO")

	def getAttributeVector(self):
		raise NotImplementedError("TODO:")

	def plotNodeCentralityCorrelations(self, figsize=(6,6)):
		cmap = seaborn.diverging_palette(220, 20, as_cmap=True)
		f, ax = plt.subplots(figsize=figsize)
		logging.info("calculating correlation heatmap")
		seaborn.corrplot(self.nodeCentralities, cmap=cmap, method="spearman")
		f.tight_layout()

	def plotNodeCentralityHistograms(self, figsize=(12,8)):
		self.nodeCentralities.hist(bins=50, figsize=figsize)
