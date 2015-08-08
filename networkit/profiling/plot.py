#
# file: threadpool.py
# author: Mark Erb
#

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import io
from urllib.parse import quote
import math


class Theme:
	""" TODO: """
	def __init__(self):
		self.set()
    
	
	@classmethod
	def RGBA2RGB(cls, color, alpha, background):
		result = (
			color[0] * alpha + background[0] * (1-alpha),
			color[1] * alpha + background[1] * (1-alpha),
			color[2] * alpha + background[2] * (1-alpha),
			1
		)
		return result		
	
	
	def set(self, style="light", color=(0, 0, 1)):
		optionsStyle = ["light", "system"]
		if style not in optionsStyle:
			raise ValueError("possible style options: " + str(optionsStyle))
		if len(color) != 3:
			raise ValueError("(r,g,b) tuple required")
				
		if style == "system":
			self.__rcParams = mpl.rcParams
			raise ValueError("not implemented, yet")
			
		if style == "light":
			self.__defaultColor = (0, 0, 0)
			self.__defaultWidth = 1
			self.__backgroundColor = (1, 1, 1)
			self.__plotColor = Theme.RGBA2RGB(color, 0.6, self.__backgroundColor)
			self.__plotWidth = 3
			self.__faceColor = (color[0], color[1], color[2], 0.2)
			self.__faceColorGray = "lightgray"
			self.__edgeColor = (color[0], color[1], color[2], 0.6)
			self.__edgeColorGray = (0, 0, 0)
			self.__edgeWidth = 2
			self.__gridColor = "lightgray"
			self.__fontColor = (0, 0, 0)
			self.__fontSize = 10

		self.__color = color
		self.__style = style
    
	
	def get(self):
		return (self.__style, self.__color)

		
	def getRcParams():
		return self.__rcParams

		
	def getDefaultColor(self):
		return self.__defaultColor
	def getDefaultWidth(self):
		return self.__defaultWidth
	def getPlotColor(self):
		return self.__plotColor
	def getPlotWidth(self):
		return self.__plotWidth
	def getFaceColor(self):
		return self.__faceColor
	def getFaceColorGray(self):
		return self.__faceColorGray
	def getEdgeColor(self):
		return self.__edgeColor
	def getEdgeColorGray(self):
		return self.__edgeColorGray
	def getEdgeWidth(self):
		return self.__edgeWidth
	def getBackgroundColor(self):
		return self.__backgroundColor
	def getGridColor(self):
		return self.__gridColor
	def getFontSize(self):
		return self.__fontSize
	def getFontColor(self):
		return self.__fontColor

		
class Measure:
	""" TODO: """
	def __init__(self, name, params):
		self.__name = name
		self.__params = params

	def getName(self):
		return self.__name

	def getType(self):
		return "Plot.Measure"

	def run(self):
		(index, stat, label, theme) = self.__params
		plt.ioff()

		
		def funcSpace(min, max):
			result = 0.1
			if min < max:
				result = (max - min) * 0.04
			return result


		def funcTicks(min, max, numberOfTicks):
			result = []
			if numberOfTicks > 0:
				value = min    
				step = (max - min) / numberOfTicks
				for i in range(numberOfTicks):
					result.append(min + step*i)
					value += step
				result.append(max)
			return result


		def funcPlotEnd(fig, ax, theme, width, height, drawAxis=True):
			ax.patch.set_facecolor(theme.getBackgroundColor())
			if drawAxis:
				axisColor = theme.getGridColor()
			else:
				axisColor = theme.getBackgroundColor()
			ax.spines["bottom"].set_color(axisColor)
			ax.spines["top"].set_color(axisColor) 
			ax.spines["right"].set_color(axisColor)
			ax.spines["left"].set_color(axisColor)
			ax.tick_params(axis="x", colors=theme.getGridColor(), which="both", labelsize=theme.getFontSize())
			ax.tick_params(axis="y", colors=theme.getGridColor(), which="both", labelsize=theme.getFontSize())
			ax.xaxis.label.set_color(theme.getFontColor())
			ax.yaxis.label.set_color(theme.getFontColor())
			[x_ticklabel.set_color(theme.getFontColor()) for x_ticklabel in ax.get_xticklabels()]
			[y_ticklabel.set_color(theme.getFontColor()) for y_ticklabel in ax.get_yticklabels()]
			fig.set_size_inches(width, height)
			

		def funcPlotBox(ax, x_numberOfTicks, x_showTickLabels, showGrid):
			q1 = stat["Location"]["1st Quartile"]
			q3 = stat["Location"]["3rd Quartile"]
			median = stat["Location"]["Median"]
			outlier_lower = stat["Location"]["Outlier (Lower)"][0]
			outlier_upper = stat["Location"]["Outlier (Upper)"][0]
			whisker_lower = stat["Location"]["Outlier (Lower)"][1]
			whisker_upper = stat["Location"]["Outlier (Upper)"][1]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			space = funcSpace(x_min, x_max)
			ticks = funcTicks(x_min, x_max, x_numberOfTicks)
			ax.scatter(
				[x_min, x_max],
				[0.5, 0.5],
				color = theme.getEdgeColor(),
				s = 35
			)
			ax.scatter(
				[outlier_lower, outlier_upper],
				[0.5, 0.5],
				color = theme.getDefaultColor(),
				marker = 'x',
				s = 50
			)
			ax.plot(
				[median, median],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_lower, whisker_lower],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_lower, q1],
				[0.5, 0.5],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_upper, whisker_upper],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[whisker_upper, q3],
				[0.5, 0.5],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.plot(
				[median, median],
				[0.2, 0.8],
				color = theme.getDefaultColor(),
				linestyle = '-',
				linewidth = theme.getDefaultWidth()
			)
			ax.add_patch(patches.Rectangle(
				(q1, 0.2),
				width = q3-q1,
				height = 0.6,
				facecolor = theme.getFaceColor(),
				linestyle = "solid",
				linewidth = theme.getEdgeWidth(),
				edgecolor = theme.getEdgeColor()
			))			
			ax.set_xlim([x_min-space, x_max+space])
			ax.set_ylim([0, 1])
			# ax.set_xticks(ticks)
			ax.set_yticks([])
			if not x_showTickLabels:
				ax.set_xticklabels([])
			if showGrid:
				ax.grid(showGrid, which="both", color=theme.getGridColor(), linestyle="-")
			return ax
			

		def funcPlotPDF(ax, x_numberOfTicks, y_numberOfTicks, x_showTickLabels, y_showTickLabels, showGrid):
			numberOfBins = stat["Binning"]["Number Histogram"]
			intervals = stat["Binning"]["Intervals Histogram"]
			absoluteFrequencies = stat["Binning"]["Absolute Frequencies Histogram"]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			y_min = 0
			y_max = stat["Binning"]["Mode"][1]
			x_space = funcSpace(x_min, x_max)
			y_space = funcSpace(y_min, y_max)
			x_ticks = funcTicks(x_min, x_max, x_numberOfTicks)
			for i in range(numberOfBins):
				ax.add_patch(patches.Rectangle(
					(intervals[i], 0),
					width = intervals[i+1]-intervals[i],
					height = absoluteFrequencies[i],
					facecolor = theme.getFaceColor(),
					linestyle = "solid",
					linewidth = theme.getEdgeWidth(),
					edgecolor = theme.getEdgeColor(),
				))
			ax.set_xlim([x_min-x_space, x_max+x_space])
			ax.set_ylim([0, y_max+y_space])
			# ax.set_xticks(x_ticks)
			if not x_showTickLabels:
				ax.set_xticklabels([])
			if not y_showTickLabels:
				ax.set_yticklabels([])
			if showGrid:
				ax.grid(showGrid, which="both", color=theme.getGridColor(), linestyle="-")
			return ax


		def funcPlotCDF(ax, x_numberOfTicks, y_numberOfTicks, x_showTickLabels, y_showTickLabels, showGrid):
			numberOfBins = stat["Binning"]["Number CDF"]
			intervals = stat["Binning"]["Intervals CDF"]
			comulativeRelativeFrequencies = stat["Binning"]["Relative Frequencies CDF"]
			x_min = stat["Location"]["Min"]
			x_max = stat["Location"]["Max"]
			x_space = funcSpace(x_min, x_max)
			y_space = funcSpace(0, 1)
			x_ticks = funcTicks(x_min, x_max, x_numberOfTicks)
			y_ticks = funcTicks(0, 1, y_numberOfTicks)
			for i in range(numberOfBins):
				ax.plot(
					[intervals[i], intervals[i+1]],
					[comulativeRelativeFrequencies[i], comulativeRelativeFrequencies[i]],
					color = theme.getPlotColor(),
					linestyle = "-",
					linewidth = theme.getEdgeWidth()
				)
				ax.plot(
					[intervals[i], intervals[i]],
					[0 if i==0 else comulativeRelativeFrequencies[i-1], comulativeRelativeFrequencies[i]],
					color = theme.getDefaultColor(),
					linestyle = "dotted",
				linewidth = theme.getDefaultWidth()
				)
			ax.set_xlim([x_min-x_space, x_max+x_space])
			ax.set_ylim([0, 1+y_space])
			# ax.set_xticks(x_ticks)
			ax.set_yticks(y_ticks)
			if not x_showTickLabels:
				ax.set_xticklabels([])
			if not y_showTickLabels:
				ax.set_yticklabels([])
			if showGrid:
				ax.grid(showGrid, which="both", color=theme.getGridColor(), linestyle="-")
			return ax
    
	
		def funcPlotPie(ax):
			numberOfTooSmallSubsets = stat["Binning"]["Pie"][1]
			relativeFrequencies = stat["Binning"]["Pie"][0]
			radius = 2
			accumulator = 0
			
			for i in range(len(relativeFrequencies)):
				value = relativeFrequencies[i]
				alpha = 360 * value
				label = "{:1.1f}%".format(value * 100)
				labelRadius = radius * 1.1
				if i == 0:
					label = str(numberOfTooSmallSubsets) + " Subsets\n" + label
				else:
					t = accumulator + alpha/2
					if (value < 0.022 and (
						t <= 30 and i%2 == 0 or
						t >= 150 and t <= 180 and i%2 == 1 or
						t >= 180 and t <= 210 and i%2 == 0 or
						t >= 330 and i%2 == 1
					)):
						labelRadius = radius * 1.19
				if accumulator + alpha/2 > 180:
					ha = "left"
				else:
					ha = "right"
				
				if i == 0:
					ax.add_patch(patches.Wedge(
						(math.cos(math.pi/180 * (90 + alpha/2)) * radius * 0.1, 
						math.sin(math.pi/180 * (90 + alpha/2)) * radius * 0.1),
						radius,
						90,
						90 + accumulator + alpha,
						facecolor = theme.getFaceColorGray(),
						edgecolor = theme.getDefaultColor()
					))
					labelRadius = radius * 1.18
				else:
					scale = 1-relativeFrequencies[1]/value
					faceColor = (
						theme.getPlotColor()[0] * scale,
						theme.getPlotColor()[1] * scale,
						theme.getPlotColor()[2] * scale,
						theme.getEdgeColor()[3]
					)
					ax.add_patch(patches.Wedge(
						(0, 0),
						radius,
						90 + accumulator,
						90 + accumulator + alpha,
						facecolor = faceColor,
						edgecolor = theme.getDefaultColor()
					))
				plt.text(
					math.cos(math.pi/180 * (90 + accumulator + alpha/2)) * labelRadius, 
					math.sin(math.pi/180 * (90 + accumulator + alpha/2)) * labelRadius,
					s = label,
					ha = ha,
					#family = 'sans-serif',
					size = theme.getFontSize()
				)
				accumulator += alpha
			
			ax.set_xlim([-3.2, 3.2])
			ax.set_ylim([-2.7, 2.7])
			ax.set_xticks([])
			ax.set_yticks([])
			return ax
	
	
		if index == 0:
			fig = plt.figure()
			
			ax1 = plt.subplot2grid((40, 8), (0, 0), colspan=8, rowspan=3)
			ax1.set_ylabel('Box')
			ax1.yaxis.set_label_position("right")
			funcPlotBox(
				ax = ax1,
				x_numberOfTicks = 5,
				x_showTickLabels = False,
				showGrid = False
			)
			funcPlotEnd(
				fig = fig,
				ax = ax1,
				theme = theme,
				width = 4,
				height = 0.2
			)
			
			ax2 = plt.subplot2grid((40, 8), (3, 0), colspan=8, rowspan=20)
			ax2.set_ylabel("PDF (absolute)")
			ax2.yaxis.set_label_position("right")
			funcPlotPDF(
				ax = ax2,
				x_numberOfTicks = 5,
				y_numberOfTicks = 5,
				x_showTickLabels = False,
				y_showTickLabels = True,
				showGrid = True
			)
			funcPlotEnd(
				fig = fig,
				ax = ax2,
				theme = theme,
				width = 4,
				height = 3
			)

			ax3 = plt.subplot2grid((40, 8), (23, 0), colspan=8, rowspan=17)
			ax3.set_xlabel(label)
			ax3.set_ylabel('CDF (relative)')
			ax3.yaxis.set_label_position("right")
			funcPlotCDF(
				ax = ax3,
				x_numberOfTicks = 5,
				y_numberOfTicks = 5,
				x_showTickLabels = True,
				y_showTickLabels = True,
				showGrid = True
			)
			funcPlotEnd(
				fig = fig,
				ax = ax3,
				theme = theme,
				width = 4,
				height = 3
			)
			fig.set_size_inches(6, 6)
		
		elif index == 1:
			fig = plt.figure()
			ax = fig.gca()
			
			ax.set_xlabel(label)
			# ax.xaxis.set_label_position("top")
			ax.set_ylabel("PDF (absolute)")
			ax.yaxis.set_label_position("right")
			funcPlotPDF(
				ax = ax,
				x_numberOfTicks = 5,
				y_numberOfTicks = 5,
				x_showTickLabels = True,
				y_showTickLabels = True,
				showGrid = True
			)
			funcPlotEnd(
				fig = fig,
				ax = ax,
				theme = theme,
				width = 4,
				height = 2.5
			)
		
		elif index == 2:
			fig = plt.figure()
			ax = fig.gca()
			
			funcPlotPie(
				ax = ax,
			)
			funcPlotEnd(
				fig = fig,
				ax = ax,
				theme = theme,
				width = 6.4*1.5,
				height = 5.4*1.5,
				drawAxis = False
			)
	
		fig.tight_layout()
		imgdata = io.StringIO()
		fig.savefig(imgdata, format='svg')
		plt.close(fig)
		plaintext = imgdata.getvalue()
		plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
		encoded = quote(plaintext, safe='');
		return (index, encoded)

		
class Scatter:
	""" TODO: """
	def __init__(self, name, params):
		self.__name = name
		self.__params = params

	def getName(self):
		return self.__name

	def getType(self):
		return "Plot.Scatter"

	def run(self):
		nameA = self.__name
		(nameB, sample_1, sample_2) = self.__params
		plt.ioff()

		
		def hexbin(ax, x, y, color, **kwargs):
			# cmap = sns.light_palette(color, as_cmap=True)
			ax.hexbin(x, y, gridsize=32, bins="log", **kwargs)
			# ax = plt.hexbin(x, y, gridsize=32, bins="log", cmap=cmap, **kwargs)
			ax.set_xlabel(nameA)
			# ax.xaxis.set_label_position("top")
			ax.set_ylabel(nameB)
			# ax.yaxis.set_label_position("right")
		
		
		fig = plt.figure()
		ax = fig.gca()	
		
		hexbin(ax, sample_1, sample_2, "#000070")
		
		fig.set_size_inches(4, 3.75)

		fig.tight_layout()
		imgdata = io.StringIO()
		fig.savefig(imgdata, format='svg')
		plt.close(fig)
		plaintext = imgdata.getvalue()
		plaintext = " ".join(plaintext[plaintext.find("<svg "):].split())
		encoded = quote(plaintext, safe='');
		return (nameB, encoded)