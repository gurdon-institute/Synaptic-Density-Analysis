# Density and co-localisation of pre- and post-synaptic markers.
# by Richard Butler, Gurdon Institute Imaging Facility, University of Cambridge
#
# Copyright 2021, 2022 Richard Butler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>


import math as maths
import itertools

from ij import IJ, ImagePlus, ImageStack
from ij.process import ImageProcessor, Blitter, StackStatistics, AutoThresholder
from ij.plugin import Duplicator, GaussianBlur3D, Filters3D
from ij.plugin.filter import ThresholdToSelection
from ij.gui import Overlay, ShapeRoi, PointRoi
from ij.measure import ResultsTable

from java.awt import Color, Rectangle


factor = 2*maths.sqrt(2*maths.log(2)) # factor to get standard deviation from the FWHM of a Gaussian function
r2p = maths.sqrt(2*maths.pi)

dotW = 0.2	# µm, the estimated width of the dots to be detected
sigma = dotW/factor

kdsigma = dotW*2	# µm, the standard deviation of the density kernel
kdr = kdsigma*factor


class Dot:	# stores voxel and calibrated coordinates for a dot
	def __init__(self, x,y,z):
		self.x = x
		self.y = y
		self.z = z
		self.xCal = self.x * cal.pixelWidth
		self.yCal = self.y * cal.pixelHeight
		self.zCal = self.z * cal.pixelDepth
		
	def distance(self,other):
		dist = maths.sqrt( 	((self.xCal-other.xCal)*(self.xCal-other.xCal)) +
							((self.yCal-other.yCal)*(self.yCal-other.yCal)) +
							((self.zCal-other.zCal)*(self.zCal-other.zCal)) )
		return dist

	def equals(self,other):
		return self.x==other.x and self.y==other.y and self.z==other.z


def LoG3D(image, c, sigma):	# returns the 3D Laplacian of Gaussian of image channel c approximated by the difference of Gaussians with K = 1.4
	dup = Duplicator()
	log = dup.run(image, c,c, 1,imp.getNSlices(), 1,1)
	sub = dup.run(image, c,c, 1,imp.getNSlices(), 1,1)

	sigmaXY = sigma/cal.pixelWidth	# anisotropic standard deviations from voxel dimensions
	sigmaZ = sigma/cal.pixelDepth
	
	GaussianBlur3D.blur(log, sigmaXY, sigmaXY, sigmaZ)
	GaussianBlur3D.blur(sub, sigmaXY*1.4, sigmaXY*1.4, sigmaZ*1.4)

	stack = ImageStack(bounds.width, bounds.height)
	for z in range(1,image.getNSlices()+1):
		ip0 = log.getStack().getProcessor(z)
		ip1 = sub.getStack().getProcessor(z)
		ip0.copyBits(ip1, 0,0, Blitter.SUBTRACT)
		stack.addSlice(ip0)
	
	return ImagePlus(imp.getTitle()+"_LoG", stack)


def getDots3D(image, threshold): # returns a list of Dots from a 1 channel z-stack
	stack = image.getStack()
	vxy = dotW/cal.pixelWidth
	vz = dotW/cal.pixelDepth
	maxStack = Filters3D.filter(stack, Filters3D.MAX, vxy, vxy, vz)
	dots = []
	for z in range(stack.getSize()):
		for y in range(stack.getHeight()):
			for x in range(stack.getWidth()):
				if stack.getVoxel(x,y,z) == maxStack.getVoxel(x,y,z) >= threshold: # this is a local maximum
					dots.append( Dot(bounds.x+x, bounds.y+y, z) )

	keep = [dot for dot in dots]
	for p in dots:
		edx = min(p.x-bounds.x,bounds.x+bounds.width-p.x)*cal.pixelWidth
		edy = min(p.y-bounds.y,bounds.y+bounds.height-p.y)*cal.pixelHeight
		edz = min(p.z,D-p.z)*cal.pixelDepth
		if edx < dotW or edy < dotW or edz < dotW:	# this dot is on an edge
			if p in keep:
				ded = PointRoi(p.x, p.y)
				ded.setStrokeColor(Color.RED)
				ded.setPosition(0,int(p.z)+1,0)
				ol.add(ded)
				keep.remove(p)

	for p,q in itertools.combinations(dots,2):
		if p.distance(q) <= dotW:
			merge = Dot((p.x+q.x)/2., (p.y+q.y)/2., (p.z+q.z)/2.)
			if p in keep:
				keep.remove(p)
			if q in keep:
				keep.remove(q)
			keep.append(merge)

	return keep


def getThreshold(image, method):	# returns the automatic threshold from a 1 channel z-stack using the specified AutoThresholder.Method
	stats = StackStatistics(image)
	hist = [int(v) for v in stats.getHistogram()]
	thresh = AutoThresholder().getThreshold(method, hist)
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
	return thresh

# get the current ImagePlus and area to be analysed (projected through Z)
imp = IJ.getImage()
userRoi = imp.getRoi()
if userRoi is not None:
	bounds = userRoi.getBounds()
else:
	bounds = Rectangle(0,0, imp.getWidth(),imp.getHeight())
D = imp.getNSlices()
cal = imp.getCalibration()
ol = Overlay()

procC2 = LoG3D(imp, 2, sigma)
procC3 = LoG3D(imp, 3, sigma)

# find Dots and add them to the Overlay
threshC2 = getThreshold(procC2, AutoThresholder.Method.Otsu)
dotsC2 = getDots3D(procC2, threshC2)
for dot in dotsC2:
	roi = PointRoi(dot.x, dot.y)
	roi.setStrokeColor(Color.CYAN)
	roi.setPosition(0, int(dot.z)+1, 0)
	ol.add(roi)

threshC3 = getThreshold(procC3, AutoThresholder.Method.Otsu)
dotsC3 = getDots3D(procC3, threshC3)
for dot in dotsC3:
	roi = PointRoi(dot.x, dot.y)
	roi.setStrokeColor(Color.MAGENTA)
	roi.setPosition(0, int(dot.z)+1, 0)
	ol.add(roi)

# output results for this image in a new table
rt = ResultsTable()
rt.showRowNumbers(False)
rt.setPrecision(6)
totalV = bounds.width*cal.pixelWidth * bounds.height*cal.pixelHeight * D*cal.pixelDepth
ncoloc = 0
for i,dc2 in enumerate(dotsC2):

	rt.setValue("C2 Dot", i, i)
	rt.setValue("X", i, dc2.xCal)
	rt.setValue("Y", i, dc2.yCal)
	rt.setValue("Z", i, dc2.zCal)

	C2kde = 0.
	C2minD = float("inf")
	for dc2b in dotsC2:
		if dc2b is dc2:
			continue
		C2dist = dc2.distance(dc2b)
		C2kde += (1/(kdsigma*r2p))*maths.exp( -(C2dist*C2dist)/(2*kdsigma*kdsigma) )
		C2minD = min(C2dist, C2minD)
	rt.setValue("Nearest C2 Dot Distance", i, C2minD)
	d2 = (C2kde/len(dotsC2))
	rt.setValue("C2 Density", i, d2)

	C3count = 0
	C3kde = 0.
	C3minD = float("inf")
	for dc3 in dotsC3:
		C3dist = dc2.distance(dc3)
		if C3dist <= dotW:
			C3count += 1
			ncoloc += 1
		C3kde += (1/(kdsigma*r2p))*maths.exp( -(C3dist*C3dist)/(2*kdsigma*kdsigma) )
		C3minD = min(C3dist, C3minD)
	rt.setValue("n C3 Dots <= "+str(dotW)+" "+u'\u00b5'+"m", i, C3count)
	rt.setValue("Nearest C3 Dot Distance", i, C3minD)
	d3 = (C3kde/len(dotsC3))
	rt.setValue("C3 Density", i, d3)

rt.show(imp.getTitle()+"_Dots")
imp.setOverlay(ol)

# add a row for this image to the main ImageJ ResultsTable
mainrt = ResultsTable.getResultsTable()
row = mainrt.getCounter()
mainrt.setValue("Image", row, imp.getTitle())
mainrt.setValue("Volume ("+u"\u00b5"+"m"+u"\u00b3"+")", row, totalV)
mainrt.setValue("C2 Count", row, len(dotsC2))
mainrt.setValue("C3 Count", row, len(dotsC3))
mainrt.setValue("Colocalised Count", row, ncoloc)
mainrt.setValue("C2 Density", row, len(dotsC2)/totalV)
mainrt.setValue("C3 Density", row, len(dotsC3)/totalV)
mainrt.show("Results")
