import tin
import numpy as np

#reading the dataset
dataset = np.loadtxt('pts1000c.dat', delimiter = ' ')
data = tin.TIN(dataset)
#Return the terrain in 3D
data.plotting()
#Return the elevation of x point. 
#Also, putting the parameter 'True', it shows the triangulation graph.
print(data.elevation([4,5], True))
#Return if the point is a water source. Receives a point. 
print(data.peak_point(dataset[450,:2]))
#Return if the point is a water pit. Receives a point.
print(data.anti_peak_point(dataset[450,:2]))
#Plot all the relative maximum points.
data.relative_max()
#Plot all the relative minimum points.
data.relative_min()
#Plot the dual graph over the triangulation.
data.graph()