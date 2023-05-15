from FlightOperations import *

"""
Read the dataset.
"""

border = 'borders_CO.dat'
air = 'airports_CO.dat'
Col = FlightOperations(border, air)

"""
Plots the original map and its voronoi diagram
"""

Col.plot_original()

"""
Plots the map with the color scale that depends on thealtitude of each airport
"""

Col.plotting_altitude()

"""
Returns the center and the coordinates of the new airport.
If the input argument is True then it plots the new airport.
The default input is False.
"""

print(Col.find_new_airport(True))

"""
Returns the nearest airports.
"""

print(Col.nearest_airports())

"""
Returns the farthest airports.
"""

print(Col.farthest_airports())

"""
Returns the path from the airport with index 8 to airport with index 9.
The threshold is defined by 100.
If the input argument is True then it plots the path.
"""

print(Col.path(100, 8, 9, True))
