from AutonFleet import *

"""
Read the dataset.
"""
data = 'fleet100.dat'
x = AutonFleet(data)

"""
Returns and plots the s nearest ship according to the ship.
"""
s = 3
ship = 2
x.nearest_ship(s, ship)

"""
Plot the ships that are close to the ship in a radius of r.
"""
r = 0.35
x.avoiding_collisions(r, ship)

"""
Plots all the closest-pair of locations of ships.
"""
x.partner_lookup()

"""
Plots the approximate circle that contains all the fleet.
"""
x.operation_radius()
