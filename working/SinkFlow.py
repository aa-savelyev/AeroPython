import numpy
from matplotlib import pyplot

import pblocks

N = 200                                  # number of points in each direction
x_start, x_end = -4.0, 4.0               # boundaries in the x-direction
y_start, y_end = -2.0, 2.0               # boundaries in the y-direction
x = numpy.linspace(x_start, x_end, N)    # creates a 1D-array with the x-coordinates
y = numpy.linspace(y_start, y_end, N)    # creates a 1D-array with the y-coordinates
X, Y = numpy.meshgrid(x, y)              # generates a mesh grid


#####################################

# plots the grid of points
#size = 10
#pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
#pyplot.xlabel('x', fontsize=16)
#pyplot.ylabel('y', fontsize=16)
#pyplot.xlim(x_start, x_end)
#pyplot.ylim(y_start, y_end)
#pyplot.scatter(X, Y, s=10, color='#CD2305', marker='o', linewidth=0)

strength_source = 5.0                      # source strength
x_source, y_source = -1.0, 0.0             # location of the source

# computes the velocity field
u_source, v_source = pblocks.source_get_velocity(strength_source, x_source, y_source, X, Y)

# computes the stream-function
psi_source = pblocks.source_get_stream_function(strength_source, x_source, y_source, X, Y)

# plotting the streamlines
size = 10
pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.streamplot(X, Y, u_source, v_source, density=2, linewidth=1, arrowsize=2, arrowstyle='->')
pyplot.scatter(x_source, y_source, color='#CD2305', s=80, marker='o', linewidth=0)

# adds the dividing line to the figure
pyplot.contour(X, Y, psi_source, 
            levels=[-strength_source/3, +strength_source/3], 
            colors='#CD2305', linewidths=2, linestyles='solid')