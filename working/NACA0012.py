import numpy
from matplotlib import pyplot

import pblocks

N = 51                                   # number of points in each direction
x_start, x_end = -1.0, 2.0               # boundaries in the x-direction
y_start, y_end = -0.5, 0.5               # boundaries in the y-direction
x = numpy.linspace(x_start, x_end, N)    # creates a 1D-array with the x-coordinates
y = numpy.linspace(y_start, y_end, N)    # creates a 1D-array with the y-coordinates
X, Y = numpy.meshgrid(x, y)              # generates a mesh grid


#####################################

x_source = numpy.loadtxt('resources\\NACA0012_x.txt')
y_source = numpy.loadtxt('resources\\NACA0012_y.txt')
strength_source = numpy.loadtxt('resources\\NACA0012_sigma.txt')

U = 1.0             # freestream

# computes stream-function and velocity field
psi, u, v = 0, 0, 0
for i in range(len(x_source)):
    psi += pblocks.source_get_stream_function(strength_source[i], x_source[i], y_source[i], X, Y)
    u_s, v_s = pblocks.source_get_velocity(strength_source[i], x_source[i], y_source[i], X, Y)
    u += u_s
    v += v_s
    
psi += pblocks.freestream_get_stream_function(U, 0., X, Y)
u_stream, v_stream = pblocks.freestream_get_velocity(U, 0., X, Y)
u += u_stream
v += v_stream

# computes the pressure coefficient field
cp = 1.0 - (u**2 + v**2)/U**2

##############################
# plotting the streamlines
size = 10
pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)

index = numpy.argmax(cp)
ix, iy = numpy.unravel_index(index, [N, N])
print(ix, iy)
print(cp[ix, iy])

## plots the stremlines
#pyplot.streamplot(X, Y, u, v, density=2, linewidth=1, arrowsize=2, arrowstyle='->')
#pyplot.scatter(x_source, y_source, color='#CD2305', s=10, marker='o', linewidth=0)
#pyplot.plot(x_source, y_source, linewidth=2, color='k')

# plots the pressure coefficient field
contf = pyplot.contourf(X, Y, cp, levels=numpy.linspace(-0.5, 0.5, 50), extend='both')
cbar = pyplot.colorbar(contf)
cbar.set_label('$C_p$', fontsize=16)
cbar.set_ticks([-0.5, 0.0, 0.5])
pyplot.plot(x_source, y_source, linewidth=1, color='k')
pyplot.scatter(X[0, iy], Y[ix, 0], color='k', s=10, marker='o', linewidth=0)

# adds the dividing line to the figure
#pyplot.contour(X, Y, psi, 
#            levels=[-strength_source/3, +strength_source/3], 
#            colors='#CD2305', linewidths=2, linestyles='solid')
    