import numpy
import math
from matplotlib import pyplot

import pblocks

def Jouk_transform(c, z):
    '''Joukowski transformation of z'''
    
    return z + c/z

dx = dy = 0.1
c = 1.0
#R = ((c-xc)**2 + yc**2)**0.5
R = 1.15
xc = -0.15
yc = 0.
j = (-1)**0.5

N_r = 100
N_th = 145
r_start, r_end = R, 5.0                  # boundaries in the r-direction
theta_start, theta_end = -0.0, 2*math.pi # boundaries in the y-direction
r = numpy.linspace(r_start, r_end, N_r)    # creates a 1D-array with the x-coordinates
theta = numpy.linspace(theta_start, theta_end, N_th)    # creates a 1D-array with the y-coordinates
R_mesh, Theta = numpy.meshgrid(r, theta) 

X = xc + R_mesh*numpy.cos(Theta)
Y = yc + R_mesh*numpy.sin(Theta)
Z = X + Y*j
#pyplot.scatter(X, Y, s=1)

Z2 = Jouk_transform(c, Z)
X2 = Z2.real
Y2 = Z2.imag
#pyplot.scatter(X2, Y2, s=1)

u_inf = 1.0
strength = R**2 * u_inf * 2*math.pi
print(strength)
u_stream, v_stream = pblocks.freestream_get_velocity(u_inf, 0., X, Y)
u_dublet, v_dublet = pblocks.doublet_get_velocity(strength, xc, yc, X, Y)

u = u_stream + u_dublet
v = v_stream + v_dublet
#u = u_dublet
#v = v_dublet
cp = 1 - (u**2 + v**2)/u_inf
U = u - v*j
U2 = U / (1 - (c/Z)**2)
u2 = U2.real
v2 = U2.imag
cp2 = 1 - (u2**2 + v2**2)/u_inf

print(numpy.shape(u))
print(u2[-62, 0])
print(v2[-62, 0])
print(numpy.min(cp[:,0]))
print(numpy.min(cp2[:,0]))

##############################
# plotting the streamlines
x_start, x_end = -3.0, 3.0
y_start, y_end = -3.0, 3.0
size = 5.
pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.quiver(X2, Y2, u2, v2, cp2, scale = 100.)
#pyplot.streamplot(X, Y, u, v, density=2, linewidth=1, arrowsize=2, arrowstyle='->')

#contf = pyplot.contourf(X2, Y2, cp2, levels=numpy.linspace(-0.8, 0.8, 50), extend='both')
#cbar = pyplot.colorbar(contf)
#cbar.set_label('$C_p$', fontsize=16)
#cbar.set_ticks(numpy.linspace(-0.8,0.8,5))


