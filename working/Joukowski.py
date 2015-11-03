import numpy
import math
from matplotlib import pyplot
import pblocks

j = (-1)**0.5

# plotting parameters
x_start, x_end = -5.0, 5.0
y_start, y_end = -5.0, 5.0
size = 5.
pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)

def Jouk_transform(c, z):
    '''Joukowski transformation of z'''
    
    return z + c/z

def rotate(xc, yc, AoA, z):
    ''' Rotation'''
    
    return (z - (xc + j*yc)) * numpy.exp(-j*AoA)

dx = dy = 0.1
c = 1.0
#R = ((c-xc)**2 + yc**2)**0.5
R = 1.15
xc = -0.15
yc = 0.0

N_r = 100
N_th = 145
r_start, r_end = R, 5.0                  # boundaries in the r-direction
theta_start, theta_end = 0.0, 2*math.pi  # boundaries in the y-direction
r = numpy.linspace(r_start, r_end, N_r)    # creates a 1D-array with the x-coordinates
theta = numpy.linspace(theta_start, theta_end, N_th)    # creates a 1D-array with the y-coordinates
R_mesh, Theta = numpy.meshgrid(r, theta) 

X = xc + R_mesh*numpy.cos(Theta)
Y = yc + R_mesh*numpy.sin(Theta)
Z = X + Y*j
#pyplot.scatter(X, Y, s=1)

# Joukowski transformation
Z_c = Jouk_transform(c, Z)
X_c = Z_c.real
Y_c = Z_c.imag
#pyplot.scatter(X_c, Y_c, s=1)

# freestream
u_inf = 1.0
AoA = 20./180 * math.pi

# rotate mesh
Z_dash = rotate(xc, yc, AoA, Z)
X_dash = Z_dash.real
Y_dash = Z_dash.imag
#pyplot.scatter(X_dash, Y_dash, s=1)

psi_stream = pblocks.get_stream_function_freestream(u_inf, 0., X_dash, Y_dash)
u_stream, v_stream = pblocks.get_velocity_freestream(u_inf, 0., X_dash, Y_dash)

# dublet
strength = R**2 * u_inf * 2*math.pi
print('strength = ', strength)
psi_dublet = pblocks.get_stream_function_doublet(strength, 0., 0., X_dash, Y_dash)
u_dublet, v_dublet = pblocks.get_velocity_doublet(strength, 0., 0., X_dash, Y_dash)

# vortex
gamma = 4*math.pi*R*u_inf*math.sin(AoA)
print('gamma = ', gamma)
psi_vortex = pblocks.get_stream_function_vortex(gamma, 0, 0., X_dash, Y_dash)
u_vortex, v_vortex = pblocks.get_velocity_vortex(gamma, 0, 0., X_dash, Y_dash)

psi = psi_stream + psi_dublet + psi_vortex
u_dash = u_stream + u_dublet + u_vortex
v_dash = v_stream + v_dublet + v_vortex

pyplot.plot(X_c[:,0], Y_c[:,0])
pyplot.contour(X_c, Y_c, psi, 50, colors='k', linestyles='solid')

# back transformation of velocity
U = (u_dash - j*v_dash) * numpy.exp(-j*AoA)
u = U.real
v = -U.imag

#U = u + j*v
cp = 1 - (u**2 + v**2)/u_inf
print('cp_max = ', numpy.max(cp[:,0]))

U_c = U / (1 - (c/Z)**2)
u_c = U_c.real
v_c = -U_c.imag
cp_c = 1 - (u_c**2 + v_c**2)/u_inf
print('cp_c_max = ', numpy.max(cp_c[:,0]))

print('indices = ', numpy.where(cp_c[:,0] == 1)[0])
print('u = ', u_c[91, 0])
print('v = ', v_c[91, 0])
print('cp = ', cp_c[110, 0])

#print(theta)
D, L = 0., 0.
for i in range(len(cp[:,0])):
    D +=  0.5*u_inf**2 * cp[i,0] * math.cos(theta[i])
    L += -0.5*u_inf**2 * cp[i,0] * math.sin(theta[i])
print('D = ', D * R*2*math.pi/(N_th-1))
print('L = ', L * R*2*math.pi/(N_th-1))

##############################
# plotting the contours
x_start, x_end = -5.0, 5.0
y_start, y_end = -5.0, 5.0
size = 5.
pyplot.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
#pyplot.quiver(X, Y, u, v, cp, scale = 100.)

contf = pyplot.contourf(X, Y, cp, levels=numpy.linspace(-1, 1, 100), extend='both')
cbar = pyplot.colorbar(contf)
cbar.set_label('$C_p$', fontsize=16)
cbar.set_ticks(numpy.linspace(-0.88,0.88,9))

