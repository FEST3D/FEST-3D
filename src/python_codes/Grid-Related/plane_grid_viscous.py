# This module gets a refined bumpgrid near the wall for viscous
# flow simulation

import numpy as np
import matplotlib.pyplot as plt
import pylab
from refinegrid import write_grid
import os
import sys
import getopt
from pprint import pprint


def get_distributed_axes_vector(x):
    # See page 2: http://www3.nd.edu/~gtryggva/CFD-Course/2011-Lecture-25.pdf
    a = x[0]
    b = x[-1]
    x_norm = (x-a)/(b-a)
    A = 2.5
    L = 1.0
    xc = 0.5
    x_new = a + ((L*x_norm) + A*(xc - L*x_norm)*(1-x_norm)*x_norm ) * (b-a)
    return x_new


def get_one_directional_stretching(x):
    # See page 2: http://www3.nd.edu/~gtryggva/CFD-Course/2011-Lecture-25.pdf
    b = x[0]
    a = x[-1]
    n = 3.0
    x_norm = (x-a)/(b-a)
    x_new = a + (x_norm ** (1/n))*(b - a)
    return x_new
    

def one_directional_stretching_matrix(Y):
    # Assume stretching along Y, along column

    a = Y[0, :]
    b = Y[-1, :]
    n = 2.5

    p,q = Y.shape
    for i in range(p):
        y = Y[i, :]
        a = y[-1]
        b = y[0]
        y = np.linspace(a, b, len(y))
        y_norm = (y-a)/(b-a)
        y_new = a + (y_norm ** (1/n))*(b - a)
        Y[i, :] = y_new
    
    return Y


def two_directional_stretching_matrix(Y):
    # Assume stretching along Y, along column
    # Do not use lower value than -2 for method 1
    A = -2.0
    L = 1.0
    yc = 0.5

    # Method 2: www.cfd-online/Wiki/Structured_mesh_generation
    delta = 3.5

    p,q = Y.shape
    for i in range(p):
        y = Y[i, :]
        a = y[0]
        b = y[-1]
        y = np.linspace(a, b, len(y))
        y_norm = (y-a)/(b-a)

        # Method 1
    #   y_new = a + ((L*y_norm) + A*(yc - L*y_norm)*(1-y_norm)*y_norm ) * (b-a)

        # Method 2
        y_new = a + (0.5*(1 + (np.tanh(delta*(y_norm-0.5))/np.tanh(delta/2.0))))*(b-a)
        Y[i, :] = y_new

    #   print y_new
    #   print y_norm

    #   print a, b, y_new[0], y_new[-1]
    #   plt.plot(np.ones(y_new.shape), y_new, '*')
    #   plt.show()
    
    return Y


def get_line_with_bump(R, x_center, x, y_coord, sign):
    # We are assuming that line vector is along x
    # Sign = -1 or +1 depending on which semi circle you want

    # We are using x and y specifically, but it can be generalised
    # The output will be such that 'x' will be modified. If there
    # should be a bump along the 'y' coordinate, then the receiving
    # function should swap it

    y = np.ones(x.shape) * y_coord
    imx,  = x.shape
    for i in range(imx):
        if abs(x[i] - x_center) <= R:
            y[i] = y_coord + sign*np.sqrt(R**2 - (x[i] - x_center)**2)
    
    return (x, y)

    
def get_plane_grids(bump='bottom'):
    x_limits = [0, 0.45]
    y_limits = [0, 0.075]
    R = 2
    dx = 0.003
    dy = 0.001
    plane_grids = {}
    axes_vectors = {}

    x = np.linspace(x_limits[0], x_limits[1], \
                    round((x_limits[1] - x_limits[0]) / dx + 1), endpoint=True)
    y = np.linspace(y_limits[0], y_limits[1], \
                    round((y_limits[1] - y_limits[0]) / dy + 1), endpoint=True)

    if bump == 'bottom' or bump == 'top':
        pass
      # x = get_distributed_axes_vector(x)
        #y = get_one_directional_stretching(y)
    elif bump == 'right' or bump == 'left':
        y = get_distributed_axes_vector(y)
        x = get_one_directional_stretching(x)
        
    axes_vectors['x'] = x
    axes_vectors['y'] = y    

    plane_grids['left'] = {'X':np.ones(y.shape)*x[0], 'Y':y}

    plane_grids['right'] = {'X':np.ones(y.shape)*x[-1], 'Y':y}

    plane_grids['top'] = {'X':x, 'Y':np.ones(x.shape)*y[-1]}

    plane_grids['bottom'] = {'X':x, 'Y':np.ones(x.shape)*y[0]}

  # if bump == 'bottom':
  #     X, Y = get_line_with_bump(R, 0.0, x, y[0], 1)
  #     plane_grids['bottom'] = {'X':X, 'Y':Y}
  # elif bump == 'top':
  #     X, Y = get_line_with_bump(R, 0.0, x, y[-1], -1)
  #     plane_grids['bottom'] = {'X':X, 'Y':Y}
  # elif bump == 'left':
  #     Y, X = get_line_with_bump(R, 0.0, y, x[0], 1)
  #     plane_grids['left'] = {'X':X, 'Y':Y}
  # elif bump == 'right':
  #     Y, X = get_line_with_bump(R, 0.0, y, x[-1], -1)
  #     plane_grids['right'] = {'X':X, 'Y':Y}
    
    return (plane_grids, axes_vectors)


def get_starting_grid(bump='bottom'):   
    
    # We are going to solve the laplace equation in the interior of the mesh grid

    plane_grids, axes_vectors = get_plane_grids(bump)

    x = axes_vectors['x']
    y = axes_vectors['y']

    X, Y = np.meshgrid(x, y, indexing='ij')
    m,n = X.shape

    # Fixing outer positions
    for i in range(m):
        for j in range(n):
            if X[i, j] == x[0]:
                X[i, j] = plane_grids['left']['X'][j]
                Y[i, j] = plane_grids['left']['Y'][j]
            elif X[i, j] == x[-1]:
                X[i, j] = plane_grids['right']['X'][j]
                Y[i, j] = plane_grids['right']['Y'][j]
            elif Y[i, j] == y[0]:
                X[i, j] = plane_grids['bottom']['X'][i]
                Y[i, j] = plane_grids['bottom']['Y'][i]
            elif Y[i, j] == y[-1]:
                X[i, j] = plane_grids['top']['X'][i]
                Y[i, j] = plane_grids['top']['Y'][i]

    Y = two_directional_stretching_matrix(Y)

#   if bump == 'bottom':
#       y_lim = np.amax(plane_grids[bump]['Y'])
#       print y_lim
#       Y[1:m-1, 1:n-1] = y_lim + \
#                            ( ((y_lim - y[0])/(y[-1] - y[0])) * \
#                            Y[1:m-1, 1:n-1] )

#        Y[m-2:0:-1, n-2:0:-1] = y[-1] + \
#                             ( ((y_lim - y[0])/(y[-1] - y[0])) * \
#                             Y[m-2:0:-1, n-2:0:-1] )
#   elif bump == 'top':
#       y_lim = np.amin(plane_grids[bump]['Y'])
#       Y[1:m-1, 1:n-1] = y[0] + \
#                            ( ((y_lim - y[0])/(y[-1] - y[0])) * \
#                            Y[1:m-1, 1:n-1] )
#   elif bump == 'right':
#       x_lim = np.amin(plane_grids[bump]['X'])
#       X[1:m-1, 1:n-1] = x[0] + \
#                            ( ((x_lim - x[0])/(x[-1] - x[0])) * \
#                            X[1:m-1, 1:n-1] )
#   elif bump == 'left':
#       x_lim = np.amax(plane_grids[bump]['X'])
#       X[m-2:0:-1, n-2:0:-1] = x[-1] + \
#                            ( ((x_lim - x[0])/(x[-1] - x[0])) * \
#                            X[m-2:0:-1, n-2:0:-1] )

    return (X, Y)


def get_elliptical_grid(bump):
    # Now that the proper meshgrid is obtained, time to iterate towards 
    # elliptical grid.
    X, Y = get_starting_grid(bump)
    m, n = X.shape

    no_iterations = 70
    
    for itr in range(no_iterations):
        # Interior cells only
        print 'Iteration: ', itr
        for i in range(1, m-1):
            for j in range(1, n-1):
                X[i, j] = (X[i-1, j] + X[i+1, j] + \
                           X[i, j-1] + X[i, j+1]) / 4.0
                Y[i, j] = (Y[i-1, j] + Y[i+1, j] + \
                           Y[i, j-1] + Y[i, j+1]) / 4.0

    return (X, Y)


def parse_grid(filename):
    f = open(filename)
    d = f.read()
    f.close()
    d = d.splitlines()

    x = []
    y = []

    gridsize = d.pop(0)  # The first line of the file should contain imx, jmx and kmx
    gridsize = gridsize.strip()
    gridsize = gridsize.split()
    gridsize = [int(i) for i in gridsize]

    for point in d:
        p = point
        p = p.strip()
        t = p.split()
        x.append(float(t[0]))
        y.append(float(t[1]))

    x = np.array(x)
    y = np.array(y)
    
    grid = {
        'x': x,
        'y': y,
    }
    
    return (grid, gridsize)


def writevtk_2D(grid, gridsize, data, filename, comment=None):
    # This function writes ONLY the grid file. Used for testing the grid
    '''Writes the grid and data in vtk format - 2D files 
    See www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf for 
    information regarding the format of vtk files.  '''

    # We need to write the point list in the order of changing i, then j
    # and then k. But this is the default order given in the grid file

    # Attemping structured grid

    f = open(filename + '.part', 'w')

    imx, jmx = gridsize
    num_points = imx * jmx
    num_cells = (imx - 1) * (jmx - 1)

    # Write Header
    f.write('# vtk DataFile Version 3.1\n')
    if comment:
        f.write(comment)
    else:
        f.write('Glomar response.')
    f.write('\n')
    f.write('ASCII\n')
    f.write('DATASET STRUCTURED_GRID\n')
    f.write('\n')
    f.write('DIMENSIONS ' + str(imx) + ' ' + str(jmx) + ' 1')
    f.write('\n')

    # Write Pointdata
    f.write('POINTS ')
    f.write(str(num_points))
    f.write(' FLOAT\n')
    for i in range(num_points):
        f.write(str(grid['x'][i]))
        f.write(' ')
        f.write(str(grid['y'][i]))
        f.write(' ')
        f.write('0') # z coordinate
        f.write('\n')
    f.write('\n')
    f.close()

    # Rename the file: remove the .part from the end.
    os.rename(filename + '.part', filename)


def translate_fortran_to_vtk(gridfile, datafile, opfilename, filecomment):
    (grid, gridsize) = parse_grid(gridfile)
    data = None
    writevtk_2D(grid, gridsize, data, opfilename, filecomment)


def main_function():
    # Could not figure out a better name
    gridfile = 'plane-viscous.txt'
    datafile = None
    opfilename = 'plane-viscous.vtk'
    filecomment = 'Testing 2D grid for viscous flows'
    bump = 'bottom'

    print 'Getting the spherical grid'
#    X, Y = get_elliptical_grid(bump)
    X, Y = get_starting_grid(bump)

#    plt.scatter(X, Y)
#    plt.show()

    print 'Write gridfile as text'
    write_grid(gridfile, X, Y)

    print 'Converting to vtk'
    translate_fortran_to_vtk(gridfile, datafile, opfilename, filecomment)


if __name__ == '__main__':
    main_function()
