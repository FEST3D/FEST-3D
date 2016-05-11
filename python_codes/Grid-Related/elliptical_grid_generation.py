# This module attempts to create a elliptical grid from a given set of
# edges or faces

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
import pylab
from refinegrid import write_grid_3D
import os
import sys
import getopt
from pprint import pprint

def get_grid_layer_with_spherical_bump(R, X_offset, y, z, dy, dz):  
    # Usage: get_grid_layer_with_spherical_bump(R, delta_X, delta_Y, dx, X_end_wall)

    centre = (X_offset, 0.5*(y[0] + y[-1]), 0.5*(z[0] + z[-1]))
    x = np.array([X_offset])
    X, Y, Z = np.meshgrid(x, y, z, indexing = 'ij')

    m, n, o = X.shape
    global sphere_indices
    sphere_indices = []

    for j in range(n):
        for k in range(o):
            r = np.sqrt((Y[0, j, k] - centre[1]) ** 2. + \
                        (Z[0, j, k] - centre[2]) ** 2.)
            if (r ** 2.) <= (R ** 2.):
                X[0, j, k] = X_offset - np.sqrt(R**2. - r**2.)
                sphere_indices.append([0, j, k])

    return (X, Y, Z)


def new_plot(fig, X, Y, Z):
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    ax.view_init(32, -52)
    ax.scatter(X,Y,Z) 
    pylab.ion()
    return ax


def add_plot(fig, X, Y, Z):
    ax.scatter(X, Y, Z)



def plot_axes_vector(x):
    # Gets a stretcheddistribution of the discretized axes vector
    a = x[0]
    b = x[-1]
   
    x_norm = (x-a)/(b-a)

    # See page 2: http://www3.nd.edu/~gtryggva/CFD-Course/2011-Lecture-25.pdf

    # Working cubic function below: Use if needed
#   x_new = a + (0.5*((2*(x_norm - 0.5))**3. + 1)) * (b-a)

    # Stretching strength it seems (<0 for boundary layer like grid 
    # stretching, >0 for more points in the centre)
    # L: no idea what it means. Seems to work. Change with caution!
    # xc: Point about which the slope changes. Here, it is the centre
    A = 3.0
    L = 1.0
    xc = 0.5

#   x_new = a + ((L*x_norm) + A*(xc - L*x_norm)*(1-x_norm)*x_norm ) * (b-a)
    x_new = a + (x_norm ** (1/2.0))*(b - a)

    print x
    print x_new




    y1 = 5 * np.ones(x.shape)
    y2 = 10 * np.ones(x_new.shape)

    plt.plot(x, y1, 'r*')
    plt.plot(x_new, y2, 'b*')
    plt.show()


def get_distributed_axes_vector(x):
    # See page 2: http://www3.nd.edu/~gtryggva/CFD-Course/2011-Lecture-25.pdf
    a = x[0]
    b = x[-1]
    x_norm = (x-a)/(b-a)
    A = 3.0
    L = 1.0
    xc = 0.5
    x_new = a + ((L*x_norm) + A*(xc - L*x_norm)*(1-x_norm)*x_norm ) * (b-a)
    return x_new


def get_one_directional_stretching(x):
    # See page 2: http://www3.nd.edu/~gtryggva/CFD-Course/2011-Lecture-25.pdf
    a = x[0]
    b = x[-1]
    x_norm = (x-a)/(b-a)
    x_new = a + (x_norm ** (1/2.0))*(b - a)
    return x_new
    


def get_plane_grids():
    # Note: assumed that the bump is in the right side - or at the x_end limit
    
    x_limits = [0, 20]
    y_limits = [-20, 20]
    z_limits = [-20, 20]
    R = 4
    dx = 0.5
    dy = 0.6
    dz = 0.6

    plane_grids = {}
    axes_vectors = {}

    x = np.linspace(x_limits[0], x_limits[1], \
                    round((x_limits[1] - x_limits[0]) / dx + 1), endpoint=True)
    y = np.linspace(y_limits[0], y_limits[1], \
                    round((y_limits[1] - y_limits[0]) / dy + 1), endpoint=True)
    z = np.linspace(z_limits[0], z_limits[1], \
                    round((z_limits[1] - z_limits[0]) / dz + 1), endpoint=True)

    x = get_one_directional_stretching(x)
    y = get_distributed_axes_vector(y)
    z = get_distributed_axes_vector(z)

    axes_vectors['x'] = x
    axes_vectors['y'] = y
    axes_vectors['z'] = z

 #  fig = plt.figure()

    # i face
    X, Y, Z = np.meshgrid(np.array([x[0]]), y, z, indexing='ij')
 #  ax = new_plot(fig, X, Y, Z)
    plane_grids['left'] = {'X':X, 'Y':Y, 'Z':Z}

    # i + 1 face
    X, Y, Z = get_grid_layer_with_spherical_bump(R, x[-1], y, z, dy, dz)
 #  ax.scatter(X, Y, Z)
    plane_grids['right'] = {'X':X, 'Y':Y, 'Z':Z}

    # j face
    X, Y, Z = np.meshgrid(x, np.array([y[0]]), z, indexing='ij')
 #  ax.scatter(X, Y, Z)
    plane_grids['front'] = {'X':X, 'Y':Y, 'Z':Z}

    # j + 1 face
    X, Y, Z = np.meshgrid(x, np.array([y[-1]]), z, indexing='ij')
 #  ax.scatter(X, Y, Z)
    plane_grids['back'] = {'X':X, 'Y':Y, 'Z':Z}

    # k face
    X, Y, Z = np.meshgrid(x, y, np.array([z[0]]), indexing='ij')
 #  ax.scatter(X, Y, Z)
    plane_grids['bottom'] = {'X':X, 'Y':Y, 'Z':Z}

    # k + 1 face
    X, Y, Z = np.meshgrid(x, y, np.array([z[-1]]), indexing='ij')
 #  ax.scatter(X, Y, Z)
    plane_grids['top'] = {'X':X, 'Y':Y, 'Z':Z}

 #  plt.show()
    
    return (plane_grids, axes_vectors)
     
    
def get_starting_grid():   
    
    # We are going to solve the laplace equation in the interior of the mesh grid

    plane_grids, axes_vectors = get_plane_grids()

    x = axes_vectors['x']
    y = axes_vectors['y']
    z = axes_vectors['z']

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    m, n, o = X.shape

    # Fixing outer positions
    for i in range(m):
        for j in range(n):
            for k in range(o):
                if X[i, j, k] == x[0]:
                    X[i, j, k] = plane_grids['left']['X'][0, j, k]
                    Y[i, j, k] = plane_grids['left']['Y'][0, j, k]
                    Z[i, j, k] = plane_grids['left']['Z'][0, j, k]
                elif X[i, j, k] == x[-1]:
                    X[i, j, k] = plane_grids['right']['X'][0, j, k]
                    Y[i, j, k] = plane_grids['right']['Y'][0, j, k]
                    Z[i, j, k] = plane_grids['right']['Z'][0, j, k]
                elif Y[i, j, k] == y[0]:
                    X[i, j, k] = plane_grids['front']['X'][i, 0, k]
                    Y[i, j, k] = plane_grids['front']['Y'][i, 0, k]
                    Z[i, j, k] = plane_grids['front']['Z'][i, 0, k]
                elif Y[i, j, k] == y[-1]:
                    X[i, j, k] = plane_grids['back']['X'][i, 0, k]
                    Y[i, j, k] = plane_grids['back']['Y'][i, 0, k]
                    Z[i, j, k] = plane_grids['back']['Z'][i, 0, k]
                elif Z[i, j, k] == z[0]:
                    X[i, j, k] = plane_grids['bottom']['X'][i, j, 0]
                    Y[i, j, k] = plane_grids['bottom']['Y'][i, j, 0]
                    Z[i, j, k] = plane_grids['bottom']['Z'][i, j, 0]
                elif Z[i, j, k] == z[-1]:
                    X[i, j, k] = plane_grids['top']['X'][i, j, 0]
                    Y[i, j, k] = plane_grids['top']['Y'][i, j, 0]
                    Z[i, j, k] = plane_grids['top']['Z'][i, j, 0]

    # Making interior cells interior to the new boundary

    # Since we know that the new boundary is due to the bump:
    # Find the min X in the right boundary which turns out to be the
    # maximum allowable x limit
    new_x_max_right = np.amin(plane_grids['right']['X'])

    # Just contract the mesh along the x direction by this amount
    # This amounts to just changing the X indices!!
    X[1:m-1, 1:n-1, 1:o-1] = x[0] + \
                             ( ((new_x_max_right - x[0])/(x[-1] - x[0])) * \
                             X[1:m-1, 1:n-1, 1:o-1] )

    return (X, Y, Z)


def get_elliptical_grid():
    # Now that the proper meshgrid is obtained, time to iterate towards 
    # elliptical grid.
    X, Y, Z = get_starting_grid()

    global sphere_indices
    print 'Sphere indices: '
    pprint(sphere_indices)

    no_iterations = 70
#   no_iterations = 0
    
    m, n, o = X.shape
    print X.shape
    
    for i in range(len(sphere_indices)):
    # The required indices are imx, j, k
        sphere_indices[i][0] = m
    

    for itr in range(no_iterations):
        # Interior cells only
        print 'Iteration: ', itr
        for i in range(1, m-1):
            for j in range(1, n-1):
                for k in range(1, o-1):
                    # Update is the average of neighbouring cells
                    X[i, j, k] = (X[i-1, j, k] + X[i+1, j, k] + \
                                  X[i, j-1, k] + X[i, j+1, k] + \
                                  X[i, j, k-1] + X[i, j, k+1]) / 6.0
                    Y[i, j, k] = (Y[i-1, j, k] + Y[i+1, j, k] + \
                                  Y[i, j-1, k] + Y[i, j+1, k] + \
                                  Y[i, j, k-1] + Y[i, j, k+1]) / 6.0
                    Z[i, j, k] = (Z[i-1, j, k] + Z[i+1, j, k] + \
                                  Z[i, j-1, k] + Z[i, j+1, k] + \
                                  Z[i, j, k-1] + Z[i, j, k+1]) / 6.0


    return (X, Y, Z)


def parse_grid(filename):
    f = open(filename)
    d = f.read()
    f.close()
    d = d.splitlines()

    x = []
    y = []
    z = []

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
        if len(gridsize) == 3:
            z.append(float(t[2]))

    x = np.array(x)
    y = np.array(y)
    
    if len(z) != 0:
        z = np.array(z)
        grid = {
             'x' : x,
             'y' : y,
             'z' : z
        }
    else:    
        grid = {
            'x': x,
            'y': y,
        }
    
    return (grid, gridsize)


def writevtk_3D(grid, gridsize, data, filename, comment=None):
    # This function writes ONLY the grid file. Used for testing the grid
    '''Writes the grid and data in vtk format - 2D files 
    See www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf for 
    information regarding the format of vtk files.  '''

    # We need to write the point list in the order of changing i, then j
    # and then k. But this is the default order given in the grid file

    # Attemping structured grid

    f = open(filename + '.part', 'w')

    imx, jmx, kmx = gridsize
    num_points = imx * jmx * kmx
    num_cells = (imx - 1) * (jmx - 1) * (kmx - 1)

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
    f.write('DIMENSIONS ' + str(imx) + ' ' + str(jmx) + ' ' + str(kmx))
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
        f.write(str(grid['z'][i]))  # z coordinate
        f.write('\n')
    f.write('\n')
    f.close()

    # Rename the file: remove the .part from the end.
    os.rename(filename + '.part', filename)


def translate_fortran_to_vtk(gridfile, datafile, opfilename, filecomment):
    (grid, gridsize) = parse_grid(gridfile)
    data = None
    if len(gridsize) == 2:
        writevtk(grid, gridsize, data, opfilename, filecomment)
    elif len(gridsize) == 3:
        writevtk_3D(grid, gridsize, data, opfilename, filecomment)


def write_sphere_indices_to_file():
    global sphere_indices
    f = open('sphere-indices.txt', 'w')
    f.write(len(sphere_indices).__str__() + '\n')
    for index in sphere_indices:
        f.write(index[0].__str__() + ' ' + index[1].__str__() + ' ' + index[2].__str__() + '\n') 
    f.close()


def main_function():
    # Could not figure out a better name
    gridfile = 'spherical-bumpgrid.txt'
    datafile = None
    opfilename = 'spherical-bumpgrid.vtk'
    filecomment = 'Testing 3D grid for 3D bow shock'

    print 'Getting the spherical grid'
    X, Y, Z = get_elliptical_grid()

    write_sphere_indices_to_file()
    print 'Write gridfile as text'
    write_grid_3D(gridfile, X, Y, Z)

    print 'Converting to vtk'
    translate_fortran_to_vtk(gridfile, datafile, opfilename, filecomment)
