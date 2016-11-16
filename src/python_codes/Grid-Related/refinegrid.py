# This script refines the row by n2 times and columns by n1
# times.

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# The current format of the grid file supported:
# All lines contain 2 numbers separated by a space in between
# First line contains 2 integers imx and jmx which represent the number
# of rows and columns.
# The other imx*jmx lines give the coordinates of the grid points as
# 'x y'.

# TODO: This code is tested only for 2D and 3D uniform grid.

def read_grid(gridfile):
    f1 = open(gridfile, 'r')
    [imx, jmx] = [int(i) for i in f1.readline().split()]
    gridx = 10*np.zeros((imx,jmx))
    gridy = 10*np.zeros((imx,jmx))
    for j in range(jmx):
        for i in range(imx):
            [gridx[i][j], gridy[i][j]] = [float(k) for k in f1.readline().split()]
    f1.close()
    return [imx, jmx, gridx, gridy]


def read_grid_3D(gridfile):
    f1 = open(gridfile, 'r')

    [imx, jmx, kmx] = [int(i) for i in f1.readline().split()]
    gridx = 10*np.zeros((imx,jmx, kmx))
    gridy = 10*np.zeros((imx,jmx, kmx))
    gridz = 10*np.zeros((imx,jmx, kmx))
    
    print imx, jmx, kmx

    for k in range(kmx):
        for j in range(jmx):
            for i in range(imx):
                [gridx[i][j][k], gridy[i][j][k], gridz[i][j][k]] = [float(l) for l in f1.readline().split()]
    f1.close()
    return [imx, jmx, kmx, gridx, gridy, gridz]


def write_grid(gridfile, gridx, gridy):
    imx, jmx = gridx.shape
    f2 = open(gridfile, 'w')
    f2.write(str(imx) + ' ' + str(jmx) + '\n')
    for j in range(jmx):
        for i in range(imx):
            f2.write(str(gridx[i][j]) + ' ' + str(gridy[i][j]) + '\n')
    f2.close()


def write_grid_3D(gridfile, gridx, gridy, gridz):
    imx, jmx, kmx = gridx.shape
    f2 = open(gridfile, 'w')
    f2.write(str(imx) + ' ' + str(jmx) + ' ' + str(kmx)+ '\n')
    for k in range(kmx):
        for j in range(jmx):
            for i in range(imx):
                f2.write('%.17f'%(gridx[i][j][k]) + ' ' + '%.17f'%(gridy[i][j][k]) + ' ' + '%.17f'%(gridz[i][j][k])+ '\n')
    f2.close()


def save_to_vtk_3D(gridx, gridy, gridz, filename, comment='Grid Generation'):
    if filename[-4:] != '.vtk':
        filename = os.path.join(os.curdir, filename + '.vtk')
    else:
        filename = os.path.join(os.curdir, filename)
    f = open(filename, 'w')

    imx, jmx, kmx = gridx.shape

    f.write(r'# vtk DataFile Version 3.1')
    f.write('\n')
    f.write(comment + '\n')
    f.write('ASCII\n')
    f.write('DATASET STRUCTURED_GRID\n\n')

    f.write('DIMENSIONS ' + imx.__str__() + ' ' + \
                            jmx.__str__() + ' ' + \
                            kmx.__str__() + '\n')
    f.write('POINTS ' + (imx * jmx * kmx).__str__() + ' FLOAT\n') 

    for k in range(kmx):
        for j in range(jmx):
            for i in range(imx):
                f.write(gridx[i][j][k].__str__() + ' ' + \
                        gridy[i][j][k].__str__() + ' ' + \
                        gridz[i][j][k].__str__() + '\n')
    f.write('\n')        
    f.close()


def interpolate_grid(gridx, gridy):
    imx, jmx = gridx.shape
    gridx_new = np.zeros((n1*(imx-1) + 1, n2*(jmx-1) + 1))
    gridy_new = np.zeros((n1*(imx-1) + 1, n2*(jmx-1) + 1))
    # For all the points except top and right row
    for j in range(jmx-1):
        for i in range(imx-1):
            left_vec_x = np.linspace(gridx[i][j], gridx[i][j+1],\
                         n2+1, endpoint=True)
            right_vec_x = np.linspace(gridx[i+1][j], gridx[i+1][j+1],\
                         n2+1, endpoint=True)
            left_vec_y = np.linspace(gridy[i][j], gridy[i][j+1],\
                         n2+1, endpoint=True)
            right_vec_y = np.linspace(gridy[i+1][j], gridy[i+1][j+1],\
                          n2+1, endpoint=True)
            for v in range(n2):
                gridx_new[n1*i:n1*(i+1), n2*j + v] = \
                  np.linspace(left_vec_x[v], right_vec_x[v],\
                    n1, endpoint=False)
                gridy_new[n1*i:n1*(i+1), n2*j + v] = \
                  np.linspace(left_vec_y[v], right_vec_y[v],\
                    n1, endpoint=False)
    
    # Top row
    for i in range(imx-1):
        gridx_new[n1*i:n1*(i+1), n2*(jmx - 1)] = \
          np.linspace(gridx[i][jmx-1], gridx[i+1][jmx-1], \
            n1, endpoint = False)
        gridy_new[n1*i:n1*(i+1), n2*(jmx - 1)] = \
          np.linspace(gridy[i][jmx-1], gridy[i+1][jmx-1], \
            n1, endpoint = False)

    # Right row
    for j in range(jmx-1):
        gridx_new[n1*(imx - 1), n2*j:n2*(j+1)] = \
          np.linspace(gridx[imx-1][j], gridx[imx-1][j+1], \
            n2, endpoint = False)
        gridy_new[n1*(imx - 1), n2*j:n2*(j+1)] = \
          np.linspace(gridy[imx-1][j], gridy[imx-1][j+1], \
            n2, endpoint = False)

    # Top right point
    gridx_new[n1*(imx-1), n2*(jmx-1)] = gridx[imx-1, jmx-1]
    gridy_new[n1*(imx-1), n2*(jmx-1)] = gridy[imx-1, jmx-1]
    
    global newgridfilename
    newgridfilename = gridfile[:-4] + '-' + n1.__str__() + '-' + \
                      n2.__str__() + '-refined' + gridfile[-4:]
    
    return [gridx_new, gridy_new]

if __name__ == '__main__':
    # Usage: python refinegrid.py <grid_file_name> <n1> <n2>
    arg_len = len(sys.argv)
    global n1, n2
    if arg_len == 1:
        print "Not enough arguments"
        print "Usage: python refinegrid.py <grid_file_name> <n1> <n2>"
        sys.exit(1)
    try:
        gridfile = sys.argv[1]
    except IndexError: 
        print "grid file necessary"        
        print "Usage: python refinegrid.py <grid_file_name> <n1> <n2>"
        sys.exit(1)

    try:
        n1 = int(sys.argv[2])
    except IndexError:
        print "Using default value of 1"
        n1 = 1

    try:
        n2 = int(sys.argv[3])
    except IndexError:
        print "Using default value of 1"
        n2 = 1
    
    imx, jmx, gridx, gridy = read_grid(gridfile)
    print 'Grid file read'

    gridx_new, gridy_new = interpolate_grid(gridx, gridy)
    
    plt.plot(gridx.flatten(), gridy.flatten(), '*')
    print 'Got new interpolated grid coordinates'
    print 'Plotting new grid over old grid'
    plt.plot(gridx_new.flatten(), gridy_new.flatten(), '+')
    plt.show()

    print 'Writing new grid file...'
    write_grid(newgridfilename, gridx_new, gridy_new)

    print 'All done! Enjoy!'
