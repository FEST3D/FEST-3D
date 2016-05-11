# This code reads the 2D uniform grid supplied as bumpgrid and rotates
# it by required degrees. It uses the plotgrid module

from __future__ import division
from refinegrid import read_grid, write_grid
import numpy as np
import matplotlib.pyplot as plt
import sys

gridfile = 'bumpgrid.txt'

def reverse_row(gridx, gridy):
    [imx, jmx] = gridx.shape
    gridx_new = np.zeros((imx, jmx))
    gridy_new = np.zeros((imx, jmx))
    
    for i in range(imx):
        gridx_new[i, :] = gridx[imx-i-1, :]
        gridy_new[i, :] = gridy[imx-i-1, :]
    return (gridx_new, gridy_new)    


def reverse_column(gridx, gridy):
    [imx, jmx] = gridx.shape
    gridx_new = np.zeros((imx, jmx))
    gridy_new = np.zeros((imx, jmx))
    
    for j in range(jmx):
        gridx_new[:, j] = gridx[:, jmx-j-1]
        gridy_new[:, j] = gridy[:, jmx-j-1]
    return (gridx_new, gridy_new)    


def rotate_90(gridx, gridy):
    [imx, jmx] = gridx.shape
    gridx_new = np.zeros((imx, jmx))
    gridy_new = np.zeros((imx, jmx))

    # First transpose
    gridx_new = np.transpose(gridx)
    gridy_new = np.transpose(gridy)

    # Then Reverse each row
    [gridx_new, gridy_new] = reverse_row(gridx_new, gridy_new)
    return (gridx_new, gridy_new)    


def rotate_90cw(gridx, gridy):
    [imx, jmx] = gridx.shape
    gridx_new = np.zeros((imx, jmx))
    gridy_new = np.zeros((imx, jmx))

    # First transpose
    gridx_new = np.transpose(gridx)
    gridy_new = np.transpose(gridy)

    # Then Reverse each column
    [gridx_new, gridy_new] = reverse_column(gridx_new, gridy_new)
    return (gridx_new, gridy_new)    


def rotate_180(gridx, gridy):
    [imx, jmx] = gridx.shape
    gridx_new = np.zeros((imx, jmx))
    gridy_new = np.zeros((imx, jmx))
    
    print 'here'

    # Rotate by 90 twice
    [gridx_new, gridy_new] = rotate_90(gridx, gridy)
    [gridx_new, gridy_new] = rotate_90(gridx_new, gridy_new)
     
    # Stackoverflow suggests the following two lines but code works
    # ONLY without them
    # Reverse each row
    # [gridx_new, gridy_new] = reverse_row(gridx_new, gridy_new)
    # Reverse each column 
    # [gridx_new, gridy_new] = reverse_column(gridx_new, gridy_new)
    return (gridx_new, gridy_new)    


def rotate_the_array(gridx, gridy, theta):
    # This functions rotates the array as per the stackoverflow answer by
    # "dimple" in the following page: 
    # stackoverflow.com/questions/42519/how-do-you-rotate-a-two-dimensional-array
    #
    # A more graphical reasoning has been given in another answer by "tweaking"
    #
    # The basic idea is this:
    # The rotation of arrays is known for the four basic angles: +-90, +180
    # -180 is same as +180
    # Based on the input angle, the closest one of these four will be chosen

    [imx, jmx] = gridx.shape
    gridx_new = np.zeros((imx, jmx))
    gridy_new = np.zeros((imx, jmx))
 
    if abs(theta) > 360 or theta < 0:
        theta = theta%360

    if (45 < theta <= 135):
        [gridx_new, gridy_new] = rotate_90(gridx, gridy)
    elif (135 < theta <= 225):
        [gridx_new, gridy_new] = rotate_180(gridx, gridy)
    elif (225 < theta <= 315):
        [gridx_new, gridy_new] = rotate_90cw(gridx, gridy)
    elif (0 <= theta <= 45) or (315 < theta <= 360):
        # Nothing to do. Just leave it as it is
        gridx_new = gridx
        gridy_new = gridy

    return (gridx_new, gridy_new)


def transform_grid(gridfile, theta_deg):
    # Note: theta in degrees
    imx, jmx, gridx, gridy = read_grid(gridfile)

    # Transforming into radians for use in np.cos and np.sin
    theta = np.pi/180*theta_deg
    
    gridx_new = np.zeros((imx, jmx))
    gridy_new = np.zeros((imx, jmx))

    # cos and sin of -theta is taken. To find out why, draw the case
    # for theta = 90 and it will be clear. The other terms are from 
    # the usual formulation for rotation of axes. The new axes is 
    # written in terms of the old one. Another way to think about
    # this: rotating the axes has the following formula with theta
    # but rotating the body with the axes same will be equal and
    # opposite and hence -theta

    # This step just transforms the x and y co-ordinates of the grid-
    # points, NOT the array indices, which we shall do subsequently
    gridx_new = gridx*np.cos(-theta) + gridy*np.sin(-theta)
    gridy_new = -gridx*np.sin(-theta) + gridy*np.cos(-theta)

    # Rotate the array, i.e., rotate the index numbering
    gridx_new, gridy_new = rotate_the_array(gridx_new, gridy_new, theta_deg) 

    print 'Writing to file'
    if not theta_deg.is_integer():
        theta_deg = int(theta_deg)
    newgridfilename = gridfile[:-4] + '-' + theta_deg.__str__() + \
                     '-rot' + gridfile[-4:]
    write_grid(newgridfilename, gridx_new, gridy_new)
    print 'Written to file'

    plt.plot(gridx_new.flatten(), gridy_new.flatten(), 'b*')
    plt.show()
    

if __name__ == '__main__':
    # Usage: python transformgrid <grid_file_name> angle
    arg_len = len(sys.argv)
    global theta_deg, gridfile
    
    if arg_len == 1:
        print "Not enough arguments"
        print "Usage: python transformgrid <grid_file_name> angle (in degrees)"
        sys.exit(1)
    
    try:
        gridfile = sys.argv[1]
    except IndexError:
        print "Grid file necessary"
        print "Usage: python transformgrid <grid_file_name> angle (in degrees)"
        sys.exit(1)

    try:
        theta_deg = float(sys.argv[2])
    except IndexError:
        print "Using default value of 0: No rotation"
        theta = 0.0
        
    transform_grid(gridfile, theta_deg)
