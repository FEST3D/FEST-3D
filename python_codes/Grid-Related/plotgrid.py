from __future__ import division
import pdb
import numpy as np
import matplotlib.pyplot as plt
import sys

# Usage: python plotgrid <gridfilename>

# This script plots the grid file
# TODO: This code is to be tested only for 2D uniform grid.

def read_grid(gridfile):
    f1 = open(gridfile, 'r')
    [imx, jmx] = [int(i) for i in f1.readline().split()]
    gridx = np.zeros((imx,jmx))
    gridy = np.zeros((imx,jmx))
    for j in range(jmx):
        for i in range(imx):
            [gridx[i][j], gridy[i][j]] = [float(k) for k in f1.readline().split()]
    f1.close()
    return [imx, jmx, gridx, gridy]


def plot_grid(gridfilename):
    imx, jmx, gridx, gridy = read_grid(gridfilename)
    print 'Grid file read'

    print 'Plotting grid'
    plt.plot(gridx.flatten(), gridy.flatten(), 'b*')
    plt.show()


if __name__ == '__main__':
    arg_len = len(sys.argv)
    if arg_len == 1:
        print 'Not enough arguments'
        print 'Usage: python plotgrid <gridfilename>'
        sys.exit(1)

    try:
        gridfilename = sys.argv[1]
    except IndexError:
        print 'grid file necessary'
        print 'Usage: python plotgrid <gridfilename>'
        sys.exit(1)
    
    plot_grid(gridfilename)
