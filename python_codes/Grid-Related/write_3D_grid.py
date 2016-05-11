from __future__ import division
import numpy as np
from refinegrid import write_grid_3D, read_grid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def convert_to_3D(gridfile):
    imx, jmx, gridx, gridy = read_grid(gridfile)
    print 'Grid file read'
    st = 0
    end = 0.03
    dz = 0.03
    z_coord_list = np.linspace(st, end, (end - st)/dz + 1)
    num = len(z_coord_list)

    gridx_new = np.zeros((imx, jmx, num))
    gridy_new = np.zeros((imx, jmx, num))
    gridz_new = np.zeros((imx, jmx, num))

    for k in range(len(z_coord_list)):
        gridx_new[:, :, k] = gridx
        gridy_new[:, :, k] = gridy
        gridz_new[:, :, k] = z_coord_list[k]

    gridfile_new = gridfile[:-4] + '-3D' + gridfile[-4:]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal')
    ax.scatter(gridx_new.flatten(), gridy_new.flatten(), gridz_new.flatten(), c = 'b', marker = '*') 
    plt.show()

    write_grid_3D(gridfile_new, gridx_new, gridy_new, gridz_new)


if __name__ == '__main__':
    gridfile = 'Rchannel.txt'
    convert_to_3D(gridfile)
