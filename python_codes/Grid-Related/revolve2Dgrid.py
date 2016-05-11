from __future__ import division
import numpy as np
import refinegrid as rg

def rotate_2D_grid(gridx, gridy, dtheta, kmx):
    # This function revolves the grid about the x axis
    # The revolution is -dtheta to +dtheta

    imx, jmx = gridx.shape
    
    deltaY = abs(gridy[0, jmx-1] - gridy[0, 0])
    R = 0.01 * deltaY
    R = 0.0

    GridX = np.zeros((imx, jmx, kmx))
    GridY = np.zeros((imx, jmx, kmx))
    GridZ = np.zeros((imx, jmx, kmx))

    for i in range(imx):
        for j in range(jmx):
            X = gridx[i, j]
            Y = gridy[i, j]
            for k in range(kmx):
                # k is from 0 to kmx - 1 because python...
                theta = dtheta * ((2*k/(kmx - 1)) - 1)
                # X remains the same
                GridX[i, j, k] = X
                GridY[i, j, k] = ((Y + R) * np.cos(theta))
                GridZ[i, j, k] = ((Y + R) * np.sin(theta))

    return (GridX, GridY, GridZ)


if __name__ == '__main__':
    gridfilename2D = 'trap.txt'
    gridfilename3D = 'trap-3D-revolved.txt'

    dtheta = np.pi / 180 * 180.0
    kmx = 30
    
    imx, jmx, gridx, gridy = rg.read_grid(gridfilename2D)

    GridX, GridY, GridZ = rotate_2D_grid(gridx, gridy, dtheta, kmx)
    rg.save_to_vtk_3D(GridX, GridY, GridZ, gridfilename3D[:-4])
    rg.write_grid_3D(gridfilename3D, GridX, GridY, GridZ)



