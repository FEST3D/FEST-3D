import matplotlib.pyplot as plt
from refinegrid import *

# This program reads the required grid file and plots it along with
# the four corners of the grid (in terms of indices) in the following 
# form:
# 0,0 index: +
# imx, 0 index: .
# 0, jmx index: o
# imx, jmx index: *
#
# The position of these symbols can help debug if the array has 
# rotated properly

gridfile = 'bumpgrid_transformed_modified.txt'
st_index = 32
end_index = 63 + 1

imx, jmx, gridx, gridy = read_grid(gridfile)

print imx, jmx

plt.plot(gridx.flatten(), gridy.flatten(), 'y.')
plt.plot(gridx[imx-1, st_index], gridy[imx-1, st_index], 'g*')
plt.plot(gridx[imx-1, end_index], gridy[imx-1, end_index], 'g*')
plt.plot(gridx[imx-1, jmx-1], gridy[imx-1, jmx-1], 'r*')
plt.plot(gridx[imx-1, 0], gridy[imx-1, 0], 'k.')
plt.plot(gridx[0, 0], gridy[0, 0], 'k+')
plt.plot(gridx[0, jmx-1], gridy[0, jmx-1], 'ko')
plt.plot(gridx[imx-1, 32:63], gridy[imx-1, 32:63], 'k*')

plt.show()
